#ifndef EVOLUTION_H // include guard
#define EVOLUTION_H

#include "ndsolve.hpp"
#include "csolver.hpp"
#include "macros.hpp"
#include "logger.hpp"
#include <chrono>


namespace physics{
    class UEvolver{
        public:

            long perft_JJd = 0;

            std :: vector<Logger *> loggers; 

            CSolver csolver;

            UEvolver();

            UEvolver(CSolver _csolver);
            /**
             * @brief Evolves the universe forward in time
             * 
             * @param universe The universe to step forward in time.
             * @param dt The small timestep to step the universe forward.
             */
            void evolve_step(Universe &u, double dt);

            EVector get_forces(EVector &q, EVector &v,Universe &u);

            EVector get_a(EVector &q, EVector &v,Universe &u);

        private:
            //TODO : Put into class
            /**
             * @brief Adds to the 3d force vector to the 3 * particle count force vector
             * 
             * @param forces_super The forces to add to
             * @param force_to_add The forces to add
             * @param n The id of the particle
             */
            void add_to_superpos(EVector &forces_super, EVectorNd &force_to_add, int n);

            /**
             * @brief Sums up interactions for a particle adds to superposition of forces.
             * 
             * @param forces The forces acting on the particle
             * @param p1 The particle to get forces on
             * @param n The particle id
             * @param forces_super The vector to add the forces to
             */
            void superposition_force_p(vForces &forces, Particle &p1,int n, EVector &forces_super);


            /**
             * @brief Sums up interactions for a particle adds to superposition of forces.
             * 
             * @param u The universe 
             * @param p The particle to sum forces for
             * @param n The id of the particle
             * @param forces_super The vector to add the interaction forces to
             *
             */
            void superposition_interaction_p(EVector &q, EVector &v, EVector &m,
                Particle &p, v2Interactions const &interactions_all, EVector &forces_super);

            /**
             * @brief Gets a superposition of all non constraint forces and interactions
             * 
             * @param universe 
             * @return  
             */
            EVector superposition_force(EVector &q, EVector &v, EVector &m, vForces &forces,
               v2Interactions &interactions, int p_count);




            /**
             * @brief Expand m vector to DIMENSION * particle count. Returned m vector will be (m1, m1, m1, ... , mn, mn, mn)
             * 
             * @param universe
             * @return Eigen 
             */
            EVector expand_m_vector(Universe &u);

            void log(Universe &u);
    };
}



#endif
