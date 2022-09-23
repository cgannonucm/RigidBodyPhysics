#ifndef EVOLUTION_H // include guard
#define EVOLUTION_H

#include "ndsolve.h"

namespace physics{
    namespace evolution{
        //TODO : Put into class
        /**
         * @brief Adds to the 3d force vector to the 3 * particle count force vector
         * 
         * @param forces_super The forces to add to
         * @param force_to_add The forces to add
         * @param n The id of the particle
         */
        void add_to_superpos(Eigen :: VectorXd &forces_super, Eigen :: Vector3d const force_to_add, int n){
            for (int j = 0; j < Universe :: DIMENSION; j++){
                auto pos = Universe :: get_pos(n,j);
                forces_super(Universe :: get_pos(n,j)) += force_to_add[j];
            }
        }

        /**
         * @brief Sums up interactions for a particle adds to superposition of forces.
         * 
         * @param forces The forces acting on the particle
         * @param p1 The particle to get forces on
         * @param n The particle id
         * @param forces_super The vector to add the forces to
         */
        void superposition_force_p(std :: vector<Force *> const &forces, Particle const &p1, int n, Eigen :: VectorXd &forces_super){
            //Loop through forces and each interaction force
            for (Force *f : forces){
                //Evaluate force and add it's components to our force vector
                auto f_val = f->get_force(p1);
                add_to_superpos(forces_super, f_val, n);
            }
        }


        /**
         * @brief Sums up interactions for a particle adds to superposition of forces.
         * 
         * @param u The universe 
         * @param p The particle to sum forces for
         * @param n The id of the particle
         * @param forces_super The vector to add the interaction forces to
         *
         */
        void superposition_interaction_p(Universe &u, Particle const &p, int n, Eigen :: VectorXd &forces_super)
        {
            std :: vector<Interaction *> interactions = u.get_interactions(n);

            for (Interaction * interaction : interactions){
                Particle p2;
                //Get 2nd particle in interaction
                if (interaction->pid1 == n){
                    p2 = u.get_p(interaction->pid2);
                }
                else{
                    p2 = u.get_p(interaction->pid1);
                }
                //Get force calculates force on first particle
                //p is always first particle since we are calculating the force on it
                auto f_val = interaction->get_force(p,p2);
                add_to_superpos(forces_super,f_val,n);
            }
        }

        /**
         * @brief Gets a superposition of all non constraint forces and interactions
         * 
         * @param universe 
         * @return Vector3d<long double> 
         */
        Eigen :: VectorXd superposition_force(Universe &universe){
            using namespace Eigen;

            auto p_count = universe.get_p_count();
            auto forces = universe.get_forces();

            VectorXd forces_super = VectorXd :: Zero(Universe :: DIMENSION * p_count);
            
            //Loop through all particles in universe evaluate interactions and forces for all particles
            //And add to superposition-
            for (int n = 0; n < p_count; n++){
                Particle p1 = universe.get_p(n);

                //Add non interaction forces to the superposition of forces
                superposition_force_p(forces,p1,n,forces_super);


                //TODO future optimization by newtons law we only need to calculate interaction
                //force for one particle in the interaction and we know the force on the other
                //Add interaction forces to the superposition of forces
                superposition_interaction_p(universe,p1,n,forces_super);
            }

            return forces_super;
            
        };

        

        
        /**
         * @brief Expand m vector to DIMENSION * particle count. Returned m vector will be (m1, m1, m1, ... , mn, mn, mn)
         * 
         * @param universe
         * @return Eigen 
         */
        Eigen :: VectorXd expand_m_vector(Universe &u){
            using namespace Eigen;
            
            auto p_count = u.get_p_count();

            int dim = Universe :: DIMENSION * p_count;

            auto m_small = u.get_m();

            VectorXd m_large = VectorXd :: Zero(dim);
            
            for (int n = 0; n < p_count; n++){
                for (int j = 0; j < Universe :: DIMENSION;  j++){
                    m_large(Universe :: get_pos(n,j)) = m_small(n);
                }
            }

            return m_large;
            
        }


        /**
         * @brief Steps the particles forward in time based on forces.
         * 
         * @param forces The forces acting on each of the particles.
         * @param universe Universe to step forward in time.
         * @param dt The timestep to step the universe forward - smaller is better.
         */
        void evolve_step_forces(Eigen :: VectorXd &forces, Universe &universe, long double dt){
            using namespace Eigen;

            auto dim = forces.size();

            auto q0 = universe.get_q();
            auto v0 = universe.get_v();
            auto m = expand_m_vector(universe);

            VectorXd q = VectorXd :: Zero(dim);
            VectorXd v = VectorXd :: Zero(dim);

            for (int i = 0; i < dim; i++){

                auto q0_i = q0(i);
                auto v0_i = v0(i);
                auto m_i = m(i);
                auto f_i = forces(i);

                //From newton's law m * x''_i = f(x,x',t)
                //We need to use ndsolve to solve equation x''_i = 1/m * f(x_i,x'_i,t)
                auto [q_i, v_i] = math_ndsolve :: step(math_ndsolve :: EULER, f_i / m_i, dt, q0_i, v0_i);
                
                q(i) = q_i;
                v(i) = v_i;
            }

            universe.set_q(q);
            universe.set_v(v);
        }

        /**
         * @brief Evolves the universe forward in time
         * 
         * @param universe The universe to step forward in time.
         * @param dt The small timestep to step the universe forward.
         */
        void evolve_step(Universe &universe, long double dt){
            //Calculate a superposition of forces acting on all particles
            auto forces_super = superposition_force(universe);


            //Step particles position asnd velocity forward in in time based on the superposition of forces
            evolve_step_forces(forces_super,universe,dt);

            //Update the universe clock
            universe.clock += dt;
        };
    }
}



#endif
