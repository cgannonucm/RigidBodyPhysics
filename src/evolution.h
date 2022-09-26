#ifndef EVOLUTION_H // include guard
#define EVOLUTION_H

#include "ndsolve.h"
#include "constraint_solver.h"
#include <chrono>
#include "macros.h"
#include "logger.h"

namespace physics{
    class UniverseSolver{
        public:

            long perft_JJd = 0;

            std :: vector<Logger *> loggers; 

            ConstraintSolver csolver;

            UniverseSolver(){

            }

            UniverseSolver(ConstraintSolver _csolver){
                csolver = _csolver;
            }

            /**
             * @brief Evolves the universe forward in time
             * 
             * @param universe The universe to step forward in time.
             * @param dt The small timestep to step the universe forward.
             */
            void evolve_step(Universe &u, double dt){
                using namespace Eigen;
                using namespace math_ndsolve;
                using namespace std;
                using namespace std :: placeholders;

                auto q0 = u.get_q();
                auto v0 = u.get_v();

                auto f = bind(&UniverseSolver :: get_a,this,_1,_2,ref(u));
                auto [q,v] = ndstep(RK4,f,q0,v0,dt);

                u.set_q(q);
                u.set_v(v);

                //Update the universe clock
                u.clock += dt;

                log(u);

            }

            EVector get_forces(EVector &q, EVector &v,Universe &u){
                using namespace Eigen;
                typedef std::chrono::high_resolution_clock Clock;

                auto mv = u.get_m();
                auto forces = u.get_forces();
                auto interactions = u.get_interactions();
                auto constraints = u.get_constraints();


                //Calculate a superposition of external forces acting on all particles
                auto fex = superposition_force(q,v,mv,forces,interactions,u.get_p_count());

                //Solve for constraint forces
                auto fc = csolver.solve(q,v,mv,constraints,fex);

                //Sum of force is fex + fc
                VectorXd f = fex + fc;

                return f;
            }

            EVector get_a(EVector &q, EVector &v,Universe &u){
                return ConstraintSolver :: get_M_inv(u.get_m()) * get_forces(q,v,u);
            }

        private:
            //TODO : Put into class
            /**
             * @brief Adds to the 3d force vector to the 3 * particle count force vector
             * 
             * @param forces_super The forces to add to
             * @param force_to_add The forces to add
             * @param n The id of the particle
             */
            void add_to_superpos(Eigen :: VectorXd &forces_super, Eigen :: Vector3d &force_to_add, int n){
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
            void superposition_force_p(std :: vector<Force *> &forces, Particle &p1,
                int n, Eigen :: VectorXd &forces_super){
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
            void superposition_interaction_p(EVector &q, EVector &v, EVector &m,
                Particle &p, std :: vector<std :: vector<Interaction *>> const &interactions_all, EVector &forces_super)
            {
                int pid = p.id;

                std :: vector<Interaction *> interactions = Universe :: get_p_interactions(pid,interactions_all);

                for (Interaction * interaction : interactions){
                    
                    int pid2 = interaction->get_other_id(p.id);
                    auto p2 = Universe :: get_p(pid,q,v,m);
                    
                    //Get force calculates force on first particle
                    //p is always first particle since we are calculating the force on it
                    auto f_val = interaction->get_force(p,p2);
                    add_to_superpos(forces_super,f_val,pid);
                }
            }

            /**
             * @brief Gets a superposition of all non constraint forces and interactions
             * 
             * @param universe 
             * @return Vector3d<long double> 
             */
            EVector superposition_force(EVector &q, EVector &v, EVector &m, std :: vector<Force *> &forces,
                std :: vector<std :: vector<Interaction *>> &interactions, int p_count){
                    
                using namespace Eigen;

                VectorXd forces_super = VectorXd :: Zero(Universe :: DIMENSION * p_count);

                //Loop through all particles in universe evaluate interactions and forces for all particles
                //And add to superposition-
                for (int n = 0; n < p_count; n++){
                    Particle p1 = Universe :: get_p(n,q,v,m);

                    //Add non interaction forces to the superposition of forces
                    superposition_force_p(forces,p1,n,forces_super);


                    //TODO future optimization by newtons law we only need to calculate interaction
                    //force for one particle in the interaction and we know the force on the other
                    //Add interaction forces to the superposition of forces
                    superposition_interaction_p(q,v,m,p1,interactions,forces_super);
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

            void log(Universe &u){
                for (auto logptr : loggers){
                    logptr->log(u);
                }
            }
    };
}



#endif
