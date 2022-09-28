//Physics universe
#ifndef PHYSICS_UNIVERSE_H // include guard
#define PHYSICS_UNIVERSE_H

#include "force.h"
#include "interaction.h"
#include "constraint.h"
#include "macros.h"


namespace physics
{

    enum Coordinate{
        X = 0,
        Y = 1,
        Z = 2
    };

    /**
     * @brief A universe in which particles forces and interactions exist.
     * 
     */
    class Universe{
        public:
            
            double clock = 0;

            //static const int DIMENSION = 3;

            Universe(int n){
                init(n);
            }

            /**
             * @brief Sets a particle's mass velocity and r vector
             * 
             * @param n The particle to set
             * @param _r The displament vector
             * @param _v The velocity vector
             * @param _m The mass of the particle
             */
            void set_particle(EVectorNd _r, EVectorNd _v, double _m, int n){
                assert(n < p_count && n >= 0);

                //Store each position / velocity as an entry in the q and v vector
                for (int i = 0; i < UDim; i++){
                    q(UDim * n + i) = _r[i];
                    v(UDim * n + i) = _v[i];
                }

                m(n) = _m;
            }

            /**
             * @brief Sets the particle at n to have specified properties
             * 
             * @param p Tbe properties of the particle
             * @param n The id of the particle
             */
            void set_particle(Particle p, int n){
                set_particle(p.r,p.v,p.m,n);
            }

            /**
             * @brief Adds a particle with given properties.
             * Starts out at particle 0. Every particle add increments the particle position.
             * (IE) next particle added will be particle 1.
             * 
             * @param p The particle properties to add
             * @return int The particle id of the added particle
             */
            int add_particle(Particle p){
                set_particle(p,c_id);
                c_id++;
                return c_id - 1;
            }

            /**
             * @brief Gets the r coordinates of specified particle
             * 
             * @param pid The particle id to get coordinates of
             * 
             * @return EVectorNd The r coordinates
             */
            EVectorNd get_p_r(int pid){
                return Universe :: get_p_r(pid,q);
            }
            
            static EVectorNd get_p_r(int pid, EVector &_q){
                assert(pid * (UDim + 1) <= _q.size() && pid >= 0);
                EVectorNd r = EVectorNd :: Zero();
                for (int j = 0; j < UDim; j++)
                    r(j) = _q(Universe :: get_pos(pid,j));

                return r;
            }

            /**
             * @brief Gets the v coordinates of specified particle
             * 
             * @param pid The particle id to get the v coordinates of
             * 
             * @return EVectorNd The v coordinates
             */
            EVectorNd get_p_v(int pid){
                return Universe :: get_p_v(pid,v);
            }

            static EVectorNd get_p_v(int pid, EVector &_v){
                return Universe :: get_p_r(pid,_v);
            }


            /**
             * @brief Get the ptr to m coordinates of specified particle
             * 
             * @param pid The id of the particle to get the coordinates of.
             * @return double* The pointer to the mass
             */
            double get_p_m(int pid){
                return Universe :: get_p_m(pid,m); 
            }

            static double get_p_m(int pid, EVector &_m){
                assert(pid < _m.size());
                return _m(pid);
                
            }

            /**
             * @brief Add a force to the universe 
             * 
             * @param force The force to add
             */
            void add_force(Force *force){
                forces.push_back(force); 
            }


            /**
             * @brief Get all forces in the universe
             * 
             * @return vector<Force *> Forces in the universe
             */
            std :: vector<Force *> get_forces(){
                return forces;
            }

            /**
             * @brief Add an interaction to the universe acting on 2 particles
             * 
             * @param i The interaction to add
             * @param pid1 The first particle in the interaction
             * @param pid2 The second particle in the interaction
             */
            void add_interaction(Interaction *i, int pid1, int pid2){
                //Set interaction particle ids
                i->pid1 = pid1;
                i->pid2 = pid2;

                //Assign this interaction to both particle
                interactions.at(pid1).push_back(i);
                interactions.at(pid2).push_back(i);
            }

            /**
             * @brief Gets the interactions in the universe
             * 
             * @return Interactions in the universe
             * 
             */
            std :: vector<std :: vector<Interaction *>> get_interactions(){
                return interactions;
            }

            std :: vector<Interaction *> get_p_interactions(int pid){
                return get_p_interactions(pid,interactions);
            }

            static std :: vector<Interaction *> get_p_interactions(int pid, std :: vector<std :: vector<Interaction *>> _interactions){
                return _interactions.at(pid);
            }

            /**
             * @brief Adds a constraint to the universe
             * 
             * @param constraint The constraint to add
             */
            void add_constraint(Constraint * constraint)
            {
                //TODO Check to make sure universe is not overconstrained
                constraints.push_back(constraint); 
            }

            /**
             * @brief returns all constraints in the universe
             * 
             * @return vector<Constraint *> The constraints in the universe.
             */
            std :: vector<Constraint *> get_constraints()
            {
                return constraints;
            }


            static Particle get_p(int pid, EVector &_q, EVector &_v, EVector &_m){
                Particle p;
                p.r = get_p_r(pid,_q);
                p.v = get_p_v(pid,_v);
                p.m = get_p_m(pid,_m);
                p.id = pid;
                return p;
            }

            /**
             * @brief Gets all the properties of the particle.
             * 
             * @param pid The particle to get properties of.
             * @return Particle the properties of the particle.
             */
            Particle get_p(int pid){
                return get_p(pid,q,v,m);
            }

            /**
             * @brief Returns q vector containing the positions of all particles in the universe.
             * Is 3n dimension vector (n is number of particles) with the particle
             * positions stored as (x1,y1,z1,...,xn,yn,zn)
             * 
             * @return q vector
             */
            EVector get_q(){
                return q;
            }

            /**
             * @brief Gets the v vector containing the velocities of all particles in the universe.
             * Is 3n dimension vector (n is number of particles) with the particle
             * velocities stored as (v_x1,v_y1,v_z1,...,v_xn,v_yn,v_zn)
             * 
             * @return v vector
             */
            EVector get_v(){
                return v;
            }

            /**
             * @brief Gets the m vector containing the masses of all particles in the universe.
             * Is n dimension vector (n is the number of particles) with the particles
             * masses stored as (m1, m2, ... mn)
             * 
             * @return Eigen 
             */
            EVector get_m(){
                return m;
            }


            /**
             * @brief Set the q vector
             * 
             * @param _q The vector to set q to
             */
            void set_q(EVector _q){
                assert(_q.size() == UDim * p_count);
                q = _q;
            }


            /*
             * @brief Set the v vector
             * 
             * @param _v The vector to set v to
             */
            void set_v(EVector _v){
                assert(_v.size() == UDim * p_count);
                v = _v;
            }

            /**
             * @brief Gets number of particles in the universe
             * 
             * @return int Number of particles in the universe 
             */
            int get_p_count(){
                return p_count;
            }


            /**
             * @brief Returns the position in the q / v vector of the specified particle / coordinate
             * 
             * @param n The particle id
             * @param j The coordinate : 0 for x, 1 for y, 2 for z
             * @return int The position in the vector of the particle / coordinate
             */
            static int get_pos(int n, int j){
                return UDim * n + j;
            }
            
        
        private:
            int c_id = 0;

            int p_count;

            /**
             * @brief Called when creating the universe
             * 
             * @param n Number of particles in the universe
             */
            void init(int n){
                p_count = n;
                using namespace Eigen;

                q = VectorXd :: Zero(UDim * n);
                v = VectorXd :: Zero(UDim * n);
                m = VectorXd :: Ones(n);

                //The size of the interactions vector is exactly one
                interactions.resize(n);
            }


            /**
             * @brief Forces that act on all particles in the universe.
             * Forces are in no particular order.
             * 
             */
            std :: vector<Force *> forces;

            /**
             * @brief Interactions that act on two particles in the universe.
             * Interactions are in the same order as particles. 
             * IE ordered like ((I12, I13, I34,...), ...)
             * 
             */
            std :: vector<std :: vector<Interaction *>> interactions;

            /**
             * @brief Constraints acting on particles in the universe.
             * Constraints are in no particular order.
             * 
             */
            std :: vector<Constraint *> constraints;

            /**
             * @brief Vector containing the positions of all particles in the universe.
             * Is 3n dimension vector (n is number of particles) with the particle
             * positions stored as (x1,y1,z1,...,xn,yn,zn)
             * 
             */
            EVector q;

            /**
             * @brief Vector containing the velocities of all particles in the universe.
             * Is 3n dimension vector (n is number of particles) with the particle
             * velocities stored as (v_x1,v_y1,v_z1,...,v_xn,v_yn,v_zn)
             * 
             */
            EVector v;

            /**
             * @brief Vector containing the masses of all particles in the universe.
             * Is n dimension vector (n is the number of particles) with the particles
             * masses stored as (m1, m2, ... mn)
             * 
             */
            EVector m;

            /**
             * @brief Construct a 3D universe containing n particles.
             * All particles default to 0 vectors for displacement and velocity and a mass of one.
             * 
             * @param n The number of particles in the universe  
             */
            

    };


}

#endif