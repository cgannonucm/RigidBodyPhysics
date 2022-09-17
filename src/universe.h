//Physics universe
#ifndef PHYSICS_UNIVERSE_H // include guard
#define PHYSICS_UNIVERSED_H

#include "force.h"
#include "interaction.h"

using namespace std;
using namespace physics; 

namespace physics
{

    /**
     * @brief A universe in which particles exists. 
     * A universe also includes a collection of forces 
     * 
     */
    class Universe{
        public:

            /**
             * @brief Tracks how long the universe has been evolving for
             * 
             */
            long double clock = 0;

            /**
             * @brief Add a particle to Universe
             * 
             * @param p The particle to add to the universe
             * @return int The id of the particles
             */
            int add_particle(Particle *p){
                particles.insert(pair<int,Particle *>(c_id,p));
                p->id = c_id;
                c_id++;
                return p->id;
            }

            /**
             * @brief Adds an interaction between two particles with given ids
             * 
             * @param interaction The interaction to add
             */
            void add_interaction(int pid1, int pid2, Interaction *interaction){
                interaction->pid1 = pid1;
                interaction->pid2 = pid2;
                 
                interactions.insert(pair<int,Interaction *>(pid1,interaction));
                interactions.insert(pair<int,Interaction *>(pid2,interaction));
            }

            /**
             * @brief Adds a non-interaction force to the universe
             * 
             * @param force Force to add
             */
            void add_force(Force *force){
                forces.push_back(force);
            }

            /**
             * @brief Gets the particles in the universe
             * 
             * @return map<int,Particle> The particles in the universe
             */
            map<int,Particle *> get_particles(){
                return particles;
            }

            /**
             * @brief Get the interactions forces in the universe
             * 
             * @return map<int,Particle> Interaction forces in the univerce
             */
            multimap<int,Interaction *> get_interactions(){
                return interactions;
            }

            /**
             * @brief Get the forces in the univerces
             * 
             * @return vector<Force> forces in the universe
             */
            vector<Force *> get_forces(){
                return forces;
            }

            /**
             * @brief Get the particles that are participating in the interaction
             * 
             * @param interaction The interaction between the particles
             * @return vector<Particle> The particles in the interaction
             */
            vector<Particle *> get_particles_interaction(Interaction *interaction){
                return  {particles[interaction->pid1], particles[interaction->pid2]};
            }

            /**
             * @brief Construct an empty universe
             * 
             */
            Universe(){

            }
            
        private:
            /**
             * @brief The id number of the last particle added
             * 
             */
            int c_id = 0;

            /**
             * @brief The particles in the universe, labeled by an id
             * 
            */
            map<int,Particle *> particles;


            /**
             * @brief The forces in the universe that are caused by an interaction between two particles.
             * 
             */
            multimap<int,Interaction *> interactions; 


            /**
             * @brief The forces in the universe that act on all particles uniformly.
             * 
             */
            vector<Force *> forces;

    };


}

#endif