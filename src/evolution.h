#ifndef EVOLUTION_H // include guard
#define EVOLUTION_H

#include "ndsolve.h"

using namespace math_ndsolve;

namespace physics_evolution{
    
    /**
     * @brief Gets a superposition of all forces and interactions
     * 
     * @param universe 
     * @return Vector3d<long double> 
     */
    vector<Vector3d> superposition_force(Universe &universe){
        auto particles = universe.get_particles();
        auto forces = universe.get_forces();
        auto interactions = universe.get_interactions();

        vector<Vector3d> superpos_all;

        //loop through all particles
        for (auto itr = particles.begin(); itr != particles.end(); ++itr){
            Vector3d superpos(0,0,0);

            auto id = itr->first;
            auto p = itr->second;
            

            //Sum all forces and interactions

            //Sum non interaction forces
            for (auto & force : forces){
                superpos += force->get_force(p);
            }

            //Select interactions effecting our particle
            auto select = interactions.equal_range(id);

            //Sum interaction force
            //This is slower than needs to be. By newton's second law if we have force on 1 particle
            //In 2 particle interaction we have force on other.
            for (auto _itr = select.first; _itr != select.second; ++_itr){
                auto interaction = _itr->second;


                Particle *p2;

                //Find second particle in interaction
                if (id == interaction->pid1){
                   p2 = particles.find(interaction->pid2)->second;
                }
                else{
                   p2 = particles.find(interaction->pid1)->second;
                }

                auto add_force = interaction -> get_force(p,p2);
                superpos += add_force;
            }

            superpos_all.push_back(superpos);

        }

        return superpos_all;

    };

    void evolve_step(Universe &universe, long double dt){
        //Calculate a superposition of forces acting on all particles
        auto forces_super = superposition_force(universe);

        auto particles = universe.get_particles();

        //Step differential equation forward in time for each particle
        int n = 0;
        for (auto itr = particles.begin(); itr != particles.end(); ++itr){
            auto particle = itr->second;
            auto force = forces_super[n];
            
            step_force(EULER,force,dt,particle);
            n++;
        }

        //Update the universe clock
        universe.clock += dt;
    };



}



#endif
