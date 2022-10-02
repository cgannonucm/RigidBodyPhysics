#ifndef FORCE_H // include guard
#define FORCE_H

#include <Eigen/Eigen/Dense>
#include "particle.h"

namespace physics{
    
    /**
     * @brief Non interaction force applied to all particles
     * 
     */
    class Force{
        public:
            /**
             * @brief Get force applied on given particle
             * 
             * @param p The particle the force is being applied to
             * @return EVectorNd The force applied to given particle
             */
            virtual EVectorNd get_force(Particle p) = 0;
    };

    /**
     * @brief Gravitational force in the case that the particle acted on is very close to the surface of a planet
     * IE F = m * g
     * 
     */
    class ForceGravity : public Force{
        public:

            /**
             * @brief Magnitude of acceleration near the surface of the earth in units [m/s^2]
             * 
             */
            const long double MAG_EARTH = 9.81;

            /**
             * @brief The direction of the gravitational force
             * 
             */
            EVectorNd acel;

            /**
             * @brief Construct a gravitational force with a aceleration vector
             * 
             * @param dir 
             */
            ForceGravity(EVectorNd aceleration){
                acel = aceleration;
            }

            /**
             * @brief Construct a gravitational force in the minus (-) z direction with magnitude of earth's gravity (9.81 m/s^2).
             * 
             */
            ForceGravity(){
                #if UDim == 3
                EVectorNd a(0,0, - MAG_EARTH);
                #endif
                #if UDim == 2
                EVectorNd a(0,- MAG_EARTH);
                #endif
                acel = a;
            }   


            EVectorNd get_force(Particle p) override{
                //Fg_i = m * g where g is the magnitude of acelleration
                //And i is the direction of the force
                //Force is 0 in all but this direction
                return acel * p.m;
            }

    };

}

#endif
