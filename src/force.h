#ifndef FORCE_H // include guard
#define FORCE_H


#include "series.h"
#include "particle.h"

using namespace physics;

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
             * @param p1 The particle the force is being applied to
             * @return Vector3d The force applied to given particle
             */
            virtual Vector3d get_force(Particle *p1) = 0;
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
            Vector3d acel;

            /**
             * @brief Construct a gravitational force with a aceleration vector
             * 
             * @param dir 
             */
            ForceGravity(Vector3d aceleration){
                acel = aceleration;
            }

            /**
             * @brief Construct a gravitational force in the minus (-) z direction with magnitude of earth's gravity (9.81 m/s^2).
             * 
             */
            ForceGravity(){
                Vector3d a(0,0, - MAG_EARTH);
                acel = a;
            }   


            Vector3d get_force(Particle *p1) override{
                //Fg_i = m * g where g is the magnitude of acelleration
                //And i is the direction of the force
                //Force is 0 in all but this direction
                return acel * p1 -> m;
            }

    };

}

#endif