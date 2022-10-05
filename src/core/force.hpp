#ifndef FORCE_H // include guard
#define FORCE_H

#include <Eigen/Dense>
#include "particle.hpp"

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
    
            ForceGravity(EVectorNd aceleration);

            ForceGravity();

            EVectorNd get_force(Particle p) override;

        


    };

}

#endif
