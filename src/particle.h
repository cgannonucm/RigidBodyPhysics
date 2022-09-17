#ifndef PARTICLE_H // include guard
#define PARTICLE_H

#include <Eigen/Eigen/Dense>

using Eigen::Vector3d;
using Eigen::VectorXd;

namespace physics
{
    /**
     * @brief A physics particle. Described by 3 properties: displacement, velocity and mass.
     * 
     */
    class Particle{

        public:
            /**
             * @brief Particles's displacement from the origin
             * 
             */
            Vector3d r;

            /**
             * @brief Particle's velocity vector
             * 
             */
            Vector3d v;

            /**
             * @brief Particle's mass
             * 
             */
            long double m;

            /**
             * @brief Particle ID asigned when added to the universe
             * 
             */
            int id = -2; 
    };
}

#endif