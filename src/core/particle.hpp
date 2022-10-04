#ifndef PARTICLE_H // include guard
#define PARTICLE_H

#include <Eigen/Eigen/Dense>
#include "macros.hpp"

namespace physics
{
    struct Particle
    {
        EVectorNd r;
        EVectorNd v;
        double m;
        int id = -1;
    };
    
}

#endif