#ifndef PARTICLE_H // include guard
#define PARTICLE_H

#include <Eigen/Eigen/Dense>

namespace physics
{
    struct Particle
    {
        Eigen :: Vector3d r;
        Eigen :: Vector3d v;
        double m;    
    };
    
}

#endif