#ifndef MATH_NDSOLVE_H // include guard
#define MATH_NDSOLVE_H

#include <vector>
#include <functional>
#include "macros.hpp"

namespace physics{

    /**
     * @brief This namespace includes tools for numerical solution of ordinary differential equations
     * 
     */
    namespace math_ndsolve{
        //TODO put into class
        /**
         * @brief Enum with the different numerical solver methods available
         * 
         */
        enum SolveMethod{
            EULER,
            RK4
        };

        std :: array<EVector, 2> ndstep_euler_2d(funcqv f, 
            EVector q0, EVector v0, double dt);

        /**
         * @brief Steps a differential equation, x'' = f(x,x')
         * forward in time using numerical rk4 method
         * https://willbeason.com/2021/06/25/improve-your-runge-kutta-nystrom-algorithms-with-this-one-weird-trick/
         * Note there is a typo in article in many of the formulas y_0 should be y'_0
         * 
         * @param f The right side of the equation x(t)'' = f(x(t),x(t)')
         * @param dt The time to step forward. Smaller -> better
         * @param x0 The x variable
         * @param xd0 the x' variable
         * @return long double* A vector of form (x(t + dt),x(t + dt)')
         */
        std :: array<EVector, 2> ndstep_rk4_2d(funcqv f, 
            EVector q0, EVector v0, double dt);

        std :: array<EVector, 2> ndstep(SolveMethod solve_method, 
            funcqv f,
            EVector q0, EVector v0, double dt);




    }
}

#endif