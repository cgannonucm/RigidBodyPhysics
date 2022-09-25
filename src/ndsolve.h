#ifndef MATH_NDSOLVE_H // include guard
#define MATH_NDSOLVE_H

#include <vector>
#include <functional>
#include "universe.h"
#include "macros.h"

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
        EVector q0, EVector v0, double dt){ 
            
        using namespace Eigen;

        auto f0 = f(q0,v0);

        const double ONE_HALF = 1.0/2.0;
        /*
        From euler method 
        x(t + dt) = x(t) + dt * x'(t) + (1/2)(dt^2) * f(x,xd)
            and
        x'(t+dt) = x'(t) + dt * f(x,x')
        */

        VectorXd q = q0 + dt * v0 + ONE_HALF * pow(dt,2) * f0;
        VectorXd v = v0 + dt * f0; 

        return {q,v};

    }

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
        EVector q0, EVector v0, double dt){ 
            
        using namespace Eigen;

        
        const double ONE2 = 1.0/2.0;
        const double ONE6 = 1.0/6.0;

        double halfdt = dt * ONE2;
        double sixthdt = dt * ONE6;

        //Step 1
        EVector k1 = f(q0,v0);
        EVector v1 = v0 + k1 * halfdt;
        EVector q1 = q0 + ONE6 * halfdt * (halfdt*k1 + 4*v0 + 2*v1);

        //Step 2
        EVector k2 = f(q1,v1);
        EVector v2 = v0 + k2 * halfdt;
        EVector q2 = q0 + ONE6 * halfdt * (halfdt*k1 + 4*v0 + 2*v2);

        //Step 3
        //This part of the formula changes up a bit from first 2 steps
        //This is not a typo
        EVector k3 = f(q2,v2);
        EVector v3 = v0 + k3 * dt;
        EVector q3 = q0 + ONE6 * dt * (dt*k1 + 4*v0 + 2*v3);

        //Step 4
        EVector k4 = f(q3,v3);


        //Put it together
        EVector q = q0 + sixthdt * (v0 + 2 * v1 + 2 * v2 + v3);
        EVector v = v0 + sixthdt * (k1 + 2 * k2 + 2 * k3 + k4);

        return {q,v};

    }


    std :: array<EVector, 2> ndstep(SolveMethod solve_method, 
        funcqv f,
        EVector q0, EVector v0, double dt)
    {
        switch (solve_method)
        {
            case (EULER):
                return ndstep_euler_2d(f,q0,v0,dt);
                break;
            case (RK4):
                return ndstep_rk4_2d(f,q0,v0,dt);
            default:
                break;
        };
        throw std :: invalid_argument("Requested solve method not implemented");
        
    }




            
}

#endif