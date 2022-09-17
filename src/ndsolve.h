#ifndef MATH_NDSOLVE_H // include guard
#define MATH_NDSOLVE_H

#include <vector>
#include "universe.h"

using namespace std;
using namespace physics;

/**
 * @brief This namespace includes tools for numerical solution of ordinary differential equations
 * 
 */
namespace math_ndsolve{
    /**
     * @brief Enum with the different numerical solver methods available
     * 
     */
    enum SolveMethod{
        EULER
    };

    /**
     * @brief Steps a differential equation, x'' = f(x,x')
     * forward in time using numerical euler method
     * https://en.wikipedia.org/wiki/Euler_method
     * 
     * @param f The right side of the equation x(t)'' = f(x(t),x(t)')
     * @param dt The time to step forward. Smaller -> better
     * @param xxd A two dimensional vecotr of form (x(t),x(t)')
     * @return long double* A vector of form (x(t + dt),x(t + dt)')
     */
    vector<long double> step_euler_2d(long double f, long double dt, vector<long double> xxd){
            long double x0 = xxd[0];
            long double xd0 = xxd[1];
            const long double ONE_HALF = 1.0/2.0;
            /*
            From euler method 
            x(t + dt) = x(t) + dt * x'(t) + (1/2)(dt^2) * f(x,xd)
                and
            x'(t+dt) = x'(t) + dt * f(x,x')
            */
            long double x = x0 + dt * xd0 + ONE_HALF * pow(dt,2) * f;
            long double xd = xd0 + dt * f;

            return {x,xd};
    };   

    /**
     * @brief Steps the equation m*r'' = f(r,r')  forward in time for a particle using euler method
     * 
     * @param force Force acting on particle
     * @param dt Timestep - smaller = better
     * @param particle The particle to step forward in time
     */
    void step_force_euler(Vector3d force,long double dt,Particle *particle){
        //Evolve each of the compontents of the particle position and velocity forward in time using the euler method
        for(int i = 0; i < 3; i++){
            //mx''_i = F_i -> x''_i = 1/m * F_i . 
            //We need to solve the equation x''_i = 1/m * F_i
            auto xxd = step_euler_2d( (1 / particle->m) * force[i],dt,{particle->r[i],particle->v[i]});

            particle->r[i] = xxd[0];
            particle->v[i] = xxd[1];
        };
    };

    /**
     * @brief Steps the vector equation m*r'' = f(r,r')  forward in time for a particle 
     * 
     * @param force The force on the object evaluated a current time
     * @param particle The particle to step forward 
     */
    void step_force(SolveMethod solve_method,Vector3d force,long double dt,Particle *particle){
        //Choose solution method
        switch (solve_method)
        {
            case (EULER):
                step_force_euler(force,dt,particle);
                break;
            default:
                break;
        };
    };




            
}

#endif