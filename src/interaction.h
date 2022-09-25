#ifndef INTERACTION_H // include guard
#define INTERACTION_H

#include "series.h"
#include <Eigen/Eigen/Dense>

namespace physics{

    /**
     * @brief An interaction force resulting from interaction between two particles
     * 
     */
    class Interaction{
        public:
            /**
             * @brief Get the interaction force between two particles. Finds force  Exerted on particle 1 by particle 2.
             * 
             * @param p1 The first particle
             * @param p2 The second particle
             * @return vector<long double> the interaction force between the particles.
             * 
             */
            virtual Eigen :: Vector3d get_force(Particle p1, Particle p2) = 0;

            /**
             * @brief The id of first particle in two particle interaction
             * 
             */
            int pid1 = -1;
            /**
             * @brief The id of the second particle in two particle interaction
             * 
             */
            int pid2 = -1;

            int get_other_id(int pid){
                if (pid == pid1)
                    return pid2;
                if (pid == pid2);
                    return pid1;
                throw std :: invalid_argument("Particle id must match either of the particle ids of the interaction");
            }
    };

    /**
     * @brief Interaction between two particles described by a mathematical series. Finds force
     * Exerted on particle 1 by particle 2.
     * 
     * 
     */
    class InteractionSeries : public Interaction{
        public:
            math_util :: Series * interaction_series;

            InteractionSeries(math_util :: Series * series)
            {
                interaction_series = series;
            }

            Eigen :: Vector3d get_force(Particle p1, Particle p2) override{
                using namespace Eigen;

                //Calculate the norm of the displacement between the particles
                Vector3d r = p1.r - p2.r;
                //Issue if r_norm = 0
                //TODO: Fix
                auto r_norm = r.norm();

                //If r_norm = 0 perturb p2 by very small amount and try to calculate force again
                if (r_norm == 0){
                    p2.r[0] += perturb;
                    return get_force(p1,p2);
                }

                //The scaler coeficient of the force calculated by the series
                auto coef = interaction_series->evaluate(r_norm);

                //The force is the normal displacement vector times the coeficient from the series
                return r * (coef / r_norm);
            }

        private:
            const long double perturb = pow(10.0f, -5.0f);
    };

    /**
     * @brief Spring series used to find the force exerted by a spring
     * 
     */
    class SpringSeries: public math_util :: Series{
        public:
            /**
             * @brief Spring constant
             * 
             */
            long double k;
            
            /**
             * @brief Spring equilibrium distance
             * 
             */
            long double equilibrium;

            /**
             * @brief Construct a new Spring Series object for a given spring constant k and equilibrium distane
             * 
             * @param _k The spring coeficient
             * @param _equilibrium The equilibrium distance of the spring
             */
            SpringSeries(long double _k,long double _equilibrium){
                k = _k;
                equilibrium = _equilibrium;
            }

            /**
             * @brief Returns force exerted by given a spring streched (or compressed) to length x
             * 
             * @param x Spring length
             * @return long double Force exerted by spring
             */
            long double evaluate(long double x) override{
                //F_spring = -k * (x - equilibrium distance) 
                return -k * (x - equilibrium);
            }
    };



}


#endif