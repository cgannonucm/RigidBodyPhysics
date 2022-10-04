#ifndef INTERACTION_H // include guard
#define INTERACTION_H

#include "series.hpp"
#include "macros.hpp"
#include "particle.hpp"
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
            virtual EVectorNd get_force(Particle p1, Particle p2) = 0;

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

            int get_other_id(int pid);
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

            InteractionSeries(math_util :: Series * series);

            EVectorNd get_force(Particle p1, Particle p2) override;

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
            SpringSeries(long double _k,long double _equilibrium);

            /**
             * @brief Returns force exerted by given a spring streched (or compressed) to length x
             * 
             * @param x Spring length
             * @return long double Force exerted by spring
             */
            long double evaluate(long double x) override;
    };



}


#endif