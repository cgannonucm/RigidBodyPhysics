#include "interaction.hpp"

namespace physics{

    /**
     * Interaction
     * -----------------------------------------
    */

    int Interaction::get_other_id(int pid){
        if (pid == pid1)
            return pid2;
        if (pid == pid2);
            return pid1;
        throw std :: invalid_argument("Particle id must match either of the particle ids of the interaction");
    }

    /**
     * -----------------------------------------
    */

    /**
     * InteractionSeries
     * -----------------------------------------
    */

    InteractionSeries::InteractionSeries(math_util :: Series * series)
    {
        interaction_series = series;
    }

    EVectorNd InteractionSeries::get_force(Particle p1, Particle p2){
        using namespace Eigen;
        //Calculate the norm of the displacement between the particles
        EVectorNd r = p1.r - p2.r;
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

    /**
     * -----------------------------------------
    */

    /**
     * SpringSeries
     * -----------------------------------------
    */

    /**
     * @brief Construct a new Spring Series object for a given spring constant k and equilibrium distane
     * 
     * @param _k The spring coeficient
     * @param _equilibrium The equilibrium distance of the spring
     */
    SpringSeries::SpringSeries(long double _k,long double _equilibrium){
        k = _k;
        equilibrium = _equilibrium;
    }

    /**
     * @brief Returns force exerted by given a spring streched (or compressed) to length x
     * 
     * @param x Spring length
     * @return long double Force exerted by spring
     */
    long double SpringSeries::evaluate(long double x){
        //F_spring = -k * (x - equilibrium distance) 
        return -k * (x - equilibrium);
    }

    /**
     * -----------------------------------------
    */

}