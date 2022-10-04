#include <Eigen/Eigen/Dense>
#include "force.hpp"


namespace physics{


    EVectorNd ForceGravity::get_force(Particle p){
        //Fg_i = m * g where g is the magnitude of acelleration
        //And i is the direction of the force
        //Force is 0 in all but this direction
        return acel * p.m;
    }

    /**
    * @brief Construct a gravitational force with a aceleration vector
    * 
    * @param dir 
    */
    ForceGravity::ForceGravity(EVectorNd aceleration){
        acel = aceleration;
    }
    
    /**
     * @brief Construct a gravitational force in the minus (-) z direction with magnitude of earth's gravity (9.81 m/s^2).
     * 
     */
    ForceGravity::ForceGravity(){
        #if UDim == 3
        EVectorNd a(0,0, - MAG_EARTH);
        #endif
        #if UDim == 2
        EVectorNd a(0,- MAG_EARTH);
        #endif
        acel = a;
    }   

}