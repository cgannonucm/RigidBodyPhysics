#ifndef UFORCES_H // include guard
#define UFORCES_H

#include "force.h"
#include "interaction.h"
#include "constraint.h"

namespace physics{

    typedef std :: vector<std :: vector<Interaction *>> v2Interactions;
    typedef std :: vector<Interaction *> vInteractions;
    typedef std :: vector<Force *> vForces;
    typedef std :: vector<Constraint *> vConstraints;

    struct UForces
    {
        
            /**
             * @brief Forces that act on all particles in the universe.
             * Forces are in no particular order.
             * 
             */
            vForces forces;

            /**
             * @brief Interactions that act on two particles in the universe.
             * Interactions are in the same order as particles. 
             * IE ordered like ((I12, I13, I34,...), ...)
             * 
             */
            v2Interactions interactions;

            /**
             * @brief Constraints acting on particles in the universe.
             * Constraints are in no particular order.
             * 
             */
            vConstraints constraints;
    };
    

}

#endif