#ifndef VARVECTOR_H // include guard
#define VARVECTOR_H

#include <autodiff/reverse/var.hpp>
#include <array>

namespace physics{

    class VarVector3D{

        //Come back to this
        
        std :: array<autodiff :: var *, 3> r;

        public:
        
            VarVector3D(autodiff :: var *x, autodiff :: var *y, autodiff :: var *z){
                r = {x,y,z};
            }

            
    };

}



#endif