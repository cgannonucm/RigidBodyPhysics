#ifndef VARVECTOR_H // include guard
#define VARVECTOR_H

#include <autodiff/reverse/var.hpp>
#include <array>

namespace physics{

    class VarVectorND{

        //Come back to this
        
        std :: array<autodiff :: var *, 3> r;

        public:
        
            VarVectorND(autodiff :: var *x, autodiff :: var *y, autodiff :: var *z){
                r = {x,y,z};
            }

            
    };

}



#endif