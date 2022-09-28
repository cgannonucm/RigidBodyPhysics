#ifndef CONSTRAINT_H // include guard
#define CONSTRAINT_H

#include <iostream>
#include <autodiff/reverse/var.hpp>
#include "particle.h"
#include <algorithm>

using namespace autodiff;

namespace physics{
    

    /**
     * @brief Represents a holominic constraint placed on particles in the universe.
     * A holonomic constraint C(q) evaluates to 0 when satisfied.
     * Where q is the generalized coordinate vector of
     * 
     * https://en.wikipedia.org/wiki/Holonomic_constraint
     * 
     */
    class Constraint{

        public: 

            /**
             * @brief Returns the holonomic constraint.  
             * 
             * @param q Vector of all coordinates vectors of particles in the universe that have an
             * id contained in _p_ids. This vector is in the same order as _p_id.
             * Constraints must be soley expressed in terms of var and the mathematical operations of autodiff 
             * also constraints must have smooth partials up to 2nd partial derivative for all variable. 
             * 
             * @return var The contstraint function
             */
            virtual var get(std :: vector<std :: vector<var *>> &q) = 0;

            /**
             * @brief Construct a new Constraint object constraining particles with the given id
             * 
             * @param _p_ids The particle ids the constraint applies to
             */
            Constraint(std :: vector<int> _p_ids)
            {
                USING_STANDARD_NAMESPACES;
                
                p_ids = _p_ids;
                sort(p_ids.begin(),p_ids.end());
            }

            /**
             * @brief Returns a vector of the particles in the interaction in assending order
             * 
             * @return (see description)
             */
            std :: vector<int> get_pids(){
                return p_ids;
            }

            
        
        private: 

            std :: vector<int> p_ids;

    };


    /**
     * @brief A pivot constraint
     * 
     */
    class ConstraintPivot1P : public Constraint{
        public:

            ConstraintPivot1P(Particle &p, double _l, int id) : Constraint({id}){
                l = _l;
            }

            var get(std :: vector<std :: vector<var *>> &q) override{
                auto r = q[0];
                var constr = pow(*r.at(0),2) + pow(*r.at(1),2) + pow(*r.at(2),2) - pow(l,2);
                //var constr = pow(*r[0],2);
                return constr;
            };

        private:
            var l;
            
    };

    class ConstraintPivot2p : public Constraint{
        public:


            ConstraintPivot2p(Particle &p1, Particle &p2, int id1, int id2) : Constraint({id1,id2}){
                l = (p1.r - p2.r).norm();
            }

            var get(std :: vector<std :: vector<var *>> &q) override{
                auto r1 = q[0];
                auto r2 = q[1];

                var constr = pow(*r1.at(0) - *r2.at(0),2) +  pow(*r1.at(1) - *r2.at(1),2) +  pow(*r1.at(2) - *r2.at(2),2) -  pow(l,2);
                //var constr = pow(*r[0],2);
                return constr;
            };

            private: 
                var l;

    };
}


#endif