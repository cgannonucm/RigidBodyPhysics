#include "constraint.hpp"

namespace physics{

    /**
     * @brief Construct a new Constraint object constraining particles with the given id
     * 
     * @param _p_ids The particle ids the constraint applies to
     */
    Constraint::Constraint(std :: vector<int> _p_ids)
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
    std :: vector<int> Constraint::get_pids(){
        return p_ids;
    }


    /**
     * @brief Cunstructor for a pivot constraint
     * 
     */
    ConstraintPivot1P::ConstraintPivot1P(Particle &p, double _l, int id) : Constraint({id}){
        l = _l;
    }

    /**
     * Evaluate the constraint
    */
    autodiff :: var ConstraintPivot1P::get(std :: vector<std :: vector<var *>> &q){
        auto r = q[0];
        var constr = pow(*r.at(0),2) + pow(*r.at(1),2) + pow(*r.at(2),2) - pow(l,2);
        //var constr = pow(*r[0],2);
        return constr;
    }


    /**
     * Constructor for a pivot constraint between 2 particles
    */
    ConstraintPivot2P::ConstraintPivot2P(Particle &p1, Particle &p2, int id1, int id2) : Constraint({id1,id2}){
        l = (p1.r - p2.r).norm();
    }

    /**
     * Evaluate the constraint
    */
    autodiff::var ConstraintPivot2P::get(std :: vector<std :: vector<var *>> &q){
        auto r1 = q[0];
        auto r2 = q[1];
        var constr = pow(*r1.at(0) - *r2.at(0),2) +  pow(*r1.at(1) - *r2.at(1),2) +  pow(*r1.at(2) - *r2.at(2),2) -  pow(l,2);
        //var constr = pow(*r[0],2);
        return constr;
    }
}
