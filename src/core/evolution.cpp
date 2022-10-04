#include "evolution.hpp"

namespace physics{
    
    UEvolver::UEvolver(){
    }

    UEvolver::UEvolver(CSolver _csolver){
        csolver = _csolver;
    }

    /**
     * @brief Evolves the universe forward in time
     * 
     * @param universe The universe to step forward in time.
     * @param dt The small timestep to step the universe forward.
     */
    void UEvolver::evolve_step(Universe &u, double dt){
        using namespace Eigen;
        using namespace math_ndsolve;
        using namespace std;
        using namespace std :: placeholders;
        auto q0 = u.get_q();
        auto v0 = u.get_v();
        auto f = bind(&UEvolver :: get_a,this,_1,_2,ref(u));
        auto [q,v] = ndstep(RK4,f,q0,v0,dt);
        u.set_q(q);
        u.set_v(v);
        //Update the universe clock
        u.clock += dt;
        log(u);
    }


    EVector UEvolver::get_forces(EVector &q, EVector &v,Universe &u){
        using namespace Eigen;
        auto mv = u.get_m();
        auto forces = u.get_forces();
        auto interactions = u.get_interactions();
        auto constraints = u.get_constraints();
        //Calculate a superposition of external forces acting on all particles
        auto fex = superposition_force(q,v,mv,forces,interactions,u.get_p_count());
        //Solve for constraint forces
        auto fc = csolver.solve(q,v,mv,constraints,fex);
        //Sum of force is fex + fc
        VectorXd f = fex + fc;
        return f;
    }

    EVector UEvolver::get_a(EVector &q, EVector &v,Universe &u){
        return CSolver :: get_M_inv(u.get_m()) * get_forces(q,v,u);
    }

    //TODO : Put into class
    /**
     * @brief Adds to the 3d force vector to the 3 * particle count force vector
     * 
     * @param forces_super The forces to add to
     * @param force_to_add The forces to add
     * @param n The id of the particle
     */
    void UEvolver::add_to_superpos(EVector &forces_super, EVectorNd &force_to_add, int n){
        for (int j = 0; j < UDim; j++){
            auto pos = Universe :: get_pos(n,j);
            forces_super(Universe :: get_pos(n,j)) += force_to_add[j];
        }
    }

    /**
     * @brief Sums up interactions for a particle adds to superposition of forces.
     * 
     * @param forces The forces acting on the particle
     * @param p1 The particle to get forces on
     * @param n The particle id
     * @param forces_super The vector to add the forces to
     */
    void UEvolver::superposition_force_p(vForces &forces, Particle &p1,
        int n, EVector &forces_super){
        //Loop through forces and each interaction force
        for (Force *f : forces){
            //Evaluate force and add it's components to our force vector
            auto f_val = f->get_force(p1);
            add_to_superpos(forces_super, f_val, n);
        }
    }

    /**
     * @brief Sums up interactions for a particle adds to superposition of forces.
     * 
     * @param u The universe 
     * @param p The particle to sum forces for
     * @param n The id of the particle
     * @param forces_super The vector to add the interaction forces to
     *
     */
    void UEvolver::superposition_interaction_p(EVector &q, EVector &v, EVector &m,
        Particle &p, v2Interactions const &interactions_all, EVector &forces_super)
    {
        int pid = p.id;
        std :: vector<Interaction *> interactions = Universe :: get_p_interactions(pid,interactions_all);
        for (Interaction * interaction : interactions){
            
            int pid2 = interaction->get_other_id(p.id);
            auto p2 = Universe :: get_p(pid,q,v,m);
            
            //Get force calculates force on first particle
            //p is always first particle since we are calculating the force on it
            auto f_val = interaction->get_force(p,p2);
            add_to_superpos(forces_super,f_val,pid);
        }
    }


    /**
     * @brief Gets a superposition of all non constraint forces and interactions
     * 
     * @param universe 
     * @return  
     */
    EVector UEvolver::superposition_force(EVector &q, EVector &v, EVector &m, vForces &forces,
       v2Interactions &interactions, int p_count){
            
        using namespace Eigen;
        VectorXd forces_super = VectorXd :: Zero(UDim * p_count);
        //Loop through all particles in universe evaluate interactions and forces for all particles
        //And add to superposition-
        for (int n = 0; n < p_count; n++){
            Particle p1 = Universe :: get_p(n,q,v,m);
            //Add non interaction forces to the superposition of forces
            superposition_force_p(forces,p1,n,forces_super);
            //TODO future optimization by newtons law we only need to calculate interaction
            //force for one particle in the interaction and we know the force on the other
            //Add interaction forces to the superposition of forces
            superposition_interaction_p(q,v,m,p1,interactions,forces_super);
        }
        return forces_super;
    };
    /**
     * @brief Expand m vector to DIMENSION * particle count. Returned m vector will be (m1, m1, m1, ... , mn, mn, mn)
     * 
     * @param universe
     * @return Eigen 
     */
    EVector UEvolver::expand_m_vector(Universe &u){
        using namespace Eigen;
        auto p_count = u.get_p_count();
        auto m_small = u.get_m();
        VectorXd m_large = VectorXd :: Zero(UDim);
        for (int n = 0; n < p_count; n++){
            for (int j = 0; j < UDim;  j++){
                m_large(Universe :: get_pos(n,j)) = m_small(n);
            }
        }
        return m_large;
    }

    void UEvolver::log(Universe &u){
        for (auto logptr : loggers){
            logptr->log(u);
        }
    }
}



