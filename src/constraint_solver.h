#ifndef CONSTRAINT_SOLVER_H // include guard
#define CONSTRAINT_SOLVER_H

#include "universe.h"
#include "constraint.h"
#include <Eigen/Eigen/Dense>

namespace physics{
    using namespace std;
    using namespace Eigen;
    using namespace autodiff;

    class ConstraintSolver{


        public:
            vector<MatrixXd> get_J(vector<var> &c, map<int, vector<var>> &q){
                //Number of constraints and number particles
                int pn = q.size();
                int cn = c.size();

                //J = (grad_j(C_i))_ij
                //Since J is matrix of vectors break it up into 3 matricies
                //Jk = (partial_jk(C_i))_ij
                MatrixXd Jx(pn,cn);
                MatrixXd Jy(pn,cn);
                MatrixXd Jz(pn,cn);
                
                int i = 0;
                int j = 0;

                //Populate matrix with derivatives
                for (auto itr_q = q.begin(); itr_q != q.end(); ++itr_q){
                    auto q_j = &(itr_q -> second);
                    for (var c_i : c){
                        //TODO if partical is not affected by constraint we can evaluate derivatives to 0
                        auto [partial_x,partial_y,partial_z] = derivatives(c_i,wrt(q_j->at(0),q_j->at(1),q_j->at(2)));
                        
                        Jx(i,j) = partial_x;
                        Jy(i,j) = partial_y;
                        Jz(i,j) = partial_z;
                        
                        i++;
                    }
                    j++;
                }

                cout << Jx << endl;
                return {Jx,Jy,Jz};

            }
            
            vector<Vector3d> solve(vector<Vector3d> &F_ex, Universe &u){
                
                auto constraints = u.get_constraints();
                auto particles = u.get_particles();
                
                vector<var> constraint_eval;

                map<int,vector<var>> q;

                //TODO Split up into functions to do each one of the tasks

                //Get q vector
                for (auto itr = particles.begin(); itr != particles.end(); ++itr){
                    auto p = itr->second;
                    
                    auto id =  p->id;
                    auto r = p->r;

                    q[id] = {(var) r[0],(var) r[1], (var) r[2]};
                }
                
                //Evaluate constraints
                for (Constraint *cnstrt : constraints){
                    vector<vector<var *>> q_cnstrt;

                    //Get the q vectors that the particle request
                    for (int p_id : cnstrt->p_ids){

                        //auto q_id = q.find(p_id)->first;
                        vector<var> *q_id = &(q.find(p_id)->second);
                        
                        q_cnstrt.push_back({&q_id->at(0),&q_id->at(1),&q_id->at(2)});
                    }

                    //Evaluate constraint with requested q vector
                    constraint_eval.push_back(cnstrt->get(q_cnstrt));
                }

                //calculate J
                auto J = get_J(constraint_eval,q);

                //We have constraints and q vector, pass off and calculate J
                return {{* new Vector3d(0,0,0)}};

            }

    };

}

#endif