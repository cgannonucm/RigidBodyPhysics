#ifndef CONSTRAINT_SOLVER_H // include guard
#define CONSTRAINT_SOLVER_H

#include "universe.h"
#include "constraint.h"
#include <Eigen/Eigen/Dense>

namespace physics{

    class ConstraintSolver{ 
        public:
            /**
             * @brief Gets J and J dot matrix for constraint equation (Jacobian and time derivative of jacobian).
             * Both matricies are size Dimension * n by m where n is number of particles m is number of constraints.
             * 
             * @param u The universe
             * @param J Set to the Jacobian matrix of constraints
             * @param Jd Set to the time derivative of the Jacobian matrix
             * 
             * @returns {J,Jd}
             */
            std :: array<Eigen :: MatrixXd,2> get_JJd(Universe &u){
                using namespace std;
                using namespace Eigen;
                using namespace autodiff;
                
                //Get variables everything ready
                auto constraints = u.get_constraints();

                auto q = u.get_q();
                auto v = u.get_v();
                
                int cnst_count = constraints.size();
                int rows = cnst_count;
                int p_count = u.get_p_count();
                int cols = Universe :: DIMENSION * p_count;
                    

                MatrixXd J = MatrixXd :: Zero(rows,cols);
                
                //Auxilary matricies K_j, used to calculate Jd
                //K is nonstandard naming convention because I don't know what
                //it is actually called
                vector<MatrixXd> K = init_Kmats(rows,cols);

                vector<var> qvar = get_qvar(q);

                vector<var> c_eval = eval_constraints(constraints,qvar);

                //Calculate J matrix and K matricies
                //J = (partial_j(C_i))_ij
                //K_j = (partial_i(partial_k(C_j)))_ik
                for (int i = 0; i < cnst_count; i++){
                    var *c = &c_eval.at(i);

                    auto p_ids = constraints.at(i)->p_ids;

                    for (int n : p_ids){
                        
                        auto grad = grad_j(n,c,qvar);
                        add_grad_to_mat(i,n,grad,&J);

                        //Calculate K matrix
                        //We need to calculate second partial derivatives
                        //Todo since partial_i(partial_j(C)) =  partial_j(partial_i(C))
                        //IE K is symetric
                        //we can optimize by elimating these redundant calculations

                        for (int m = 0; m < Universe :: DIMENSION; m++){
                            for (int n2 : p_ids){
                                int j = Universe :: DIMENSION * n + m;
                                auto grad2 = grad_j(n2,&grad.at(m),qvar);
                                add_grad_to_mat(j,n2,grad2,&K.at(i));
                            }                 
                        }
                    }
                }

                MatrixXd Jd = get_Jd(K,v);

                cout << J << endl;
                cout << Jd << endl; 

                return {J,Jd};
            }


            /**
             * @brief Initializes mcount number of dim * dim matricies filled with 0s
             * 
             * @param mcount (see description)
             * @param dim (see description)
             * @return std (see description)
             */
            std :: vector<Eigen :: MatrixXd> init_Kmats(int mcount, int dim){
                using namespace std;
                using namespace Eigen;

                vector<MatrixXd> K;
                K.resize(mcount);
                 
                for (int i = 0; i < mcount; i++){
                    K[i] = MatrixXd :: Zero(dim, dim);
                }

                return K;
            }


            /**
             * @brief Convert from VectorXd q to vector<var> qvar to prepare for diferentiation
             * 
             * @param q Generalized position variable
             * @return std (see description)
             */
            std :: vector<autodiff :: var> get_qvar(Eigen :: VectorXd q){
                using namespace :: std;

                auto cols = q.size();

                vector<var> qvar;
                qvar.resize(cols);

                for (int i = 0; i < cols; i++){
                    qvar[i] = (var) q(i);
                }

                return qvar;
            }


            /**
             * @brief Evaluates the constraints in preperation for diferentiation
             * 
             * @param constraints The constraints to evaluate
             * @param qvar The q vector to evaluate constraints at
             * @return std The evaluated constraints
             */
            std :: vector<autodiff :: var> eval_constraints(std :: vector<Constraint *> &constraints, std :: vector<var> &qvar){
                using namespace std;
                using namespace Eigen;
                using namespace autodiff;
                
                auto cnst_count = constraints.size();

                vector<var> c_eval;
                c_eval.resize(cnst_count);

                for (int i = 0; i < cnst_count; i++){
                    Constraint *cnst = constraints.at(i);

                    //Only pass requested coordintates 
                    vector<vector<var *>> qcnst;
                    for (int pid : cnst->p_ids){
                        
                        vector<var *> q_pid;
                        q_pid.resize(Universe :: DIMENSION);

                        for (int j = 0; j < Universe :: DIMENSION; j++)
                            q_pid[j] = &qvar.at(Universe :: get_pos(pid,j));

                        qcnst.push_back(q_pid);
                    }

                    c_eval[i] = cnst->get(qcnst);
                }

                return c_eval;
            }


            /**
             * @brief Takes gradient of f(q) with respect to r_j
             *  
             * @param j (see description)
             * @param f (see description)
             * @param qvar (see description)
             * @return std array of components of gradient
             */
            std :: array<autodiff :: var,3> grad_j(int j, autodiff :: var *f, std :: vector<var> &qvar){
                using namespace autodiff;

                //Can't loop here because of the way wrt() is written

                int xpos = Universe :: get_pos(j,0);
                int ypos = Universe :: get_pos(j,1);
                int zpos = Universe :: get_pos(j,2);

                var *x = &qvar.at(xpos);
                var *y = &qvar.at(ypos);
                var *z = &qvar.at(zpos);
                
                auto d = derivativesx(*f,wrt(*x,*y,*z));

                return d;
            }

            /**
             * @brief Sets M(i,Dimension * j + 0) = grad(0), M(i,Dimension * j + 1) = grad(1), M(i,Dimension * j + 1) = grad(2)
             *        
             * 
             * @param i (see description)
             * @param j (see description)
             * @param grad (see description)
             * @param M (see description)
             */
            void add_grad_to_mat(int i, int n, std :: array<autodiff :: var,3> &grad, Eigen :: MatrixXd *M){
                for (int k = 0; k < 3; k++){
                    auto kpos = Universe :: get_pos(n,k);
                    auto d_k = (double) grad.at(k);
                    (*M)(i, kpos) = d_k;
                }
            }


            /**
             * @brief Gets Jd based on the auxilary matricies K and v
             * 
             * @param K Auxilarry matricies of second partial derivatives of C
             * @param v Generalized velocity vector
             * @return Eigen Jd
             */
            Eigen :: MatrixXd get_Jd(std :: vector<Eigen :: MatrixXd> K, Eigen :: VectorXd v){
                using namespace Eigen;
                using namespace std;

                int cnst_count = K.size();
                int cols = v.size();

                MatrixXd Jd = MatrixXd :: Zero(cnst_count,cols);

                //The ith row of Jd_i is equal to 
                // (K_i . k)^T 
                for (int i = 0; i < cnst_count; i++){
                    auto kdotv = K.at(i)*v;
                    Jd.row(i) = kdotv.transpose();
                }

                return Jd;
            }


            
            //Eigen :: VectorXd solve(Eigen :: VectorXd const &F_ex, Universe &u){
//
            //}

    };

}

#endif