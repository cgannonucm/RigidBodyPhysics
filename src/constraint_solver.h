#ifndef CONSTRAINT_SOLVER_H // include guard
#define CONSTRAINT_SOLVER_H

#include "universe.h"
#include "constraint.h"
#include <Eigen/Eigen/Dense>
#include <Eigen/Eigen/SparseCore>
#include <Eigen/Eigen/SparseCholesky>

namespace physics{

    class ConstraintSolver{ 
        public:
        


        /**
             * @brief Solves for the constraint forces acting on the particles in the universe
             * 
             * @param F_ex The external forces acting on particles
             * @param u The Universe
             * @return Eigen Constraint forces acting on particles
             */
            EVector solve(const EVector &q, const EVector &v,const EVector &m,
            vConstraints &cnst, const EVector &F_ex){
                using namespace Eigen;
                //TODO implement options for solving
                //Solve J * M ^-1 * J^T * x = -Jd * v - J * M^-1 * F_ex
                //F_c = J^T * x
                //Convert to the equation Ax = B
                //A = J * (M ^ -1) * J^T 
                //B = -(Jd * v) - ((M^-1) * F_ex)

                //Prepare to solve
                auto [J, Jd] = get_JJd(q,v,cnst);

                Mat Jt = J.transpose();
                
                Mat M_inv = get_M_inv(m);

                Mat A = J * M_inv * Jt;

                VectorXd B = -(Jd * v) - (J * M_inv * F_ex);

                //Solve
                #if !SPARSESOLVE
                VectorXd x = A.colPivHouseholderQr().solve(B);
                #endif

                #if SPARSESOLVE
                Eigen :: SimplicialLLT<Eigen :: SparseMatrix<double>> solver;
                
                auto S = A.sparseView();
                            solver.compute(S);
                                
                EVector x = solver.solve(B);
                #endif

                //Convert to constraint force
                VectorXd F_c = Jt * x;
                
                return F_c;

            }

            /**
             * @brief Get the MInv matrix size (3n, 3n) from the m vector (m1,m2,...,mn)
             * 
             * @param u Universe to get 
             * @return Eigen 
             */
            static Mat get_M_inv(EVector mv){
                using namespace Eigen;

                int n = mv.size();
                int dim = UDim * n;

                Mat MInv = Mat :: Zero(dim,dim);

                //M is diagonal matrix of whose diagonal is (m1,m1,m1,...,mn,mn,mn)
                //Therefore M inv is a matrix whose diagonal is (1/m1,1/m1,1/m1,...,1/mn,1/mn,1/mn)
                for (int i = 0; i < n; i++){
                    double m = mv(i);

                    for (int k = 0; k < UDim; k++){
                        int j = Universe :: get_pos(i,k);
                        MInv(j,j) = 1/m;
                    }
                }

                return MInv;
            }
        
        private:
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
            std :: array<Mat,2> get_JJd(EVector q, EVector v, std :: vector<Constraint *> constraints){
                USING_STANDARD_NAMESPACES;
                using namespace autodiff;
                
                int cnst_count = constraints.size();
                int rows = cnst_count;
                int cols = q.size();
                    

                Mat J = Mat :: Zero(rows,cols);
                
                //Auxilary matricies H_j, used to calculate Jd
                //H is the Hessian Tensor of C
                vector<Mat> H = init_Kmats(rows,cols);

                vector<var> qvar = get_qvar(q);

                vector<var> c_eval = eval_constraints(constraints,qvar);

                //Calculate J matrix and H matricies
                //J = (partial_j(C_i))_ij
                //H_j = (partial_i(partial_k(C_j)))_ik
                for (int i = 0; i < cnst_count; i++){
                    var *c = &c_eval.at(i);

                    auto p_ids = constraints.at(i)->get_pids();
                    
                    //Only need to take take derivatives with respect to the coordinates passed to C_i
                    for (int n : p_ids){
                        
                        auto grad = grad_jx(n,c,qvar);
                        add_grad_to_mat(i,n,grad,&J);

                        //Calculate H matrix
                        //We need to calculate second partial derivatives
                        //Todo since partial_i(partial_j(C)) =  partial_j(partial_i(C))
                        //IE H is symetric
                        //we can optimize by elimating these redundant calculations
                        

                        for (int m = 0; m < UDim; m++){
                            for (int n2 : p_ids){
                                int j = UDim * n + m;

                                //Ignore redundant calculations
                                if(n2 > n)
                                    break;

                                auto grad2 = grad_j(n2,&grad.at(m),qvar);

                                for(int l = 0; l < UDim; l++){
                                    int k = Universe :: get_pos(n2,l);
                                    //Using the fact H_ijk = H_ikj
                                    auto grad2_l = grad2.at(l);

                                    H.at(i)(j,k) = grad2_l;
                                    H.at(i)(k,j) = grad2_l; 
                                }
                            }                 
                        }
                    }
                }

                Mat Jd = get_Jd(H,v);

                //cout << J << endl;
                //cout << Jd << endl; 

                return {J,Jd};
            }


            /**
             * @brief Initializes mcount number of dim * dim matricies filled with 0s
             * 
             * @param mcount (see description)
             * @param dim (see description)
             * @return std (see description)
             */
            std :: vector<Mat> init_Kmats(int mcount, int dim){
                USING_STANDARD_NAMESPACES;

                vector<Mat> K;
                K.resize(mcount);
                 
                for (int i = 0; i < mcount; i++){
                    K[i] = Mat :: Zero(dim, dim);
                }

                return K;
            }

            /**
             * @brief Convert from VectorXd q to vector<var> qvar to prepare for diferentiation
             * 
             * @param q Generalized position variable
             * @return std (see description)
             */
            std :: vector<autodiff :: var> get_qvar(EVector &q){
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
                USING_STANDARD_NAMESPACES;
                using namespace autodiff;
                
                auto cnst_count = constraints.size();

                vector<var> c_eval;
                c_eval.resize(cnst_count);

                for (int i = 0; i < cnst_count; i++){
                    Constraint *cnst = constraints.at(i);

                    //Only pass requested coordintates 
                    vector<vector<var *>> qcnst;
                    for (int pid : cnst->get_pids()){
                        
                        vector<var *> q_pid;
                        q_pid.resize(UDim);

                        for (int j = 0; j < UDim; j++)
                            q_pid[j] = &qvar.at(Universe :: get_pos(pid,j));

                        qcnst.push_back(q_pid);
                    }

                    c_eval[i] = cnst->get(qcnst);
                }

                return c_eval;
            }


            /**
             * @brief Takes gradient of f(q) with respect to r_j, returns a variable that can be diferentiated again
             *  
             * @param j (see description)
             * @param f (see description)
             * @param qvar (see description)
             * @return std array of components of gradient
             */
            std :: array<autodiff :: var,3> grad_jx(int j, autodiff :: var *f, std :: vector<var> &qvar){
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
             * @brief Takes gradient of f(q) with respect to r_j. Returns a double that cannot be differentiated again
             *  
             * @param j (see description)
             * @param f (see description)
             * @param qvar (see description)
             * @return std array of components of gradient
             */
            std :: array<double,3> grad_j(int j, autodiff :: var *f, std :: vector<var> &qvar){
                using namespace autodiff;

                //Can't loop here because of the way wrt() is written

                int xpos = Universe :: get_pos(j,0);
                int ypos = Universe :: get_pos(j,1);
                int zpos = Universe :: get_pos(j,2);

                var *x = &qvar.at(xpos);
                var *y = &qvar.at(ypos);
                var *z = &qvar.at(zpos);
                
                auto d = derivatives(*f,wrt(*x,*y,*z));

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
            void add_grad_to_mat(int i, int n, std :: array<autodiff :: var,3> &grad, Mat *M){
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
            Mat get_Jd(std :: vector<Mat> &K, EVector &v){
                USING_STANDARD_NAMESPACES;

                int cnst_count = K.size();
                int cols = v.size();

                Mat Jd = Mat :: Zero(cnst_count,cols);

                //The ith row of Jd_i is equal to 
                // (K_i . k)^T 
                for (int i = 0; i < cnst_count; i++){
                    auto kdotv = K.at(i)*v;
                    Jd.row(i) = kdotv.transpose();
                }

                return Jd;
            }





            

    };

}

#endif