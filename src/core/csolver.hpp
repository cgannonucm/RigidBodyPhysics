#ifndef CONSTRAINT_SOLVER_H // include guard
#define CONSTRAINT_SOLVER_H

#include "universe.hpp"
#include "constraint.hpp"
#include <Eigen/Eigen/Dense>
#include <Eigen/Eigen/SparseCore>
#include <Eigen/Eigen/SparseCholesky>

namespace physics{

    class CSolver{ 
        public:
            /**
             * @brief Solves for the constraint forces acting on the particles in the universe
             * 
             * @param F_ex The external forces acting on particles
             * @param u The Universe
             * @return Eigen Constraint forces acting on particles
             */
            EVector solve(const EVector &q, const EVector &v,const EVector &m,
            vConstraints &cnst, const EVector &F_ex);

            /**
             * @brief Get the MInv matrix size (3n, 3n) from the m vector (m1,m2,...,mn)
             * 
             * @param u Universe to get 
             * @return Eigen 
             */
            static Mat get_M_inv(EVector mv);
        
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
            std :: array<Mat,2> get_JJd(EVector q, EVector v, std :: vector<Constraint *> constraints);


            /**
             * @brief Initializes mcount number of dim * dim matricies filled with 0s
             * 
             * @param mcount (see description)
             * @param dim (see description)
             * @return std (see description)
             */
            std :: vector<Mat> init_Kmats(int mcount, int dim);

            /**
             * @brief Convert from VectorXd q to vector<var> qvar to prepare for diferentiation
             * 
             * @param q Generalized position variable
             * @return std (see description)
             */
            std :: vector<autodiff :: var> get_qvar(EVector &q);


            /**
             * @brief Evaluates the constraints in preperation for diferentiation
             * 
             * @param constraints The constraints to evaluate
             * @param qvar The q vector to evaluate constraints at
             * @return std The evaluated constraints
             */
            std :: vector<autodiff :: var> eval_constraints(std :: vector<Constraint *> &constraints, std :: vector<var> &qvar);


            /**
             * @brief Takes gradient of f(q) with respect to r_j, returns a variable that can be diferentiated again
             *  
             * @param j (see description)
             * @param f (see description)
             * @param qvar (see description)
             * @return std array of components of gradient
             */
            std :: array<autodiff :: var,3> grad_jx(int j, autodiff :: var *f, std :: vector<var> &qvar);

            /**
             * @brief Takes gradient of f(q) with respect to r_j. Returns a double that cannot be differentiated again
             *  
             * @param j (see description)
             * @param f (see description)
             * @param qvar (see description)
             * @return std array of components of gradient
             */
            std :: array<double,3> grad_j(int j, autodiff :: var *f, std :: vector<var> &qvar);

            /**
             * @brief Sets M(i,Dimension * j + 0) = grad(0), M(i,Dimension * j + 1) = grad(1), M(i,Dimension * j + 1) = grad(2)
             *        
             * 
             * @param i (see description)
             * @param j (see description)
             * @param grad (see description)
             * @param M (see description)
             */
            void add_grad_to_mat(int i, int n, std :: array<autodiff :: var,3> &grad, Mat *M);


            /**
             * @brief Gets Jd based on the auxilary matricies K and v
             * 
             * @param K Auxilarry matricies of second partial derivatives of C
             * @param v Generalized velocity vector
             * @return Eigen Jd
             */
            Mat get_Jd(std :: vector<Mat> &K, EVector &v);





            

    };

}

#endif