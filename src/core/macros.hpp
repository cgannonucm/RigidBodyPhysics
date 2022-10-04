#ifndef PHYSICS_MACROS_H // include guard
#define PHYSICS_MACROS_H

#include <Eigen/Eigen/Dense>

#define UDim 3
#define USING_STANDARD_NAMESPACES using namespace Eigen; using namespace std
#define SPARSESOLVE false
typedef Eigen :: MatrixXd Mat;


#if UDim == 3
    typedef Eigen :: VectorXd EVector;
    typedef Eigen :: Vector3d EVectorNd;
    typedef std :: function<EVector (EVector &, EVector &)> funcqv;
#endif


#endif


