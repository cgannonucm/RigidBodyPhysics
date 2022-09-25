#ifndef PHYSICS_MACROS_H // include guard
#define PHYSICS_MACROS_H

#include <Eigen/Eigen/Dense>

#define EVector Eigen :: VectorXd
#define EVectorNd Eigen :: Vector3d
#define UDim 3
typedef std :: function<EVector (EVector &, EVector &)> funcqv;


#endif


