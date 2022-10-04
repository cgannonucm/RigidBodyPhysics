#ifndef MATH_UTIL_H // include guard
#define MATH_UTIL_H

#include <cmath>
#include <vector>
#include <map>

/**
 * @brief Namespace that contains tools for dealing with series.s
 * 
 */
namespace physics{

    namespace math_util{

        /**
         * @brief Abstract class representing a series.
         * 
         */
        class Series{
            public:
                virtual long double evaluate(long double x) = 0;

        };


        /**
         * @brief A class representing a Taylor series. 
         * A taylor series is in the form \sum a_n * x ^ n
         * 
         */
        class TaylorSeries : public Series{
            public:

                /**
                 * @brief Coeficients of the taylor series the coeficients are in the form (n, a_n)
                 * 
                 */
                std :: map<int,long double> coef;
                
                /**
                 * @brief Offset of the series
                 * 
                 */
                long double x0 = 0;

                /**
                 * @brief Construct a new Taylor Series object
                 * 
                 * @param _coef The coeficients of the taylor series. 
                 */
                TaylorSeries(std :: map<int,long double> _coef);

                /**
                 * @brief Construct a new Taylor Series object
                 * 
                 * @param _coef The coeficients of the taylor series. Coeficients are in the form (n, a_n)
                 */
                TaylorSeries(std :: map<int,long double> _coef,long double _x0);

                /**
                 * @brief Evaluates the taylor series at a specific point
                 * 
                 * @param x The point to evaluate the taylor series at
                 * @return long double : The evaluation of the taylor series at the given point
                 */
                long double evaluate(long double x);
        };


    }
}

#endif