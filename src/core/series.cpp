#include "series.hpp"


/**
 * @brief Namespace that contains tools for dealing with series.s
 * 
 */
namespace physics{

    namespace math_util{

        /**
         * @brief Construct a new Taylor Series object
         * 
         * @param _coef The coeficients of the taylor series. 
         */
        TaylorSeries::TaylorSeries(std :: map<int,long double> _coef){
            //Is this the right way
            coef = _coef;
        }


        /**
         * @brief Construct a new Taylor Series object
         * 
         * @param _coef The coeficients of the taylor series. Coeficients are in the form (n, a_n)
         */
        TaylorSeries::TaylorSeries(std :: map<int,long double> _coef,long double _x0){
            //Is this the right way
            coef = _coef;
            x0 = _x0;
        }

        /**
         * @brief Evaluates the taylor series at a specific point
         * 
         * @param x The point to evaluate the taylor series at
         * @return long double : The evaluation of the taylor series at the given point
         */
        long double TaylorSeries::evaluate(long double x){
            long double out = 0;
            //Evaluate \sum a_n * x ^ n
            for (auto itr = coef.begin(); itr != coef.end(); ++itr){
                int n = itr->first;
                int a_n = itr->second;
                out += a_n * pow(x - x0,n);
            };
            return out;
        }

    }
}