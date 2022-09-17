#ifndef MATH_UTIL_H // include guard
#define MATH_UTIL_H

#include <cmath>
#include <vector>
#include <map>

using namespace std;

/**
 * @brief Namespace that contains tools for dealing with series.s
 * 
 */
namespace math_series{

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
            map<int,long double> coef;
            
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
            TaylorSeries(map<int,long double> _coef){
                //Is this the right way
                coef = _coef;
            };

            /**
             * @brief Construct a new Taylor Series object
             * 
             * @param _coef The coeficients of the taylor series. Coeficients are in the form (n, a_n)
             */
            TaylorSeries(map<int,long double> _coef,long double _x0){
                //Is this the right way
                coef = _coef;
                x0 = _x0;
            };

            /**
             * @brief Evaluates the taylor series at a specific point
             * 
             * @param x The point to evaluate the taylor series at
             * @return long double : The evaluation of the taylor series at the given point
             */
            long double evaluate(long double x){

                long double out = 0;

                //Evaluate \sum a_n * x ^ n
                for (auto itr = coef.begin(); itr != coef.end(); ++itr){
                    int n = itr->first;
                    int a_n = itr->second;
                    out += a_n * pow(x - x0,n);
                };

                return out;
            };
    };


}

#endif