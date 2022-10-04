#ifndef LOGGER_H// include guard
#define LOGGER_H

#include <vector>
#include <list>
#include "universe.hpp"

namespace physics{

    class Logger{
        

        public:
            std :: vector<std :: string> colnames;

            virtual std :: vector<std :: vector<double>> get_dat() = 0;

            virtual void log(Universe &u) = 0;

            std :: string get_csv_output();

            static void add_evtostdv(EVector &add, std :: vector<double> &addto);

            int get_colcount();
    };

    //Come back and fix ugly code
    class DefaultLogger : public Logger{

        public:
            const std :: string XNAME = "X";
            const std :: string YNAME = "Y";
            const std :: string ZNAME = "Z";
            
            const std :: array<std :: string,3> COORDNAMES = {XNAME,YNAME,ZNAME};

            DefaultLogger(Universe &u);

            DefaultLogger(Universe &u, double _tstep);

            void log(Universe &u) override;

            std :: vector<std :: vector<double>> get_dat() override;

        private: 
            void init_colnames(Universe &u);

            double nextstep = 0;
            double timestep = 0;

            std :: list<EVector> qdat;
            std :: list<EVector> vdat;
            std :: list<double> tdat; 
    };


}


#endif