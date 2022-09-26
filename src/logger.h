#ifndef LOGGER_H// include guard
#define LOGGER_H

#include <vector>
#include <list>
#include "universe.h"

namespace physics{

    class Logger{
        

        public:
            std :: vector<std :: string> colnames;

            virtual std :: vector<std :: vector<double>> get_dat() = 0;

            virtual void log(Universe &u) = 0;

            std :: string get_csv_output(){
                using namespace std;

                auto dat = this->get_dat();

                string out = "";

                int colcount = colnames.size();

                for (int i = 0; i < colcount; i++){
                    string append = ",";
                    if (i == colcount - 1)
                        append = "\n";

                    out += colnames.at(i) + append;
                }

                for (auto v : dat){

                    auto vsize = v.size();
                    assert(vsize == colcount);

                    for (int i = 0; i < vsize; i++){
                        string append = ",";
                        if (i == colcount - 1)
                            append = "\n";
                        out += to_string(v.at(i)) + append;
                    }
                }

                return out;
            }

            static void add_evtostdv(EVector &add, std :: vector<double> &addto){
                for (int i = 0; i < add.size(); i++){
                    addto.push_back(add(i));
                }
            }

            int get_colcount(){
                return colnames.size();
            }
    };

    class DefaultLogger : public Logger{

        public:
            const std :: string XNAME = "X";
            const std :: string YNAME = "Y";
            const std :: string ZNAME = "Z";


            const std :: array<std :: string,3> COORDNAMES = {XNAME,YNAME,ZNAME};


            DefaultLogger(Universe &u){
                init_colnames(u);
            }

            void log(Universe &u) override{
                qdat.push_back(u.get_q());
                vdat.push_back(u.get_v());
                tdat.push_back(u.clock);
            }

            std :: vector<std :: vector<double>> get_dat() override{
                using namespace std;
                
                auto datcount = qdat.size();

                vector<vector<double>> out;

                auto itrq = qdat.begin();
                auto itrv = vdat.begin();
                auto itrt = tdat.begin();

                int colcount = get_colcount();

                for (int i = 0; i < datcount; i++){
                    
                    vector<double> entry;

                    auto t = *itrt;

                    int size = itrq->size() + itrv->size() + 1;
                    assert(size == colcount);
                    entry.push_back(*itrt);
                    add_evtostdv(*itrq,entry);
                    add_evtostdv(*itrv,entry);

                    auto sz = entry.size();

                    out.push_back(entry);

                    itrq++;
                    itrv++;
                    itrt++;
                }

                return out;
            }

        private: 
            void init_colnames(Universe &u){
                using namespace std;
                
                int pcount = u.get_p_count();


                int qlen = UDim * pcount;

                colnames.resize(2 * qlen + 1);

                colnames.at(0) = "t";

                for (int i = 0; i < u.get_p_count();i++){
                    for (int j = 0; j < UDim; j++){
                        string nameq = COORDNAMES.at(j) + to_string(i);
                        string namev = COORDNAMES.at(j) + "'" + to_string(i);
                        
                        int pos = Universe :: get_pos(i,j) + 1;

                        colnames.at(pos) = nameq;
                        colnames.at(qlen + pos) = namev;
                    }
                }
            }

            std :: list<EVector> qdat;
            std :: list<EVector> vdat;
            std :: list<double> tdat; 
    };


}


#endif