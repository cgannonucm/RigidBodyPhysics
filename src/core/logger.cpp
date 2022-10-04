#include "logger.hpp"

namespace physics{

    /**
     * Logger
     * ------------------------------------
    */

    std :: string Logger::get_csv_output(){
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
    
    void Logger::add_evtostdv(EVector &add, std :: vector<double> &addto){
        for (int i = 0; i < add.size(); i++){
            addto.push_back(add(i));
        }
    }

    int Logger::get_colcount(){
        return colnames.size();
    }

    /**
     * ------------------------------------
    */
   /**
    * Default Logger
    * -------------------------------------
    * TODO clean up ugly code
   */
    
    DefaultLogger::DefaultLogger(Universe &u){
        init_colnames(u);
    }

    DefaultLogger::DefaultLogger(Universe &u, double _tstep){
        init_colnames(u);
        timestep = _tstep;
    }

    void DefaultLogger::log(Universe &u){
        //Only record at the requested time intervals
        if(u.clock < nextstep)
            return;
        qdat.push_back(u.get_q());
        vdat.push_back(u.get_v());
        tdat.push_back(u.clock);
        
        nextstep += timestep;
    }

    std :: vector<std :: vector<double>> DefaultLogger::get_dat(){
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

    void DefaultLogger::init_colnames(Universe &u){
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
    /**
     * -------------------------------------
    */
    

}