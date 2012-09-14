using namespace Rcpp;

class Parlist {
    public:
        //Parlist():
            // check that list is initialized before accessing??
            //{}
        int i(std::string varname){
            // get parameter list element named varname, return as int
            return as<int>(list[varname]);
        };
        double d(std::string varname){
            return as<double>(list[varname]);
        };

        double add(std::string varname, double val){
            list[varname] = list[varname] + val;
        };

        //Rcpp::List& operator()( std::string name ){
            //return list[name];
        //}

        void set(Rcpp::List newlist){
            // initialize elements from a named input list with values
            list = newlist;
            N = list.size();
            names = mk_names( list.attr("names") );
            if (N != names.size()) {
                Rf_PrintValue(list);
                throw std::range_error("Error in Parlist.set: newlist must have names");
            }
        }

        void set(SEXP newlist_){
            // exposed to R
            //
            Rcpp::List newlist(newlist_);
            set(newlist);
        }

        void init_fromnames( Rcpp::List nameslist) {
            // needed??
            // init from a list of names, fill with zeros (as doubles)
            // called from C++ Events
            names = mk_names(nameslist);
            list = mk_list(names);
            N = list.size();
            // fill with zeros 
            fill(0);
        }

        void fill( double val) {
            std::fill( list.begin(), list.end(), val);
        }

        void copy( NumericVector vals ) {
            if (list.size() != vals.size()) {
                Rf_PrintValue(list);
                Rf_PrintValue(vals);
                throw std::range_error("Error in Parlist.copy: new values must have same dim as list");
            }
            std::copy(vals.begin(), vals.end(), list.begin());
        }

    public:
        // variables
        int N;
        std::vector< std::string> names;
        Rcpp::List list;

    private:
        // convenience functions
        std::vector< std::string> mk_names( Rcpp::List names_list ) {
            // function to turn list of names (from attributes) into vector suitable for list indexing
            std::vector< std::string> ret = Rcpp::as< std::vector< std::string> >(names_list);
            return ret;
        }
        Rcpp::List mk_list( std::vector< std::string> list_names ) {
            // function to turn named list from vector of names
            // initialize list by size
            Rcpp::List ret( list_names.size() );
            // set names
            ret.attr("names") = list_names;
            return ret;
        }

};



