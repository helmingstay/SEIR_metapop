using namespace Rcpp;

class Parlist {
    public:
        //Parlist():
            // check that list is initialized before accessing??
            //{}
        int i(std::string varname){
            // get parameter list element named varname, return as int
            return as<int>(parlist[varname]);
        };
        double d(std::string varname){
            return as<double>(parlist[varname]);
        };

        void set(Rcpp::List newlist){
            // initialize elements from a named input list with values
            parlist = newlist;
            N = parlist.size();
            names = mk_names( parlist.attr("names") );
            if (N != names.size()) {
                Rf_PrintValue(parlist);
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
            parlist = mk_list(names);
            N = parlist.size()
            // fill with zeros 
            fill(0);
        }

        void fill( double val) {
            fill( parlist.begin(), parlist.end(), val);
        }

        void copy( NumericVector vals ) {
            if (parlist.size() != vals.size()) {
                Rf_PrintValue(parlist);
                Rf_PrintValue(vals);
                throw std::range_error("Error in Parlist.copy: new values must have same dim as parlist");
            }
            std::copy(vals.begin(), vals.end(), parlist.begin());
        }

    public:
        // variables
        int N;
        std::vector< std::string> names;
        Rcpp::List parlist;

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



