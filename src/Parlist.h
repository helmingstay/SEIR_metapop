using namespace Rcpp;

template <typename T>
class Parlist {
    public:
        //Parlist():
            // check that list is initialized before accessing??
            //{}
        T operator()(const std::string varname){
            // get parameter list element named varname, return as int
            return as<T>(list[varname]);
        };

        void add(std::string varname, T val){
            // use the operator() to do the cast when getting the value
            list[varname] = (*this)(varname) + val;
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

        void init_fromnames( CharacterVector nameslist) {
            // init from a list of names, fill with zeros (as doubles)
            // called from C++ Events
            names = mk_names(nameslist);
            list = mk_list(names);
            N = list.size();
            // fill with zeros 
            fill(0);
        }
        void fill( T val) {
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

        T min() {
            std::vector<T> tmpvec;
            tmpvec = std::vector<T>(N);
            std::copy(list.begin(), list.end(), tmpvec.begin());
            T min_elem = *(std::min_element(tmpvec.begin(), tmpvec.end()));
            return min_elem;
        }

        arma::colvec get_colvec() {
            arma::colvec ret(N);
            for(int ii=0; ii<N; ii++){
                ret[ii] = static_cast<T>(list[ii]);
            }
            return ret;
        }

    public:
        // variables
        int N;
        std::vector< std::string> names;
        Rcpp::List list;
    //map<std::string, T> list;

    private:

        std::vector< std::string> mk_names( CharacterVector names_list ) {
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



