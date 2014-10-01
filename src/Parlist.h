using namespace Rcpp;

template <typename T>
class Parlist {
    typedef typename std::vector< std::string> STRINGVEC;
    typedef typename std::map<std::string, T> MAP;
    typedef typename std::map<std::string, T>::iterator MAPIT;

    public:
        //Parlist():
            // check that list is initialized before accessing??
            //{}
        //?? reference??
        inline T& operator()(const std::string &varname, bool debug = false){
            // give error if varname not in names?
            //bool debug = false;
            if (debug) {
                STRINGVEC::iterator it = find(names.begin(), names.end(), varname);
                if (it == names.end()) {
                    Rf_PrintValue(wrap("Attempt to access non-existent named element:"));
                    Rf_PrintValue(wrap(varname));
                    throw std::range_error(varname);
                }
            }
            // get parameter list element named varname, return as int
            return list[varname];
        };

        // !! remove this, replace with () in Pop.h
        inline T& operator[](const std::string &varname){
            // get parameter list element named varname, return as int
            return list[varname];
        };

        inline void add(const std::string &varname, T val){
            // use the operator() to do the cast when getting the value
            list[varname] = list[varname] + val;
        };

        //Rcpp::List& operator()( std::string name ){
            //return list[name];
        //}

        void set(SEXP newlist){
            std::size_t ii;
            // initialize elements from a named input list (passed as SEXP) with values
            Rcpp::List tmplist(newlist);
            CharacterVector tmpnames = tmplist.attr("names") ;
            std::string tmpname;
            // reinitialize list with zero-fill
            init_fromnames(tmpnames);
            // then fill with new values by name
            for (ii = 0; ii<N; ii++) {
               tmpname = names[ii];
                list[ tmpname ] = as<T>( tmplist[ tmpname ] );
            }
        }

        void set(Rcpp::List tmplist){
            std::size_t ii;
            // initialize elements from a named input list (passed as SEXP) with values
            CharacterVector tmpnames = tmplist.attr("names") ;
            std::string tmpname;
            // reinitialize list with zero-fill
            init_fromnames(tmpnames);
            // then fill with new values by name
            for (ii = 0; ii<N; ii++) {
               tmpname = names[ii];
                list[ tmpname ] = as<T>( tmplist[ tmpname ] );
            }
        }

        void init_fromnames( CharacterVector names_list) {
            std::size_t ii;
            // Generate zero-filled named list from vector of names
            names = Rcpp::as< STRINGVEC >(names_list);
            N = names_list.size();
            if (N != names.size()) {
                Rf_PrintValue(wrap(list));
                throw std::range_error("Error in Parlist.set: newlist must have names");
            }
            for(ii = 0; ii<N; ii++){
              list[ names[ii] ] = static_cast<T>(0);
            }
        }

        inline void fill( T val) {
            MAPIT it;
            // fill list with a single value
            for (it = list.begin(); it != list.end(); it++){
                it->second = val;
            }
            //std::fill( list.begin(), list.end(), val);
        }

        void copy( const std::vector<T> &vals ) {
            std::size_t ii;
            // copy vector of values into the list by position
            if (list.size() != vals.size()) {
                Rf_PrintValue(wrap(list));
                Rf_PrintValue(wrap(vals));
                throw std::range_error("Error in Parlist.copy: new values must have same dim as list");
            }
            //std::copy(vals.begin(), vals.end(), list.begin());
            // map iterator doesn't seem to order by first fill...
            // fill by name instead
            for (ii = 0; ii < N; ii++) {
                list[ names[ii] ] = vals[ii];
            }
        }

        arma::Col<T> get_colvec() {
            std::size_t ii;
            arma::Col<T> ret(N);
            // copy list values into the column for return
            for (ii = 0; ii < N; ii++) {
                ret[ii] = list[ names[ii] ];
            }
            return ret;
        }

    public:
        // variables
        std::size_t N;
        STRINGVEC names;
        MAP list;

};



