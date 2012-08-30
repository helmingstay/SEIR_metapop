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
        // why does this not work??
        //NumericVector NV(std::string varname){
            //SEXP tmp = pairlist[varname];
            //NumericVector ret(tmp);
            //return ret;
        //};
        void set(SEXP newlist){
            Rcpp::List tmppar(newlist);
            parlist = tmppar;
        }

    private:
        Rcpp::List parlist;
};
