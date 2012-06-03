#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

class SEIR {
    public:
        SEIR( SEXP blocksize_, SEXP nstates_, SEXP nstep_, SEXP schooltype_):
          //spos(spos_), 
          blocksize(as<int>(blocksize_)), //maxday(as<int>(blocksize_)), 
          nstates(as<int>(nstates_)), day(0), repday(0), 
          nstep(as<int>(nstep_)), 
          schooltype(as<int>(schooltype_)), Pi(4*atan(1)) {
            state = arma::mat(nstates, blocksize );
            // initialize to 0 to avoid surprises...
            state.zeros();
            maxday = blocksize * nstep;
        }

    private:    // variable defs
        //arma::colvec ret;        // used?  
        arma::mat state;          // the state matrix - "final" output - grow as needed
        //NumericVector spos;         // start index of within-state compartments. needed??
        unsigned int blocksize;  // grow state by this much when needed
        unsigned int maxday;     // current ncol of state - must grow if more are needed
        Rcpp::List pars;         // List of model parameters -- access by name
        unsigned int nstates;     // number of states == number of rows in state matrix
        unsigned int day;                    // current day of the model
        unsigned int repday;                    // current reporting day
        double Pi;  //pi
        double rbirth, rsigma , rgamma , rdeltaR , rR0 , rpobs , rbetaforce, rbeta0, beta_now, deltat;
        NumericVector rschooldays;  //-1 for vacation, 1 for school
        int nstep;
        int schooltype;

    public:
        // model-dependent variables and methods
        // these need to be easily accessed by metapop
        // use doubles avoid confusion when multiplying by params
        double S, E, Eobs, I, Ieff, Iobs, R, N;  

        void setpars(SEXP pars_) {
            // this is essentially part of the model definition
            List tmp(clone(pars_));
            checkpars(tmp);
            pars = tmp;
            // params for basic SEIRN functionality
            nstep = as<double>(pars["nstep"]);
            deltat = as<double>(pars["deltat"]);
            rbirth = as<double>(pars["birth"]);
            rsigma = as<double>(pars["sigma"]);
            rgamma = as<double>(pars["gamma"]);
            rdeltaR = as<double>(pars["deltaR"]);
            rR0 = as<double>(pars["R0"]);
            rpobs = as<double>(pars["pobs"]);
            rbetaforce = as<double>(pars["betaforce"]);
            rbeta0 = as<double>(pars["beta0"]);
            rschooldays = pars["schooldays"];
            /*
            if (schooltype==1) {
            // See keeling et al 2001, "seasonally forced disease dynamics",
            // physica D, eqn 3.
            // R0 = (beta0/gamma)*(1+force)^((schooldays-holidays)/365)
            //    = beta_hat/gamma
            // and beta_hat = rR0*rgamma
                beta0 = (rR0*rgamma) / pow(1.0+rbetaforce, (2.0*rschooldays - 365.0)/365.0);
            }
            */
                
        }
        
        void setstate(SEXP init_) {
            // this is essentially part of the model definition
            // init state of current day
            NumericVector tmp(init_);
            checkstate(tmp);
            S = tmp(0);
            E = tmp(1);
            Eobs = tmp(2);
            I = tmp(3);
            Ieff = tmp(4);
            Iobs = tmp(5);
            R = tmp(6);
            N = tmp(7);
        }

    private:
        // model-dependent variables and methods
        void step(void) {
            // core model definition
            // modify state variables in-place, accumulate for reporting
            int ndeltaR;
            // births adjusted for infant mortality
            // deaths and migrations included in deltaR
            // no deaths from SEI
            if (schooltype == 0) {
                //sin forcing
                beta_now = 
                  rR0*rgamma*
                  (1.0-rbetaforce*
                    cos(2.0*Pi*(day)/365.0));
            };
            if (schooltype == 1) {
                int doy = day % 365;
                beta_now = rbeta0*pow(1.0+rbetaforce, rschooldays( doy ));
            };
            int nbirth = Rf_rpois( N*rbirth );
            int nlatent = Rf_rpois( (beta_now*S*Ieff)/(N));
            //int nlatent = Rf_rpois( (beta_now*S*Ieff)/(N+Ieff-I));
            int ninfect = Rf_rpois( rsigma*E );
            int nrecover = Rf_rpois( rgamma*I );
            // if deltaR is negative, need a double negative
            if ( rdeltaR < 0 ) { 
                rdeltaR *= -1;
                ndeltaR = Rf_rpois( N * rdeltaR); 
                ndeltaR *= -1;
            }
            else { ndeltaR = Rf_rpois( N *  rdeltaR) ; };
            // checks to prevent state < 0
            while ( S + nbirth < nlatent ) nlatent--;
            while ( E + nlatent < ninfect ) ninfect--;
            while ( I + ninfect < nrecover ) nrecover--;
            S +=  (nbirth - nlatent);
            E +=  (nlatent - ninfect);
            Eobs += nlatent;
            // instantaneous number infected
            I +=  (ninfect - nrecover);
            // prevalence -- all infected in this reporting period (nstep)
            Iobs +=  ninfect;
            // Effective I is computed at the metapop level for this timestep
            R += (nrecover + ndeltaR);
            N += (nbirth + ndeltaR);
            //Rf_PrintValue(wrap(day));
            if ( day % nstep == 1 ) {
                arma::colvec ret(nstates);
                ret.zeros();
                ret(0) = S;
                ret(1) = E;
                ret(2) = Eobs;
                ret(3) = I;
                ret(4) = Ieff;
                ret(5) = Rf_rbinom(Iobs, rpobs);
                ret(6) = R;
                ret(7) = N;
                //ret(4) = Rf_rbinom(I, rpobs);
                state.col(repday) = ret;
                // set observeds to zero
                Eobs = 0;
                Iobs = 0;
                repday++;
            }
            //Rf_PrintValue(wrap(day));
        };
        
    public:
        // model-independent public methods
        void rewind(int nday){
            // rewind the model this many days
            day -= nday;
        }

        SEXP steps(const unsigned int days) {
            // advance the model this many days
            RNGScope scope; // when to call this?
            BEGIN_RCPP  //check to make sure we don't fall off the end
            for ( unsigned int i =0; i<days; i++) {
                if ( repday +1 > maxday) {
                    throw std::range_error("took too many steps -- try a larger blocksize");
                }
                step();             // move forwards one day
                day++;
                //Rf_PrintValue(wrap(day));
            }
            END_RCPP
        }
    
        NumericMatrix get(void) {
            // get the full state matrix
            arma::mat tmp(arma::trans(state.cols(1,day)));
            return wrap(tmp);
        }
        
        NumericVector getday(unsigned int i) {
            // get the full state matrix
            arma::rowvec tmp(state.col(i));
            return wrap( tmp );
        }


        void fillstate( int fillval ) {
            // fill the state with int (typically 0?)
            state.fill(fillval);
        } 
        
        List report(void){
            List ret;
            ret["pars"] = pars;
            //ret["spos"] = spos;
            ret["day"] = day;
            ret["maxday"] = maxday;
            ret["blocksize"] = blocksize;
            NumericMatrix tmp = wrap(arma::trans(state));
            ret["state"] = tmp;
            return ret;
        }

    private:
        // model-independent private methods
        SEXP checkpars(List pars_) {
            // TODO -- add checking for names...
            BEGIN_RCPP
            if ( pars_.size() < 3 ) {
                throw std::range_error("pars list too small");
            };
            NumericVector dummy(1);
            return dummy;
            END_RCPP
        }
            
        SEXP checkstate(NumericVector init_) {
            // check that init is the right length
            BEGIN_RCPP
            if ( init_.size() != nstates )  {
                throw std::range_error("init vector has incorrect dimensions");
            };
            //if ( init_.min() < 0 ) {
                //throw std::range_error("init vector cannot contain negative values");
            //};
            // check that all are integers!!
            NumericVector dummy(1);
            return dummy;
            END_RCPP
        }


};

/*
RCPP_MODULE(seirmod){
	using namespace Rcpp ;
    class_<SEIR>("SEIR")
    
        //SEIR( SEXP blocksize_, SEXP pars_, SEXP init_, SEXP spos_): 
        //SEIR( SEXP blocksize_, SEXP nstates_): 
    .constructor<SEXP, SEXP, SEXP, SEXP>("args: blocksize (number of days to preallocate), pars (named list containing reals), init (integer vector of initial values), startposition (vector of within-state compartment start indices, currently not used).")
    .constructor<SEXP, SEXP>("args: blocksize (number of days to preallocate), nstates.")
    .method("steps", &SEIR::steps, "args: nday (int).  Advance model by nday steps" )
    .method("get", &SEIR::get, "args: none. Return full model state-matrix so far.")
    .method("getday", &SEIR::getday, "args: day (int).  Return model state vector for specified day." )
    .method("setpars", &SEIR::setpars, "args: pars (named list).  Post-initialization change of model parameters.")
    .method("setstate", &SEIR::setstate, "args: statevec (integer vector of state values). Change model state for current day.")
    .method("rewind", &SEIR::rewind, "args: nday (int).  Rewind model by this many days.")
    .method("report", &SEIR::report, "args: none.  Return named list containing all current variables.")
    ;
}        
*/
        

class Metapop {
    public:
        Metapop(SEXP npop_, SEXP xymat_, SEXP cmat_, SEXP blocksize, SEXP nstates_, SEXP nstep, SEXP schooltype) : 
          nstates(as<int>(nstates_)), npop(as<int>(npop_)), xymat(xymat_), couplemat(cmat_), thiscity(0)
        {
            SEIR tmpseir(blocksize, nstates_, nstep, schooltype);
            pops = std::vector<SEIR>(npop, tmpseir);
            //pops.resize(npop);
            //for (int ii = 0; ii<npop; ii++) {
            //  pops[ii].fillstate(0);
            //}
            // check that nrow xymat == npop!!
        }

        List cityreport(int n) {
            Rcpp::List ret = pops[n-1].report();
            ret["city"] = n-1;
            return ret; 
        }

        List report() {
            Rcpp::List ret;
            for (int ii=0; ii<npop; ii++) {
                ret.push_back( pops[ii].report());
            };
            return ret; 
        }
        
        SEXP setcity( int n_) {
            // just returns exception
            // check that n_ within bounds
            BEGIN_RCPP
            if ( n_-1 > npop-1 || n_-1 < 0 )  {
                throw std::range_error("index of city in setcity out of range");
            };
            thiscity = n_-1;
            END_RCPP
        }
        
        SEXP setstates( SEXP citystates) {
            NumericMatrix tmpmat(citystates);
            BEGIN_RCPP
            if ( (tmpmat.rows() != npop)  || (tmpmat.ncol() != nstates) ) {
                throw std::range_error("improper dimensions of state matrix -- npop rows and nstate cols");
            };
            for ( int ii = 0; ii < npop; ii++) {
               pops[ii].setstate( wrap(tmpmat(ii,_))) ;
            } 
            END_RCPP
        }


        SEXP setpars( SEXP pars) {
            List tmppars(pars);
            BEGIN_RCPP
            if ( tmppars.size() != npop)  {
                throw std::range_error(" the paramater list of lists should have npop elements");
            };
            for ( int ii = 0; ii < npop; ii++) {
               pops[ii].setpars( tmppars(ii)) ;
            } 
            END_RCPP
        }


        void steps( int n ) {
            for (int ii = 0; ii<n; ii++) {
                step();
            }
        }

    private:
        //SEIR &thiscity;
        unsigned int nstates;
        unsigned int npop;
        std::vector<SEIR> pops;
        NumericMatrix xymat;  //distance matrix between cities
        NumericMatrix couplemat;  //metapop coupling matrix between cities
        int thiscity;

        void step( ) {
            // do processing of migration here,
            // then take a step for each city
            unsigned int ithis, iother;
            // first, get and save effective I for each pop
            for ( ithis = 0; ithis < npop; ithis++) {
                //reset
                pops[ithis].Ieff = 0;
               for ( iother = 0; iother < npop; iother++) {
                    // couplemat isn't symmetric -- step through this city's row
                   pops[ithis].Ieff += couplemat(ithis, iother) * pops[iother].I;
                }
            }
            // then take the next step for all pops
            for ( ithis = 0; ithis < npop; ithis++) {
               pops[ithis].steps( 1 );
            } 
        }
};

RCPP_MODULE(seirmod){
	using namespace Rcpp ;
    class_<Metapop>("Metapop")
    
    .constructor<SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP>("args: npop (int, number of cities), xymat (npop*2 matrix of city locations), couplemat (between-city coupling matrix), blocksize (number of days to preallocate), nstates (number of state variables), nstep (number of days per observation), schooltype (int, sin=0, force=1")
    /*
    .constructor<SEXP, SEXP>("args: blocksize (number of days to preallocate), nstates.")
    .method("steps", &SEIR::steps, "args: nday (int).  Advance model by nday steps" )
    .method("get", &SEIR::get, "args: none. Return full model state-matrix so far.")
    .method("getday", &SEIR::getday, "args: day (int).  Return model state vector for specified day." )
    .method("setpars", &SEIR::setpars, "args: pars (named list).  Post-initialization change of model parameters.")
    .method("setstate", &SEIR::setstate, "args: statevec (integer vector of state values). Change model state for current day.")
    .method("rewind", &SEIR::rewind, "args: nday (int).  Rewind model by this many days.")
    */
    .method("report", &Metapop::report, "args: none.  Return ordered list of named lists, each containing all current variables for 1 city.")
    .method("cityreport", &Metapop::cityreport, "args: none.  Return named list containing all current variables. for this city")
    .method("setcity", &Metapop::setcity, "args: int citynumber")
    .method("setstates", &Metapop::setstates, "args: npopXnstate matrix of states for current timestep")
    .method("setpars", &Metapop::setpars, "args: list-of-list of length npop, each containing that city's parameters")
    .method("steps", &Metapop::steps, "args: number of steps to advance all cities")
    ;
}        
        

