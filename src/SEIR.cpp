#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

class SEIR {
    public:
        SEIR( SEXP blocksize_, SEXP nstates_, SEXP nstep_):
          //spos(spos_), 
          blocksize(as<int>(blocksize_)), //maxday(as<int>(blocksize_)), 
          nstates(as<unsigned int>(nstates_)), day(0), repday(0), 
          nstep(as<int>(nstep_)), 
          Pi(4*atan(1)), 
          minparlen(14) {
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
        int nstep;
        double Pi;  //pi
        // this should get changed when we add states
        unsigned int minparlen;
        int schooltype, rschoollag, rimportmethod;
        double rbirth, rsigma , rgamma , rdeltaR , rR0 , rpobs , rbetaforce, rbeta0, rtheta, rimports, rSresid, rSfrac, beta_now, deltat, rbetasd_ratio;
        NumericVector rschooldays;  //-1 for vacation, 1 for school

    public:
        // model-dependent variables and methods
        // these need to be easily accessed by metapop
        // use doubles avoid confusion when multiplying by params
        double S, E, Eobs, I, Ieff, Iobs, R, N, Shidden;

        void setpars(SEXP pars_) {
            // this is essentially part of the model definition
            List tmp(clone(pars_));
            // change this if adding pars!!
            checkpars(tmp);
            pars = tmp;
            // params for basic SEIRN functionality
            // update minparlen if adding params
            // ??
            deltat = as<double>(pars["deltat"]);
            // percap per day birth rate 
            rbirth = as<double>(pars["birth"]);
            rsigma = as<double>(pars["sigma"]);
            rgamma = as<double>(pars["gamma"]);
            // percap per day change in R
            // includes death and migration
            rdeltaR = as<double>(pars["deltaR"]);
            // R0 used to calculate mean beta if schooltype == 0 (sin) !!fixme??
            rR0 = as<double>(pars["R0"]);
            // binomial probability of observing accumulated I (over nsteps)
            rpobs = as<double>(pars["pobs"]);
            // mean beta, modulated by schooltime
            // only used if schooltype == 1 (termtime) !!fixme??
            rbeta0 = as<double>(pars["beta0"]);
            // modulate beta0 by this much
            rbetaforce = as<double>(pars["betaforce"]);
            // add daily normal forcing to beta with SD as a ratio of current value of beta
            rbetasd_ratio = as<double>(pars["betasd_ratio"]);
            // modulates Reff
            rtheta = as<double>(pars["theta"]);
            // equilib proportion of susceptibles in hidden S class
            // set to 0 to disable hidden S
            rSfrac = as<double>(pars["Sfrac"]);
            // equilib residency time of susceptibles in hidden S class
            rSresid = as<double>(pars["Sresid"]);
            // 0 is sin, 1 is termtime
            schooltype = as<int>(pars["schooltype"]);
            // vector, only used if schooltype == 1
            rschooldays = pars["schooldays"];
            // lag of school term or sin
            rschoollag = as<int>(pars["schoollag"]);
            // rate of imports, meaning depends on importmethdo
            rimports = as<double>(pars["imports"]);
            // 1-4 different equations to incorporate imports into  (beta S I) / N
            rimportmethod = as<int>(pars["importmethod"]);
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
            Shidden = tmp(8);
        }

    private:
        #include <algorithm>
        // rpois that returns an int no larger than second arg
        int myrpois( double rate, int mymax) {
            int result = Rf_rpois( rate );
            return std::min(result, mymax);
        }
    
        // model-dependent variables and methods
        void step(void) {
            // core model definition
            // modify state variables in-place, accumulate for reporting
            double latent_rate;
            int ndeltaR;
            // births adjusted for infant mortality
            // deaths and migrations included in deltaR
            // no deaths from SEI
            ////////////////////////////
            // many in-place, sequential, possible modifications of beta
            /////////////////////////
            if (schooltype == 0) {
                //sin forcing
                beta_now = 
                  rR0*rgamma*
                  (1.0-rbetaforce*
                    cos(2.0*Pi*(day-rschoollag)/365.0));
            };
            if (schooltype == 1) {
                //termtime forcing, schoold scedule passed in as vector of 0/1??
                // add 365 to ensure doy is always positive
                int doy = (day-rschoollag+365) % 365;
                beta_now = rbeta0*pow(1.0+rbetaforce, rschooldays( doy ));
            };
            if (rbetasd_ratio != 0 ) {
                // normal noise perterbation of beta, with sd given as proportion of beta
                beta_now = Rf_rnorm(beta_now, beta_now*rbetasd_ratio);
            }
            //
            // old
            // theta = 0, density dependent, theta=1, freq depend.
            // int nlatent = Rf_rpois( (beta_now*S*Ieff)/floor(pow(N,rtheta)));
            //
            // new -- rtheta = 0, orig model
            // rtheta > 0, small cities have higher beta
            // rtheta < 0, small cities have lower beta
            // rtheta << max city size
            // 
            if (rtheta != 0 ) {
                beta_now = beta_now * ( 1.0+ (rtheta/(abs(rtheta) +N)));
            };
            ////////////////////////////
            // end modifications of beta
            /////////////////////////
            if ( rSfrac != 0 ) {
                // A "hidden" S class, doesn't seem to do much,
                // can only respond if residence time is very small
                // dSh = + a*S - b*Sh;  resid = 1/b, a/b=frac;  a = frac/resid, b = 1/resid
                // do this before anything else
                int S_tohidden = myrpois( (S*rSfrac)/rSresid, S );
                int S_fromhidden = myrpois( Shidden/rSresid, Shidden );
                S += (S_fromhidden - S_tohidden);
                Shidden += (S_tohidden - S_fromhidden);
            }; 
            int nbirth = Rf_rpois( N*rbirth );
            // Effective I is computed at the metapop level for this timestep
            //
            // multiple ways to do latent/imports...??
            if (rimports == 0 ) {
                // scaled for N at metapop level
                latent_rate = (beta_now*S*(Ieff)); 
            } else {
                Rcpp::List debug_report;
                switch ( rimportmethod ) {
                    case 0:
                        // no imports
                        latent_rate = S*(beta_now*I)/N;
                        break;
                    // no connection, all are *S
                    // try 0, 
                    case 1:
                        // constant random imports, no influence of pop 
                        latent_rate = S*((beta_now*I)/N + rimports);
                        if (0) {
                            debug_report["S"] = S;
                            debug_report["beta"] = beta_now;
                            debug_report["I"] = I;
                            debug_report["N"] = N;
                            debug_report["imports"] = rimports;
                            debug_report["latent_rate"] = latent_rate;
                            //if( day % 1000 == 0)  Rf_PrintValue(debug_report);
                        };
                        break;
                    case 2:
                        //  New try
                        //  like 2, only moved N inside
                        //  pulsed inmports, no pop 
                        // divide rimports by rbeta0 so comparable between all
                        latent_rate = S*(beta_now*(I/N + rimports/rbeta0));
                        break;
                    case 3:
                        // constant random imports proportional/modified by to pop
                        //  1.5 is constant param??
                        latent_rate = S*((beta_now*I) + (rimports*pow(N, 1.5)))/N; 
                        break;
                    case 4:
                        // constant random imports proportional/modified by to pop
                        latent_rate = S*(beta_now*(I+ ((rimports/rbeta0)*pow(N, 1.5))))/N; 
                        break;
                    default:
                        throw std::range_error("importmethod must be 0-4");
                        // constant random imports 
                        // doesn't make sense?
                        // latent_rate = S*((beta_now*I) + rimports)/N;
                        /*
                        case 2:
                            // pulsed random imports
                            // divide rimports by rbeta0 so comparable between all
                            // doesn't make sense???
                            latent_rate = S*(beta_now*(I + rimports/rbeta0))/N;
                            break;
                        */
                }
            };
                        
            //int latent_rate = (beta_eff*S*(Ieff+rimports))/N 
            int nlatent = myrpois( latent_rate, S + nbirth);
            int ninfect = myrpois( rsigma*E, E + nlatent );
            int nrecover = myrpois( rgamma*I, I + ninfect );
            // since rpois needs positive rate, use absolute value of rdeltaR
            // to compute total number of events
            ndeltaR = Rf_rpois( N * fabs(rdeltaR)); 
            if ( rdeltaR < 0 ) { 
                // if deltaR is negative, rhen change events to negative
                ndeltaR *= -1;
            }
            S +=  (nbirth - nlatent);
            E +=  (nlatent - ninfect);
            // instantaneous number infected
            I +=  (ninfect - nrecover);
            Eobs += nlatent;
            // prevalence -- all infected in this reporting period (nstep)
            Iobs +=  ninfect;
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
                ret(8) = Shidden;
                // check for negative values
                // Throw an exception if found
                int statemin = min(ret);
                if (statemin < 0) throw_negative_state();
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
            // shouldn't be (1,day)??
            arma::mat tmp(arma::trans(state.cols(1,day)));
            return wrap(tmp);
        }


        NumericMatrix getstate(int whichstate) {
            // get the full history of the given state 
            arma::colvec tmp(arma::trans(state.row(whichstate)));
            return wrap(tmp);
        }
        
        
        NumericVector getday(unsigned int i) {
            // get the full state matrix for this day
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
            if ( pars_.size() < minparlen ) {
                throw std::range_error("pars list has incorrect length");
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

        SEXP throw_negative_state(void) {
            // check that init is the right length
            BEGIN_RCPP
            throw std::range_error("State vector has negative values");
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
        Metapop(SEXP npop_, SEXP xymat_, SEXP cmat_, SEXP blocksize, SEXP nstates_, SEXP nstep) : 
          nstates(as<unsigned int>(nstates_)), npop(as<unsigned int>(npop_)), xymat(xymat_), couplemat(cmat_), thiscity(0)
        {
            SEIR tmpseir(blocksize, nstates_, nstep);
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
            for (unsigned int ii=0; ii<npop; ii++) {
                ret.push_back( pops[ii].report());
            };
            return ret; 
        }
        
        SEXP setcity( unsigned int n_) {
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
            for ( unsigned int ii = 0; ii < npop; ii++) {
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
                   pops[ithis].Ieff += couplemat(ithis, iother) * (pops[iother].I/pops[iother].N);
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
    
    .constructor<SEXP, SEXP, SEXP, SEXP, SEXP, SEXP>("args: npop (int, number of cities), xymat (npop*2 matrix of city locations), couplemat (between-city coupling matrix), blocksize (number of days to preallocate), nstates (number of state variables), nstep (number of days per observation)")
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
        

