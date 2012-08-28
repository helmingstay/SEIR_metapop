#include <Rcpp.h>       
#include <RcppArmadillo.h>
#include "transition.h"
#include "SEIR.h"
using namespace Rcpp; 

class Metapop { 
    public:
        Metapop(SEXP npop_, SEXP xymat_, SEXP cmat_, SEXP blocksize, SEXP nstates_, SEXP obs_nstep) : 
          nstates(as<unsigned int>(nstates_)), npop(as<unsigned int>(npop_)), xymat(xymat_), couplemat(cmat_), thiscity(0)
        {
            SEIR tmpseir(blocksize, nstates_, obs_nstep);
            pops = std::vector<SEIR>(npop, tmpseir);
            //pops.resize(npop);
            //for (int ii = 0; ii<npop; ii++) {
            //  pops[ii].fillstate(0);
            //}
            // check that nrow xymat == npop!!
        }

        NumericMatrix get_metapop_state(int n) { 
            // takes R-style 1-based index of state to retrieve
            // returns matrix of observed history of that state for each pop
            int nobs = pops[1].iobs;
            arma::mat  ret(nobs, npop);
            for (unsigned int ii=0; ii<npop; ii++) {
                ret.col(ii) = pops[ii].getstate(n-1);
            };
            return wrap (ret);
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
            // should initial values of states be set at metapop initialization?
            NumericMatrix tmpmat(citystates);
            BEGIN_RCPP
            if ( (tmpmat.rows() != npop)  || (tmpmat.ncol() != 4) ) {
                throw std::range_error("improper dimensions of state matrix -- npop rows and 4 cols (SEIR)");
            };
            for ( int ii = 0; ii < npop; ii++) {
               pops[ii].setstate( wrap(tmpmat(ii,_))) ;
            } 
            END_RCPP
        }

        SEXP setpars( SEXP pars) {
            Rcpp::List tmppars(pars);
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
    
    .constructor<SEXP, SEXP, SEXP, SEXP, SEXP, SEXP>("args: npop (int, number of cities), xymat (npop*2 matrix of city locations), couplemat (between-city coupling matrix), blocksize (number of observations to preallocate), nstates (number of state variables), obs_step (number of steps per observation)")
    /*
    .method("steps", &SEIR::steps, "args: nday (int).  Advance model by nday steps" )
    .method("get", &SEIR::get, "args: none. Return full model state-matrix so far.")
    .method("getday", &SEIR::getday, "args: day (int).  Return model state vector for specified day." )
    .method("setpars", &SEIR::setpars, "args: pars (named list).  Post-initialization change of model parameters.")
    .method("setstate", &SEIR::setstate, "args: statevec (integer vector of state values). Change model state for current day.")
    .method("rewind", &SEIR::rewind, "args: nday (int).  Rewind model by this many days.")
    */
    .method("get_metapop_state", &Metapop::get_metapop_state, "args: int, which state column to return (1-based indexing).  Return matrix of all cities, given state")
    .method("setcity", &Metapop::setcity, "args: int citynumber")
    .method("setstates", &Metapop::setstates, "args: npop*4 (SEIR) matrix of states for current timestep")
    .method("setpars", &Metapop::setpars, "args: list-of-list of length npop, each containing that city's parameters")
    .method("steps", &Metapop::steps, "args: number of steps to advance all cities")
    ;
}        
