#include <Rcpp.h>       
#include <RcppArmadillo.h>
#include "Parlist.h"
#include "SEIR.h"
using namespace Rcpp; 

class Metapop { 
    public:
        Metapop(SEXP npop_, SEXP xymat_, SEXP cmat_, SEXP blocksize,  SEXP obs_nstep) : 
           npop(as<unsigned int>(npop_)), xymat(xymat_), couplemat(cmat_)
        {
            // initialize a single population
            SEIR tmpseir(obs_nstep);
            // then make a vector of them
            pops = std::vector<SEIR>(npop, tmpseir);
            // check that nrow xymat == npop!!??
        }

        void init(SEXP initstates_, SEXP transmat_, SEXP states_, SEXP events_, SEXP accum_, SEXP obsall_, SEXP nobs_) {
            // ?? check that dimensions of initstates matches 
            // states as rows, cities as cols
            NumericMatrix initstates(initstates);
            // create a temp event and then add it to each population
            Events tmpevents = Events(transmat_, states_, events_, accum_, obsall_, nobs_);
            for (unsigned int ii=0; ii<npop; ii++) {
                pops[ii].events = tmpevents;
                // grab the column for this city
                NumericVector tmpinit = initstates( _, ii);
                pops[ii].events.setStates( tmpinit);
                pops[ii].ready();
            };
        };

        NumericMatrix get_metapop_state(int n) { 
            // takes R-style 1-based index of state to retrieve
            // returns matrix of observed history of that state for each pop
            // cities as cols, rows as time
            int nobs = pops[npop].events.iobs;
            arma::mat  ret(nobs, npop);
            for (unsigned int ii=0; ii<npop; ii++) {
                ret.col(ii) = pops[ii].getstate(n-1);
            };
            return wrap (ret);
        } 

        SEXP setpars( SEXP pars) {
            // Takes a list of lists, each sublist contains that city's parlist
            // Loop through and set for each pop
            Rcpp::List tmppars(pars);
            BEGIN_RCPP
            if ( tmppars.size() != npop)  {
                throw std::range_error(" the paramater list of lists should have npop elements");
            };
            for ( unsigned int ii = 0; ii < npop; ii++) {
               pops[ii].events.setpars( tmppars(ii)) ;
            } 
            END_RCPP
        }

        void steps( int n ) {
            // advance the model forward n steps
            for (int ii = 0; ii<n; ii++) {
                prestep();
                step();
                poststep();
            }
        }

    private:
        unsigned int npop;
        std::vector<SEIR> pops;
        NumericMatrix xymat;  //distance matrix between cities
        NumericMatrix couplemat;  //metapop coupling matrix between cities

        void prestep(){
            // functions that need to be completed outside of individual pops
            // pre-intrapop step
            //
            // wrap this into a function??
            // or remove dependence on Ieff -- can we grab statename from Events/SEIR?
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
        };

        void poststep(){
            // stub of anything that needs to be completed outside of individual pops
        };
 
        void step( ) {
            //  intra-population functions
            // take the next step for each pop
            for ( ithis = 0; ithis < npop; ithis++) {
               pops[ithis].steps( 1 );
            } 
        }
};


RCPP_MODULE(seirmod){
	using namespace Rcpp ;
    class_<Metapop>("Metapop")
    
    .constructor<SEXP, SEXP, SEXP, SEXP, SEXP, SEXP>("args: 
                        npop (int, number of cities), 
                        xymat (npop*2 matrix of city locations), 
                        couplemat (between-city coupling matrix), 
                        blocksize (number of observations to preallocate),  
                        obs_step (number of steps per observation)"
    )
    .method("get_metapop_state", &Metapop::get_metapop_state, "args: 
                    int (which state column to return,1-based indexing).  
                    Return matrix of all cities, given state"
    )
    .method("setpars", &Metapop::setpars, "args: 
                    list-of-list of length npop, each containing that city's parameters"
    )
    .method("steps", &Metapop::steps, "args: number of steps to advance all cities")
    .method("init", &Metapop::init, "args:
            initstates_ (nstate*ncity matrix of initial states), 
            transmat    (nstate * nevent matrix, number of transitions per event for each state),
            states_     (named list of states, values don't matter)
            events_     (named list of events, values don't matter), 
            accum_,     (named list as above),
            obsall_     (bool, yes observes states and accum, )
            nobs_       ( int, max number of observations)"
    ) 
    ;
}        
