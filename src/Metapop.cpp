#include <Rcpp.h>          
#include <RcppArmadillo.h>
#include "Parlist.h"
#include "Pop.h"


using namespace Rcpp; 

class Metapop { 
    public:
        Metapop(SEXP npop_, SEXP initstates_, SEXP transmat_, SEXP poplist_) : 
           npop(as<unsigned int>(npop_)), istep(0), metapars(), ready(false)
        {
            // initialize a single population
            Pop tmppop( transmat_, poplist_);
            // then make a vector of them
            pops = std::vector<Pop>(npop, tmppop);
            // check that nrow xymat == npop!!??
            //  
            // Set states of each pop
            // ?? check that dimensions of initstates matches 
            // for initstates, states are rows, pops are cols
            NumericMatrix initstates(initstates_);
            for (unsigned int ii=0; ii<npop; ii++) {
                // grab the column for this city
                NumericVector tmpinit = initstates( _, ii);
                pops[ii].states.copy( tmpinit);
            };
        };


        NumericMatrix get_metapop_state(int n) { 
            // takes R-style 1-based index of state to retrieve
            // returns matrix of observed history of that state for each pop
            // cities as cols, rows as time
            int nobs = pops[npop].iobs;
            arma::mat  ret(nobs, npop);
            for (unsigned int ii=0; ii<npop; ii++) {
                ret.col(ii) = pops[ii].getstate(n-1);
            };
            return wrap(ret);
        } 

        SEXP setpars( SEXP metapars_, SEXP pars_) {
            metapars.set(metapars_);
            // Takes a list of lists, each sublist contains that city's parlist
            // Loop through and set for each pop
            Rcpp::List tmppars(pars_);
            BEGIN_RCPP
            if ( tmppars.size() != npop)  {
                throw std::range_error(" the paramater list of lists should have npop elements");
            };
            for ( unsigned int ii = 0; ii < npop; ii++) {
               Rcpp::List thispars = tmppars(ii);
               pops[ii].pars.set( thispars ) ;
            } 
            END_RCPP
            ready = true;
        }

        void steps( int n ) {
            if (!ready) {
                throw std::runtime_error("Tried to run model before parameter initialization");
            }
            // advance the model forward n steps
            RNGScope scope; // when to call this??
            for (int ii = 0; ii<n; ii++) {
                prestep();
                step();
                poststep();
            }
        }

    private:
        // does parlist need to change in-flight? 
        Parlist metapars;
        unsigned int npop;
        std::vector<Pop> pops;
        int istep; // index of current step, separate from iobs
        bool ready;  // toogle this to true once pars have been set

        void prestep(){
            // access matrices via metapars
            // NumericMatrix xymat;  //distance matrix between cities
            // NumericMatrix couplemat;  //metapop coupling matrix between cities
            // FIXME!!
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
                //pops[ithis].Ieff = 0;
               //for ( iother = 0; iother < npop; iother++) {
                    // couplemat isn't symmetric -- step through this city's row
                   //pops[ithis].Ieff += couplemat(ithis, iother) * (pops[iother].I/pops[iother].N);
               //}
            }
        };

        void poststep(){
            // stub of anything that needs to be completed outside of individual pops
        };
 
        void step( ) {
            //  intra-population functions
            // take the next step for each pop
            for ( int ithis = 0; ithis < npop; ithis++) {
               pops[ithis].step(istep);
            } 
        }
};


RCPP_MODULE(seirmod){
	using namespace Rcpp ;
    class_<Metapop>("Metapop")
    .constructor<SEXP, SEXP, SEXP, SEXP>("args: \
            npop (int, number of cities), \
            metapars (list of parameters required by metapop model, might contain for example): \
            {\
                xymat (npop*2 matrix of city locations, currently unused), \
                couplemat (between-city coupling matrix), \
            }, \
            initstates_ (nstate*ncity matrix of initial states), \
            transmat_    (nstate * nevent matrix, number of transitions per event for each state),\
            poplist_ (list containing the following elements):\
            {\
                accum_,     (character vector naming events to accumulate, must be contained in colnames of transmat),\
                obsall_     (bool, yes observes states and accum), \
                nobs_       ( int, max number of observations), \
                obs_nstep (number of steps per observation), \
                deltat (timestep, implement!! ), \
            }" 
    )
    .method("get_metapop_state", &Metapop::get_metapop_state, "args: \
                    int (which state column to return,1-based indexing).  \
                    Return matrix of all cities, given state"
    )
    .method("setpars", &Metapop::setpars, "args: \
                    two lists:  \
            * first, list of metapop pars   \
            * 2nd:  list-of-list of length npop, each containing that city's parameters\
            NOTE this must be run before model steps can be taken"
    )
    .method("steps", &Metapop::steps, "args: number of steps to advance all cities")
    ;
}        
