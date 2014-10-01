#include <RcppArmadillo.h>
#include "Parlist.h"
#include "Pop.h"


//using namespace Rcpp; 

class Metapop { 
    public:
        Metapop(SEXP npop_, SEXP initstates_, SEXP transmat_, SEXP poplist_) : 
           npop(as<std::size_t>(npop_)), istep(0), metapars(), ready_pop(false), ready_metapop(false)
        {
            std::size_t ii;
            // initialize a single population
            Pop tmppop( transmat_, poplist_);
            // then make a vector of them
            pops = std::vector<Pop>(npop, tmppop);
            // should error-checking be here or in R??
            // check: nrow xymat == npop??
            // check: dimensions of initstates matches 
            //  
            // Set initial states of each pop
            // states are rows, pops are cols
            IntegerMatrix initstates(initstates_);
            for (ii=0; ii<npop; ii++) {
                // grab the column for this city
                IntegerVector tmpvec = initstates( _, ii);
                std::vector<unsigned int> tmpinit = as<std::vector<unsigned int> >(tmpvec);
                pops[ii].states.copy( tmpinit);
            };
        };

        IntegerMatrix get_metapop_state(std::size_t n) { 
            // takes R-style 1-based index of state to retrieve
            // returns matrix of observed history of that state for each pop
            // cities as cols, rows as time
            std::size_t ii;
            std::size_t nobs = pops[npop-1].iobs;
            arma::umat  ret(nobs, npop);
            for (ii=0; ii<npop; ii++) {
                ret.col(ii) = pops[ii].getstate(n-1);
            };
            return wrap(ret);
        }   

        // populate metapopulation parameter list, mark as ready
        void set_metapop( SEXP metapars_) {
            metapars.set(metapars_);
            ready_metapop = true;
        }

        void set_pop( SEXP pars_) {
            std::size_t ii;
            // Takes a list of lists, each sublist contains that city's parlist
            // Loop through and set for each pop
            Rcpp::List tmppars(pars_);
            //BEGIN_RCPP
            if ( tmppars.size() != npop)  {
                throw std::range_error(" the paramater list of lists should have npop elements");
            };
            for ( ii = 0; ii < npop; ii++) {
               Rcpp::List thispars = tmppars(ii);
               pops[ii].pars.set( thispars ) ;
            } 
            ready_pop = true;
            //END_RCPP
        }


        void set_school( SEXP tmp_) {
            std::size_t ii;
            // Takes a list of vectors, each vector contains that city's school term
            // Loop through and set for each pop
            Rcpp::List tmp(tmp_);
            //BEGIN_RCPP
            if ( tmp.size() != npop)  {
                throw std::range_error(" the schoo list should have npop elements");
            };
            for ( ii = 0; ii < npop; ii++) {
               IntegerVector thisschool = tmp(ii);
               pops[ii].schoolterm =  thisschool  ;
            } 
            //END_RCPP
        }

        void steps( std::size_t n ) {
            std::size_t ii;
            if (!(ready_pop && ready_metapop)) {
                throw std::runtime_error("Tried to run model before parameter initialization");
            }
            // advance the model forward n steps
            RNGScope scope; // when to call this??
            for (ii = 0; ii<n; ii++) {
                prestep();
                step();
                poststep();
            }
        }

    private:
        std::size_t metaI, metaN ; // Metapop sums of I and N
        // does parlist need to change in-flight? 
        Parlist<double> metapars;
        std::size_t npop;
        std::vector<Pop> pops;
        std::size_t istep; // index of current step, separate from iobs
        bool ready_pop, ready_metapop;  // toogle this to true once pars have been set

        inline void prestep(){
            std::size_t ii;
            // pre-intrapop step
            // includes functions that need to be completed outside of individual pops
            //
            // compute total metapop population and incidence
            // for migration
            //reset
            metaI, metaN = 0;
            for ( ii = 0; ii < npop; ii++) {
                   metaI +=  pops[ii].states("I");
                   metaN +=  pops[ii].states("N");
               // spatially-explicit migration would look something like this: 
               //pops[ithis].Ieff = 0;
               //for ( iother = 0; iother < npop; iother++) {
                    // couplemat isn't symmetric -- step through this city's row
                   //pops[ithis].Ieff += couplemat(ithis, iother) * (pops[iother].I/pops[iother].N);
               //}
            }
        };

        inline void poststep(){
            // stub of anything that needs to be completed outside of individual pops
        };
      
        inline void step( ) {
            std::size_t ithis;
            //  intra-population functions
            // take the next step for each pop
            for ( ithis = 0; ithis < npop; ithis++) {
               pops[ithis].step(istep, metaI, metaN);
            }
            istep++;
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
                deltat (timestep, implement?? ), \
            }" 
    )
    .method("get_metapop_state", &Metapop::get_metapop_state, "args: \
                    int (which state column to return,1-based indexing).  \
                    Return matrix of all cities, given state"
    )
    .method("set_pop", &Metapop::set_pop, "args: \
            list-of-list of length npop, each containing that city's parameters\
            NOTE this must be run before model steps can be taken"
    )
    .method("set_school", &Metapop::set_school, "args: \
            list of school term vectors for each city"
    )
    .method("set_metapop", &Metapop::set_metapop, "args: \
            list of metapop pars   \
            NOTE this must be run before model steps can be taken"
    )
    .method("steps", &Metapop::steps, "args: number of steps to advance all cities")
    ;
}        
