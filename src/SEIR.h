using namespace Rcpp;

class SEIR {
    public:
        SEIR(  SEXP obs_nstep_):
          //spos(spos_), 
          obs_nstep(as<int>(obs_nstep_)), 
          istep(0), iobs(0), isready(false)
          {
          }

    private:    // private variable defs
        unsigned int istep;         // current step of the model
        int obs_nstep;              // number of steps per observation
        bool isready;           // flip this switch when Events get filled, ok to run model

    private:
        void ready() {
            isready = true;
        }

        void step(void) {
            if( !isready ) {
                throw_not_initialize();
            };
            calcEvents(istep);
            accumEvents(istep);
            updateStates(istep);
            if ( istep % obs_nstep == 1 ) {
                observe(istep);
            }
        }

        
    public:
        Events events;
        
        SEXP steps(const unsigned int nsteps) {
            // advance the model this many steps
            RNGScope scope; // when to call this?
            BEGIN_RCPP  //check to make sure we don't fall off the end
            for ( unsigned int i =0; i<nsteps; i++) {
                step();             // move forwards one step
                istep++;
            }
            END_RCPP
        }
    
        arma::colvec getstate(int whichstate) {
            // get the full history of the given obsvariable 
            // only return weeks that have been recoded
            arma::colvec tmp(arma::trans(events.obsmat.row(whichstate)));
            tmp.resize(iobs);
            return tmp;
        }
        
};

