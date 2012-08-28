using namespace Rcpp;

class SEIR {
    public:
        SEIR( SEXP blocksize_, SEXP nstates_, SEXP obs_nstep_):
          //spos(spos_), 
          blocksize(as<int>(blocksize_)), 
          nstates(as<unsigned int>(nstates_)), istep(0), iobs(0), 
          obs_nstep(as<int>(obs_nstep_)), 
          Pi(4*atan(1)) 
          {
            state = arma::mat(nstates, blocksize );
            // initialize to 0 to avoid surprises...
            state.zeros();
            maxsteps = blocksize * obs_nstep;
          }

    private:    // private variable defs
        Parlist pars;
        arma::mat state;            // the state matrix - "final" output - grow as needed
        unsigned int blocksize;     // grow state by this much when needed
        unsigned int maxsteps;      // current ncol of state - must grow if more are needed
        unsigned int nstates;       // number of states == number of rows in state matrix
        unsigned int istep;         // current step of the model
        int obs_nstep;              // number of steps per observation
        double Pi;  //pi

    public:
        int iobs;               // index (and number) of observations so far, 
                                // accessed by metapop to determine output size
        // model-dependent variables and methods
        // these need to be easily accessed by metapop
        // use doubles avoid confusion when multiplying by params
        // fundamental state variables
        double S, E, I, R, N;
        // observations or derived state variables
        double Estep, Ieff, Istep, Eimport;
        // derived parameters
        double beta_now;

        void setpars(SEXP pars_) {
            // pass a list of parameters in to change the parlist
            pars.set(pars_);
        }
        
        void setstate(SEXP init_) {
            // this is essentially part of the model definition
            // init state of current step
            NumericVector tmp(init_);
            checkstate(tmp);
            S = tmp(0);
            E = tmp(1);
            I = tmp(2);
            R = tmp(3);
            N = S+E+I+R;
            Ieff = 0;
            // observaton variables, initialize to zero
            Estep = 0;
            Istep = 0;
            Eimport = 0;
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
            // separate internal from external (imports)
            double latent_rate, import_rate;
            int ndeltaR, nbirth, nimport, nlatent, ninfect, nrecover;
            // births adjusted for infant mortality
            // deaths and migrations included in deltaR
            // no deaths from SEI
            ////////////////////////////
            // many in-place, sequential, possible modifications of beta
            /////////////////////////
            if (pars.i("schooltype") == 0) {
                beta_now = pars.d("R0") * pars.d("gamma"); 
                //sin forcing
                beta_now = 
                  pars.d("R0") * pars.d("gamma") * 
                  ( 1.0-pars.d("betaforce")* cos(2.0*Pi*(istep-pars.i("schoollag"))/365.0));
            };
            //// checkme!!
            //if (pars.i("schooltype")== 1) {
                //termtime forcing, schoold scedule passed in as vector of 0/1??
                // add 365 to ensure doy is always positive
                //int doy = (istep-pars.i("schoollag")+365) % 365;
                //beta_now = pars.d("beta0")*pow(1.0+pars.d("betaforce"), pars.NV("schooldays")( doy ));
            //};
            /*if (rbetasd_ratio != 0 ) {
                // normal noise perterbation of beta, with sd given as proportion of beta
                beta_now = Rf_rnorm(beta_now, beta_now*rbetasd_ratio);
            }
            */
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
            //if (rtheta != 0 ) {
                //beta_now = beta_now * ( 1.0+ (rtheta/(abs(rtheta) +N)));
            //};
            ////////////////////////////
            // end modifications of beta
            /////////////////////////
            /*if ( rSfrac != 0 ) {
                // A "hidden" S class, doesn't seem to do much,
                // can only respond if residence time is very small
                // dSh = + a*S - b*Sh;  resid = 1/b, a/b=frac;  a = frac/resid, b = 1/resid
                // do this before anything else
                int S_tohidden = myrpois( (S*rSfrac)/rSresid, S );
                int S_fromhidden = myrpois( Shidden/rSresid, Shidden );
                S += (S_fromhidden - S_tohidden);
                Shidden += (S_tohidden - S_fromhidden);
            }; 
            */
            nbirth = Rf_rpois( N * pars.d("birth") );
            // Effective I is computed at the metapop level for this timestep
            //
            // multiple ways to do latent/imports...??
            // everyone gets the same internal dynamics
            latent_rate = (beta_now*S*I)/N; 
            if (pars.d("imports")== 0 ) {
                // Metapop!
                // Ieff is scaled for N at metapop level
                // changed so self-connect == 0
                import_rate = beta_now*S*Ieff; 
            } else {
                // manual imports, multiply by S at the end
                switch ( pars.i("importmethod")) {
                    case 0:
                        // no imports
                        import_rate = 0;
                        break;
                    // no connection, all are *S
                    // try 0, 
                    case 1:
                        // constant random imports, no influence of pop 
                        import_rate = pars.d("imports");
                        if (0) {
                            // example of how to print conditional debugging information
                            // make a list, add desired elements, and then print
                            Rcpp::List debug_report;
                            debug_report["S"] = S;
                            debug_report["beta"] = beta_now;
                            debug_report["I"] = I;
                            debug_report["N"] = N;
                            debug_report["latent_rate"] = latent_rate;
                            //if( istep % 1000 == 0)  Rf_PrintValue(debug_report);
                        };
                        break;
                    case 2:
                        //  like 2, only moved N inside
                        //  pulsed inmports, no pop 
                        // divide rimports by rbeta0 so comparable between all
                        import_rate = beta_now*(pars.d("imports")/pars.d("beta0"));
                        break;
                    case 3:
                        // constant random imports proportional/modified by to pop
                        //  1.5 is constant param??
                        import_rate  = (pars.i("imports")*pow(N, 1.5))/N;
                        break;
                    case 4:
                        // pulsed random imports proportional/modified by to pop
                        import_rate = (beta_now*(I+ ((pars.i("imports")/pars.i("beta0"))*pow(N, 1.5))))/N; 
                        break;
                    default:
                        throw std::range_error("importmethod must be 0-4");
                            break;
                } // end switch
                // multiply case 0-4 imports by S at the end
                import_rate *= S;
            };
                        
            // imports are not constrained by max # events, e.g. not mass-balance
            nimport = Rf_rpois( import_rate);
            // number of new exposed rrom local dynamics
            nlatent = myrpois( latent_rate, S + nbirth) ;
            // add in imports
            nlatent += nimport;
            ninfect = myrpois( pars.d("sigma")*E, E + nlatent );
            nrecover = myrpois( pars.d("gamma")*I, I + ninfect );
            // since rpois needs positive rate, use absolute value of rdeltaR
            // to compute total number of events
            ndeltaR = Rf_rpois( N * fabs(pars.d("deltaR"))); 
            if ( pars.d("deltaR")< 0 ) { 
                // if deltaR is negative, then make sure there are enough R to remove
                if (ndeltaR > R ) {
                    // try to stay as close to mass balance as possible
                    // without negative states
                    // e.g. in extinction, we run out of R, so move as many 
                    // from S to R as needed
                    int diff = (ndeltaR -R);
                     R += diff;
                     if(S > diff) {
                        S -= diff;
                    }
                }
                //  then change events to negative
                ndeltaR *= -1;
            }
            S +=  (nbirth - nlatent);
            E +=  (nlatent - ninfect);
            // instantaneous number infected
            I +=  (ninfect - nrecover);
            // prevalence -- all infected in this reporting period (iobs)
            Estep += nlatent;
            Eimport += nimport;
            Istep +=  ninfect;
            R += (nrecover + ndeltaR);
            N += (nbirth + ndeltaR);
            if ( istep % obs_nstep == 1 ) {
                iobs++;
                arma::colvec ret(nstates);
                ret.zeros();
                ret(0) = S;
                ret(1) = E;
                ret(2) = Estep;
                ret(3) = I;
                ret(4) = Istep;
                ret(5) = Rf_rbinom(Istep, pars.d("pobs"));
                ret(6) = R;
                ret(7) = N;
                ret(8) = Eimport;
                //ret(4) = Rf_rbinom(I, rpobs);
                state.col(iobs) = ret;
                // set step sums to zero
                Estep = 0;
                Eimport = 0;
                Istep = 0;
                // check for negative values
                // Throw an exception if found
                int statemin = min(ret);
                if (statemin < 0) {
                    Rf_PrintValue(wrap(istep));
                    Rf_PrintValue(wrap(ret));
                    throw_negative_state();
                }
            }
        };
        
    public:
        // model-independent public methods
        void rewind(int nstep){
            // rewind the model this many steps
            istep -= nstep;
        }

        SEXP steps(const unsigned int nsteps) {
            // advance the model this many steps
            RNGScope scope; // when to call this?
            BEGIN_RCPP  //check to make sure we don't fall off the end
            for ( unsigned int i =0; i<nsteps; i++) {
                if ( istep +1 > maxsteps) {
                    //impelemnt resize of state arma matrix here
                    throw std::range_error("took too many steps -- try a larger blocksize");
                }
                step();             // move forwards one step
                istep++;
            }
            END_RCPP
        }
    
        arma::colvec getstate(int whichstate) {
            // get the full history of the given state 
            // only return weeks that have been recoded
            arma::colvec tmp(arma::trans(state.row(whichstate)));
            tmp.resize(iobs);
            return tmp;
        }
        
        
        NumericVector getstep(unsigned int i) {
            // get the full state matrix for this step
            arma::rowvec tmp(state.col(i));
            return wrap( tmp );
        }
        
        void fillstate( int fillval ) {
            // fill the state with int (typically 0?)
            state.fill(fillval);
        } 

    private:
        SEXP checkstate(NumericVector init_) {
            // check that init is the right length
            BEGIN_RCPP
            if ( init_.size() != 4 )  {
                throw std::range_error("init vector has incorrect dimensions -- SEIR");
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

