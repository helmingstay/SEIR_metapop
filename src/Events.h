using namespace Rcpp;

class Events {
    public:
        Events( SEXP transmat_, SEXP accum_, SEXP obsall_, SEXP nobs_):
            // prep states, events, and accum as named lists in R
            // events are cols, states are rows
            transmat(transmat),
            //name variables to accumulate
            accum(accum_),
            obsall(as<bool>(obsall_)), 
            nobs(as<int>(nobs_)),
            iobs(0),
            Pi(arma::datum::pi)
            {
                // pull out row and column names from transmat,
                // get vectors of the names and size
                Rcpp::List transmat_dimnames = transmat_attr.attr("dimnames");
                // states
                state_vars = mk_names( transmat_dimnames[0] );
                states = mk_list(state_vars);
                nstates = state_vars.size();
                // events
                event_vars = mk_names( transmat_dimnames[1] );
                events = mk_list(event_vars);
                nevents = event_vars.size();
                // accum is handed in as a named list
                 accum_vars = mk_names( accum.attr("names") );
                naccum = accum.size();
                if ( obsall ) {
                    // observe all states plus all accumulators
                    nobsvars = naccum + nstates;
                } else {
                    // just observe accumulators
                    nobsvars = naccum;   
                }
                // initialize observation matrix with zero-fill
                obsmat = arma::mat(nobsvars, nobs);
                obsmat.zeros();
            }
        void setpars(SEXP pars_) {
            // pass a list of parameters in to change the parlist
            pars.set(pars_);
        }
        Rcpp::List states, events, accum, iobs;
        arma::mat obsmat;      // the observation matrix - "final" output - grow as needed

    private:
        // transition as columns, state variables as rows
        IntegerMatrix transmat;
        bool obsall;
        std::vector< std::string>  state_vars, event_vars, accum_vars;
        unsigned int blocksize;     // grow observation by this much when needed
        unsigned int maxsteps;      // current ncol of observation - must grow if more are needed
        unsigned int nobsvars;       // number of obsvars == number of rows in obs matrix
        unsigned int nstates, nevents, naccum, nobs;
        Parlist pars;
        
        // required local variables
        double beta_now, latent_rate, import_rate; 
        double Pi;

        #include <algorithm>
        int myrpois( double rate, int mymax) {
            // rpois that returns an int no larger than second arg
            int result = Rf_rpois( rate );
            return std::min(result, mymax);
        }
        std::vector< std::string> mk_names( Rcpp::List names_list ) {
            // function to turn list of names (from attributes) into vector suitable for list indexing
            std::vector< std::string> ret = Rcpp::as< std::vector< std::string> >(names_list);
            return ret;
        }
        Rcpp::List mk_list( std::vector< std::string> list_names ) {
            // function to turn named list from vector of names
            // initialize list by size
            Rcpp::List ret( list_names.size() );
            // set names
            ret.attr("names") = list_names;
            return ret;
        }

    public:
        void setStates(NumericVector newstates) {
            // check range
            if ( newstates.size() != states.size() ){
                stop  // throw exception, FIXME!!
            }
            std::copy(newstates.begin(), newstates.end(), states.begin());
        }

        // begin model definition
        };                

        void calcEvents(int istep) {
            // this is the *only* user-modifiable function
            // in combination with the transition matrix, this is the model
            // should be a private function that takes pars, states, and istep
            // FIXME!!
            // all other vars are local
            
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
            //
            //
            // Effective I is computed at the metapop level for this timestep
            //
            // multiple ways to do latent/imports...??
            // everyone gets the same internal dynamics
            latent_rate = (beta_now*states["S"]*states["I"])/states["N"]; 
            if (pars.d("imports")== 0 ) {
                // Metapop!
                // Ieff is scaled for N at metapop level
                // changed so self-connect == 0
                import_rate = beta_now*states["S"]*states["Ieff"]; 
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
                        //if (0) {
                            // example of how to print conditional debugging information
                            // make a list, add desired elements, and then print
                            //Rcpp::List debug_report;
                            //debug_report["S"] = S;
                            //debug_report["beta"] = beta_now;
                            //debug_report["I"] = I;
                            //debug_report["N"] = N;
                            //debug_report["latent_rate"] = latent_rate;
                            //if( istep % 1000 == 0)  Rf_PrintValue(debug_report);
                        //};
                        break;
                    case 2:
                        //  like 2, only moved N inside
                        //  pulsed inmports, no pop 
                        // divide rimports by rbeta0 so comparable between all
                        import_rate = beta_now*(pars.d("imports")/pars.d("beta0"));
                        break;
                    // implement
                    //case 3:
                        // constant random imports proportional/modified by to pop
                        // import_rate  = (pars.i("imports")*pow(states["N"], pars.d("import_power")))/states["N"];
                     //   break;
                    //case 4:
                        // pulsed random imports proportional/modified by to pop
                        //import_rate = (beta_now*(states["I"]+ 
                        //              ((pars.i("imports")/pars.i("beta0"))
                        //                  * pow(states["N"], pars.d("import_power")))))/states["N"]; 
                        //break;
                    default:
                        throw std::range_error("importmethod not implemented");
                            break;
                } // end switch
                // multiply imports by S at the end
                import_rate *= states["S"];
            };


            // done with prep, now compute events 
            events["birth"] = myrpois( states["N"] * pars.d("birth"), states["N"]) ;
            // imports aren't limited / no mass balance
            if (import_rate != 0 ) {
                events["imports"] = Rf_rpois( import_rate );
            };
            events["latent"] = myrpois( latent_rate, states["S"]+ events["birth"]) ;
            events["infect"] = myrpois( pars.d("sigma") * states["E"], states["E"] + events["latent"]);
            events["recover"] = myrpois( pars.d("gamma") * states["I"], states["I"] + events["infect"]);
            if ( pars.d("deltaR")< 0 ) { 
                // if deltaR is negative, only remove R max 
                events["deltaR"] = myrpois( states["N"]* fabs(pars.d("deltaR")), states["R"]); 
                // then set negative
                events["deltaR"] *= -1;
            } else {
                // otherwise we're adding
                events["deltaR"] = Rf_rpois( states["N"]* fabs(pars.d("deltaR"))); 
            };
        };

        void accumEvents(int istep) {
            // this doesn't need to be user-modifiable
            // can we use the names of accum to do this programmatically?
            for (int ii = 0; ii < naccum; ii++) {
                std::vector myname = accum_vars[ii];
                // e.g. latents, import, infect
                accum[myname] += events[myname];
            }
        }

        void updateStates(int istep) {
            int ievent, istate;
            std::string event_name, state_name; 
            for ( ievent = 0; ievent < nevents; ievent++) {
                event_name = event_names[ievent];
                for ( istate = 0; istate < nstates; istates++) {
                    state_name = state_names[istate];
                    // update states based on the cell of the transmission matrix 
                    // times the number of events
                    states[state_name] +=  ( transmat( istate, ievent) * events[event_name]);
                }
            }
            // then check for negative values
            // Throw an exception if found
            int statemin = min(states);
            if (statemin < 0) {
                Rf_PrintValue(wrap(istep));
                Rf_PrintValue(wrap(states));
                throw_negative_state();
            }
        }

        void observe() {
            if ( iobs +1 > nobs) {
                //TODO?? impelemnt resize of state arma matrix here
                // Does nobs argument supercede this?
                throw std::range_error("took too many observations");
            }
            if (obsall) {
                arma::colvec ret(nevents + naccum);
                // get accumulated and actual, concatenate
                for (int i = 0; i< naccum; i++) {
                    ret[i] = accum[i];
                }
                for (int i = 0; i< nevents; i++) {
                    ret[naccum+i] = nevents[i];
                }
            } else {
                // otherwise only return accumulated
                arma::colvec ret = as<arma::colvec>(naccum);
            }
            // reset the accumulated 
            std::fill(accum.begin(), accum.end(), 0);
            obsmat.col(iobs) = ret;
            iobs++;
        };

        SEXP throw_negative_state(void) {
            // check that init is the right length
            BEGIN_RCPP
            throw std::range_error("State vector has negative values");
            NumericVector dummy(1);
            return dummy;
            END_RCPP
        }
};
