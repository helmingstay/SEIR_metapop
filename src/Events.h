using namespace Rcpp;

class Events {
    public:
        Events( SEXP transmat_, SEXP states_, SEXP events_, SEXP accum_):
            // prep states, events, and accum as named lists in R
            // events are cols, states are rows
            transmat(transmat),
            states(states_),
            events(events_),
            //name variables to accumulate
            accum(accum_)
            {
                nstates = states.size();
                nevents = events.size();
                naccum = accum.size();
            }
        void setpars(SEXP pars_) {
            // pass a list of parameters in to change the parlist
            pars.set(pars_);
        }
        Rcpp::List states, events, accum;
        int nstates, nevents, naccum;

    private:
        // transition as columns, state variables as rows
        IntegerMatrix transmat;
        Parlist pars;
        
        // required local variables
        double beta_now, latent_rate, import_rate; 

        #include <algorithm>
        // rpois that returns an int no larger than second arg
        int myrpois( double rate, int mymax) {
            int result = Rf_rpois( rate );
            return std::min(result, mymax);
        }

    public:

        // begin model definition
        void prepEvents(int istep) {
            // function to update the rates
            // takes time
            // Model definition
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
        };                

        void calcEvents(int istep) {
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
            // can we use the names of accum to do this programmatically?
            accum["latent"] += events["latent"]
            accum["imports"] += events["imports"]
            accum["infect"] += events["infect"]
        }
        
        arma::colvec observe(bool all) {
            if (all) {
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
            return ret;
        };
};
