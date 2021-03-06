using namespace Rcpp;

class Pop {
    public:
        Pop( SEXP transmat_, SEXP poplist_):
            // prep states, events, and accum as named lists in R
            // events are cols, states are rows
            transmat(transmat_),
            //named list, naming variables to accumulate
            accum(), 
            pars(),
            // these are set from transmat dimnames
            states(), events(),
            iobs(0)
            {
                // list containing requisite elements for pop, can only pass 6 args to constructor
                // do error-checking in R to make sure we have all of: {obsall, nobs, obs_nstep, accum}
                Rcpp::List poplist(poplist_);
                // should states and observations be reported, or just accum variables
                obsall = as<bool>(poplist["obsall"]);
                nobs = as<std::size_t>(poplist["nobs"]);
                obs_nstep = as<std::size_t>(poplist["obs_nstep"]);
                // !! check that all accum names are in events names
                //
                // pull out row and column names from transmat,
                // get vectors of the names and size
                Rcpp::List transmat_dimnames = transmat.attr("dimnames");
                // init states and events from dimnames lists
                CharacterVector tmprows = transmat_dimnames[0];
                CharacterVector tmpcols = transmat_dimnames[1];
                states.init_fromnames( tmprows );
                events.init_fromnames( tmpcols );
                // accum_ is handed in as a character vector
                // need to coerce to Rcpp::List??
                CharacterVector accum_ = poplist["accum"];
                accum.init_fromnames( accum_ );
                if ( obsall ) {
                    // observe all states plus all accumulators
                    nobsvars = accum.N + states.N;
                } else {
                    // just observe accumulators
                    nobsvars = accum.N;   
                }
                // initialize observation matrix with zero-fill
                obsmat = arma::umat(nobsvars, nobs);
                obsmat.zeros();
            }
        //Rcpp::List states, events, accum;
        std::size_t iobs;  // index of current observation, public so metapop can get final dimensions
        arma::umat obsmat;      // the observation matrix - "final" output 
                                // should be simple to implement grow-as-needed, 
                                // just need a resize method taking int
                                // that updates nobs and resizes obsmat

    private:
        // transition as columns, state variables as rows
        IntegerMatrix transmat;
        bool obsall;
        std::size_t nobsvars;       // number of obsvars == number of rows in obs matrix
        std::size_t nobs;  // total number of observations, == ncol of output obs matrix
        std::size_t obs_nstep;  // timesteps per observations
        // replace by .N and .names, respectively
        //unsigned int nstates, nevents, naccum, nobs;
        //std::vector< std::string>  state_vars, event_vars, accum_vars;
        

    public:
        IntegerVector schoolterm;
        // making the Parlists public exposes their methods to metapop
         // keeping them private may be cleaner, but requires a lot of overhead to set each individually
        Parlist<double> pars; 
        // as unsigned, these will never be negative and therefore can wrap
        Parlist<unsigned int> states, events, accum; // should be unsigned??

        void step(std::size_t istep, std::size_t metaI, std::size_t metaN) {
           calcEvents(istep, metaI, metaN); 
           accumEvents(istep);
           updateStates(istep);
           if ( istep % obs_nstep == 1 ) {
               observe(istep);
           }
        }

        arma::ucolvec getstate(std::size_t whichstate) {
            // get the full history of the given obsvariable 
            // only return weeks that have been recoded
            arma::ucolvec tmp(arma::trans(obsmat.row(whichstate)));
            tmp.resize(iobs);
            return tmp;
        }

        void accumEvents(std::size_t istep) {
            std::size_t ii;
            // use the names of accum list to accumulate events for each timestep
            // accum is zeroed on observation
            for (ii = 0; ii < accum.N; ii++) {
                std::string myname = accum.names[ii];
                // e.g. latents, import, infect
                accum.add( myname, events(myname));
            }
        }

        void updateStates(std::size_t istep) {
            // loop through every event and state, 
            // add events to state according to coefficient
            // 
            // index vars
            std::size_t ievent, istate;
            // coefficient telling state delta resulting from events
            std::size_t coef;
            // names for list access
            std::string event_name, state_name; 
            // debugging 
            //Rf_PrintValue(wrap(istep));
            //Rf_PrintValue(wrap("before"));
            //Rf_PrintValue(wrap(states.get_colvec().t()));
            //Rf_PrintValue(wrap(events.get_colvec().t()));
            //?? iterate over names??
            for ( ievent = 0; ievent < events.N; ievent++) {
                // pull out name as string for accessing lists
                event_name = events.names[ievent];
                // for each state, update if needed
                for ( istate = 0; istate < states.N; istate++) {
                    coef = transmat( istate, ievent);
                    if ( coef == 0 ) continue;  // nothing to do
                    state_name = states.names[istate];
                    // update states based on the cell of the transmission matrix 
                    // times the number of events
                    states.add(state_name, coef * events(event_name));
                }
            }
            //Rf_PrintValue(wrap("after"));
            //Rf_PrintValue(wrap(states.get_colvec().t()));
            //Rf_PrintValue(wrap(events.get_colvec().t()));
        }

        void observe(std::size_t istep) {
            // write accumulated states to persistent matrix
            // and point-in-time states
            // reset accumulated sums 
            arma::ucolvec ret( nobsvars );
            if ( iobs +1 > nobs) {
                // this should never happen, right?
                // see above for implementing dynamic resizing
                throw std::range_error("took too many observations");
            }
            if (obsall) {
                // copy both accum and states
                arma::ucolvec accumvec = accum.get_colvec();
                arma::ucolvec statesvec = states.get_colvec();
                std::copy(accumvec.begin(), accumvec.end(), ret.begin());
                // also copy states
                std::copy(statesvec.begin(), statesvec.end(), ret.begin()+accum.N);
            } else {
                ret = accum.get_colvec();
            }
            // reset the accumulated 
            accum.fill(0);
            // add results to observation matrix, increment observation index
            obsmat.col(iobs) = ret;
            iobs++;
        };

        SEXP throw_negative_state(void) {
            // check that init is the right length
            //BEGIN_RCPP
            throw std::range_error("State vector has negative values");
            NumericVector dummy(1);
            return dummy;
            //END_RCPP
        }

        // these should go in model!!
        #include <algorithm>
        std::size_t myrpois( double rate, std::size_t mymax) {
            // rpois that returns an int no larger than second arg
            std::size_t result = Rf_rpois( rate );
            return std::min(result, mymax);
        }
            

        // begin model definition
        // also see prestep in Metapop
        void calcEvents(std::size_t istep, std::size_t metaI, std::size_t metaN) {
            // this is the *only* user-modifiable function
            // in combination with the transition matrix, this is the model
            // should be a private function that takes pars, states, and istep
            // FIXME!!
            // all other vars are local
            //
            //

            // required local variables
            // these should go in model!!
            double beta_now, latent_rate, import_rate; 
            double Pi = arma::datum::pi;

            ////////////////////////////
            // many in-place, sequential, possible modifications of beta
            /////////////////////////
            if (pars("schooltype") == 0) {
                beta_now = pars("R0") * pars("gamma"); 
                //sin forcing
                beta_now = 
                  pars("R0") * pars("gamma") * 
                  ( 1.0-pars("betaforce")* cos(2.0*Pi*istep/365.0));
            };
            //// checkme!!
            if (pars("schooltype")== 1) {
                //termtime forcing, school scedule passed in as vector of 1/-1
                // add 365 to ensure doy is always positive
                std::size_t doy = (istep+1) % 365;
                beta_now = pars("beta0")*pow(1.0+pars("betaforce"), schoolterm( doy ));
            };
            //
            //
            // Effective I is computed at the metapop level for this timestep
            //
            // multiple ways to do latent/imports...??
            // everyone gets the same internal dynamics
            latent_rate = (beta_now*states("S")*states("I"))/states("N"); 
            // manual imports, multiply by S at the end
            switch ( static_cast<std::size_t>(pars("importmethod"))) {
                //metapop
                case -1:
                    import_rate = (pars("imports")*(beta_now*(metaI-states("I"))))/(metaN * pow(states("N"), pars("pop.pow")));
                    break;
                // no imports
                case 0:
                    import_rate = 0;
                    break;
                // no connection, all are *S
                // try 0, 
                case 1:
                    // constant random imports, no influence of pop 
                    import_rate = pars("imports");
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
                    import_rate = beta_now*(pars("imports")/pars("beta0"));
                    break;
                // implement
                //case 3:
                    // constant random imports proportional/modified by to pop
                    // import_rate  = (pars("imports")*pow(states["N"], pars("import_power")))/states["N"];
                 //   break;
                //case 4:
                    // pulsed random imports proportional/modified by to pop
                    //import_rate = (beta_now*(states["I"]+ 
                    //              ((pars("imports")/pars("beta0"))
                    //                  * pow(states["N"], pars("import_power")))))/states["N"]; 
                    //break;
                default:
                    throw std::range_error("importmethod not implemented");
                        break;
            } // end switch

                // multiply imports by S at the end
                import_rate *= states("S");


            // done with prep, now compute events 
            // why N??
            //events.list["birth"] = myrpois( states("N") * pars("birth"), states("N")) ;
            events["dS"] = Rf_rpois( states("N") * pars("dS"));
            // imports aren't limited / no mass balance
            if (import_rate != 0 ) {
                events["imports"] = Rf_rpois( import_rate );
            };
            // Should current events be included in ceiling?
            //events.list["latent"] = myrpois( latent_rate, states("S") + events("birth"));
            //events.list["infect"] = myrpois( pars("sigma") * states("E"), states("E") + events("latent"));
            //events.list["recover"] = myrpois( pars("gamma") * states("I"), states("I") + events("infect"));
            //  
            //  This formulation is order-independent
            events["latent"] = myrpois( latent_rate, states("S"));
            events["infect"] = myrpois( pars("sigma") * states("E"), states("E"));
            events["recover"] = myrpois( pars("gamma") * states("I"), states("I"));
            if ( pars("dR")< 0 ) { 
                // if deltaR is negative, only remove R max 
                events["dR"] = myrpois( states("N")* fabs(pars("dR")), states("R")); 
                // then set negative
                events["dR"] = events("dR") * -1;
            } else {
                // otherwise we're adding
                events["dR"] = Rf_rpois( states("N")* fabs(pars("dR"))); 
            };
        };


};
