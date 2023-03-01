/*
NOTES: Currently doing risk grouping crudely for all-event predictions (not in standsurv)
28 June 2021: Adding stgenbrier.ado to calculate brier scores [SM]
*/

/* Calibration plot for competing risks */
program define pmcalplotcr, rclass
        syntax varlist(min=2 max=10 numeric) ///
					[, 	km aj plothide /// select how to generate observed risks
					models(string) ///
					events(numlist integer) /// k_list, must be specified in same order as models
					eventvar(varname numeric) ///
					allpred(varname numeric)  ///
					allnullpred(varname numeric)  /// required for IPA and KM opt
					datavar(string) /// must be coded 0 = derivation 1 = validation; possibly not needed
					timevar(varname numeric)  ///
					Timepoint(int 1) ///
					lp(varlist numeric min=2 max=10) /// 
					Bin(int 10) /// number of risk groups
					RCSDF(int 0) /// df for rcsgen
					RCSPercentiles(numlist ascending) /// percentiles for rcsgen
					RCSKnots(numlist ascending) /// knots for rcsgen
					RCSRmsknots(int 5) /// knots using rmsknots
					binprstub(string) /// to return the risk groups of LP
                    CUTpoints(numlist >=0 <=1 sort) ///
					/*atvar(string)*/ /// same as atvar() in standsurv
					link(string) ///
					spikeplot(string) ///
					riskgroups(varlist numeric min=2 max=10) ///
					cens(int 0) ///
					GENerate(string) ///
					densaxis(real 0.05) ///
					TITLEplot(string) ///
					RANGEplot(string) ///
					histwidth(real 0.001) ///
					EVENTLABel(string)													/// label for events
					Ipcw(varlist)                                                    	/// varlist for IPCW
					GENBrier(string)                                                    /// generate brier score
					GENIpa(string)                                                   	/// generate IPA score
					stpm2df(numlist integer) scale(string) stpm2knots(numlist ascending)		/// stpm2 model options
					ipa(varlist)													    /// calculate ipa option
					brier ///
					grtext ///
					cicurve ///
					]
		
	/* Set up graph style */
	grstyle clear

	set scheme s2color

	grstyle init
	grstyle set plain, horizontal grid noextend
	grstyle set legend 4, nobox
	grstyle set horizontal
	grstyle set margin "0mm 0mm 0mm 0mm": twoway
	grstyle set legend 1, inside
	grstyle color background white // get rid of background shading

	grstyle set color d3 10
	
	/* Error Checks */
	
/*
	if `rcsdf' == 0 {
		local rcsdf 
	}
*/
	
	// Check right options specified for spikeplot()
	if "`spikeplot'" == "" {
		local spikeplot single
	}
	
	local count = wordcount("`spikeplot'")
	if `count' > 1 {
		di in red "Only one of either: dual, single, or none can be specified in spikeplot() option"
		exit 198
	}
	
	if "`spikeplot'" != "dual" & "`spikeplot'" != "single" & "`spikeplot'" != "none" {
		di in red "spikeplot() option must contain one of either: dual, single, or none"
		exit 198
	}
	
	// Check stpsurv/stpci is installed
	capture which stpsurv
	if _rc >0 {
			display in red "You need to install the commands stpci and stpsurv. This can be installed using,"
			display in red ". {net install http://www.stata-journal.com/software/sj15-3/st0202_1}"
			exit  198
	}		
			
	// Check standsurv is installed
	capture which standsurv
	if _rc >0 {
			display in red "You need to install the command standsurv. This can be installed using,"
			display in red ". {net from https://www.pclambert.net/downloads/standsurv }"
			exit  198
	}
	
	// Check stpm2 is installed
	capture which stpm2
	if _rc >0 {
			display in red "You need to install the command stpm2. This can be installed using,"
			display in red ". {ssc install stpm2 }"
			exit  198
	}
	
	if "`models'" == "" {
		di in red "Must specify stored estimates for cause-specific models e.g. models(model1 model2)"
		exit 198
	}
		
	if "`km'" != "" & "`allpred'" == "" {
		di in red "If option km is specified, user must also provide model predictions for the naive 1 - KM for all-causes in option allpred()."
		di in red "This can be done after fitting an all-cause model."
		if "`brier'" != "" & "`ipa'" != "" & "`allnullpred'" == "" {
			di in red "Similarly, predictions for the all-cause null model must be provided in allnullpred() for IPA"
		}
		exit 198
	}
	
	if "`brier'" != "" & "`ipa'" != "" & "`allnullpred'" == "" & "`km'" != "" {
			di in red "If option km is specified and IPA is also required, user must also provide null model predictions for the naive 1 - KM for all-causes in option allnullpred()."
			di in red "This can be done after fitting an all-cause null model."
			exit 198
		}

	if "`km'" == "" & "`aj'" == "" {
		di in red "Use must specify either km or aj as an option to choose whether to produce calibration plots on cause-specific CIFs, or cause-specific NPDs (1-KM)"
		exit 198
	}
	
	if "`km'" != "" & "`aj'" != "" {
		di in red "Either km or aj option must be specified (not both)"
		exit 198
	}
	
	
	// Need to identify deriv/valid data in long format - maybe not needed?
	if "`datavar'" == "" {
		di in red "Need variable that identifies derivation/validation data. Ensure data is in long format where 0 = derivation data and 1 = validation data."
		exit 198
	}
	else {
		local count = wordcount("`datavar'")
		if `count' > 2 {
			di in red "REH REH (sound when something is wrong). Need to specify datavars() with variable that contains indicator for derivation/validation data and numeric value of interest. E.g. datavar(data 1)"
		exit 198
		}
		else {
			tokenize `datavar'
			local datavar_name `1'
			local datavar_num `2'
		}
	}
	

	if "`generate'" == "" {
		local newvar _pred
	}
	else {
		local newvar `generate'
	}
	
	if "`rangeplot'" == "" {
		local rangeplot 0(0.2)1
	}
	

	* how many competing events?
	local K = wordcount("`models'")
		
	* split varlist
	local count = wordcount("`varlist'")
	tokenize `varlist'
	forvalues i = 1/`count' {
		local pred_c`i' ``i''
	}
	
	* split varlist
	local count = wordcount("`lp'")
	tokenize `lp'
	forvalues i = 1/`count' {
		local lp_c`i' ``i''
	}
	
	* split events
	local count = wordcount("`events'")
	
	if "`K'" != "`count'" {
		di in red "Length of number list in events() must be equal to number of models in models()"
		exit 198
	}
	
	tokenize `events'
	forvalues i = 1/`count' {
		local cause`i'_num ``i''
		local cause`i'_name = word("`models'", `cause`i'_num')
		di in green "`cause`i'_name' = `cause`i'_num'"
		* list of competing risks (exc. cause of interest)
		local complist _`i'
		forvalues j = 1/`K' {
			if `i' != `j'{
				local complist_`i' `complist_`i'' ``j'' 
			}
		}
	}


	
	* split risk groups
	local count = wordcount("`riskgroups'")
	tokenize `riskgroups'
	forvalues i = 1/`count' {
		local riskgrp_c`i' ``i''
	}
	qui summ `riskgrp_c1'
	local bin = `r(max)'
	
	if "`csh'" != "" & "`aj'" != "" {
		di in yellow "Caution: Ensure that risk-groups are calculated over predicted RISKS i.e. cause-specific CIF, rather than the linear predictor for CSH models"
	}

 	
	* Generate time to predict risk predictions at
	tempvar _T
	qui gen `_T' = `timepoint' in 1
	
	if "`binprstub'" == "" {
		local binpvar rgrp
	}
	* Generate predictions within risk groups
	forvalues k = 1/`K' {
		qui gen `newvar'_c`k'_`binpvar' = .
		forvalues i=1/`bin' {
			qui summ `pred_c`k'' if `riskgrp_c`k'' == `i' & `datavar_name' == `datavar_num' 
			qui replace `newvar'_c`k'_`binpvar' = r(mean) in `i'
		}
	}
		
	* all-events
	if "`allpred'" != "" {
		local pred_all `allpred'
	}
	if "`allpred'" == "" & "`aj'" != "" {
		qui gen `newvar'_pred_all = 0
		forvalues k = 1/`K' {
			qui replace `newvar'_pred_all = `newvar'_pred_all + `pred_c`k''
		}
		local pred_all `newvar'_pred_all
	}
	
	tempvar cuts_all riskgrp_all 
	pctile `cuts_all' = `pred_all', nquantiles(`bin')
	xtile `riskgrp_all' = `pred_all', cutpoints(`cuts_all')
// 	/qui gen `newvar'_pred_all = `pred_all_setup' + 0

	qui gen `newvar'_all_`binpvar' = .
	forvalues i=1/`bin' {
		qui summ `pred_all' if `riskgrp_all' == `i' & `datavar_name' == `datavar_num' 
		qui replace `newvar'_all_`binpvar' = r(mean) in `i'
	}
	
	// brier score for all-causes
	if "`brier'" != "" {
			
		* split null preds and combine for all-cause
		local count = wordcount("`ipa'")
		tokenize `ipa'
		forvalues i = 1/`count' {
			local nullpred_c`i' ``i''
		}
		
		if "`aj'" != "" {
			if "`allnullpred'" == "" {
				tempvar _nullpred_all
				qui gen `_nullpred_all' = 0
				forvalues k = 1/`K' {
					qui replace `_nullpred_all' = `_nullpred_all' + `nullpred_c`k''
				}
				local ipa_all `_nullpred_all'
			}
			if "`allnullpred'" != "" {
			    local ipa_all `allnullpred'
			}
		}
/*
		if "`km'" != "" {
			local ipa_all `allnullpred'
		}
*/
		if "`aj'" != "" {
			stgenbrier `pred_all', eventlabel(all) btime(`timepoint') ipcw(`ipcw') df(`stpm2df') ipa(`ipa_all') knots(`stpm2knots')

			local brier_null_all : di %7.2f r(b_null_c1)*100
			local brier_model_all  : di %7.2f r(b_model_c1)*100
			local ipa_all  : di %7.2f r(ipa_c1)*100
			
			if "`grtext'" != "" {
				local briernote_all text(0.30 0.4 "{bf:Overall Performance (%):}" 0.21 0.4 "IPA (Brier{subscript:scaled}) = `ipa_all'" 0.12 0.4 "Brier Score (Model) = `brier_model_all'" 0.03 0.4 "Brier Score (Null) = `brier_null_all'", size(vsmall) placement(east) justification(left))
			}
			if "`grtext'" == "" {
				local briernote_all note("{bf:Overall Performance (%):}" " " "IPA (Brier{subscript:scaled}) = `ipa_all'" " " "Brier Score (Model) = `brier_model_all'" " " "Brier Score (Null) = `brier_null_all'", size(vsmall) position(4) ring(0.1) linegap(0.7) justification(right) alignment(middle) fcolor(white))
			}
		}

	}	
		
	***
	if "`km'" != "" { // Calculate predicted net risks (1 - KM)
	
		di in yellow "Caution: Calculating observed/predicted risks on net risk i.e. 1 - KM. Ensure this is what was intended. For observed/predicted risks on cause-specific CIFs, use AJ option."
		
			// Check if we want brier scores
			if "`brier'" != "" {
				
				* split genipa/brier:
				tokenize `genbrier'
				forvalues i = 1/`=wordcount("`genbrier'")' {
					local genbrier_c`i' ``i''
				}
				tokenize `genipa'
				forvalues i = 1/`=wordcount("`genipa'")' {
					local genipa_c`i' ``i''
				}
								
				forvalues k = 1/`K' {
					qui stset `_dta[st_bt]', failure(`eventvar' = `cause`k'_num') scale(`_dta[st_bs]') /// 
												id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
												time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')
					stgenbrier `pred_c`k'', eventlabel(`cause`k'_name') btime(`timepoint') ipcw(`ipcw') genbrier(`genbrier_c`k'')  genipa(`genipa_c`k'') df(`stpm2df') knots(`stpm2knots') scale(`scale') ipa(`nullpred_c`k'')
					di "stgenbrier `pred_c`k'', eventlabel(`cause`k'_name'') btime(`timepoint') ipcw(`ipcw') genbrier(`genbrier_c`k'')  genipa(`genipa_c`k'') df(`stpm2df') scale(`scale') ipa(`nullpred_c`k'')"
					local brier_null_c`k' : di %7.2f r(b_null_c1)*100
					local brier_model_c`k'  : di %7.2f r(b_model_c1)*100
					local ipa_c`k'  : di %7.2f r(ipa_c1)*100
					
					if "`grtext'" != "" {
						local briernote_c`k' text(0.30 0.4 "{bf:Overall Performance (%):}" 0.21 0.4 "IPA (Brier{subscript:scaled}) = `ipa_c`k''" 0.12 0.4 "Brier Score (Model) = `brier_model_c`k''" 0.03 0.4 "Brier Score (Null) = `brier_null_c`k''", size(vsmall) placement(east) justification(left))
					}
					if "`grtext'" == "" {
						local briernote_c`k' note("{bf: Overall Performance (%):}" " " "IPA (Brier{subscript:scaled}) = `ipa_c`k''" " " "Brier Score (Model) = `brier_model_c`k''" " " "Brier Score (Null) = `brier_null_c`k''", size(vsmall) position(4) ring(0.1) linegap(0.7) justification(right) alignment(middle) fcolor(white))
					}
				}
			}	
		
		* Generate predictions for Naive KM	
			
		/* generate pseudovalues within risk groups (mixture KM) */
		
		forvalues k = 1/`K' {
		di in green "Generating pseudovalues for event `cause`k'_name' within `bin' risk groups"
			forvalues i = 1/`bin' {
					qui stset `_dta[st_bt]', failure(`eventvar' = `cause`k'_num') scale(`_dta[st_bs]') /// 
												id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
												time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')
					//di in green "Generating pseudovalues for event `cause`k'_name' within risk group `i'"
					qui stpsurv if `datavar_name' == `datavar_num' & `riskgrp_c`k'' == `i', at(`timepoint') gen(`newvar'_c`k'_pvkmgrp_temp) failure
					qui sum `newvar'_c`k'_pvkmgrp_temp
					qui gen `newvar'_c`k'_pvkmgrp`i' = r(mean)
					cap drop `newvar'_c`k'_pvkmgrp_temp
				}
			* store predictions in long rather than wide - change for PVs
			set varabbrev off
			foreach var in `newvar'_c`k'_pvkmgrp {
				qui gen `var' = .
			}
			forvalues i=1/`bin' {
				foreach var in `newvar'_c`k'_pvkmgrp {
					qui replace `var' = `var'`i'[1] in `i'
					qui drop `var'`i'
				}
			}
			set varabbrev on
		}
		
/*
		* do for all-causes
		di in green "Generating pseudovalues for all events within `bin' risk groups"
		forvalues i = 1/`bin' {
			qui stset `_dta[st_bt]', failure(`eventvar' = "`events'") scale(`_dta[st_bs]') /// 
										id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
										time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')
			//di in green "Generating pseudovalues for all events within risk group `i'"
			qui stpsurv if `datavar_name' == `datavar_num' & `riskgrp_all' == `i', at(`timepoint') gen(`newvar'_all_pvkmgrp_temp) failure
			qui sum `newvar'_all_pvkmgrp_temp
			qui gen `newvar'_all_pvkmgrp`i' = r(mean)
			cap drop `newvar'_all_pvkmgrp_temp
		}
		* store predictions in long rather than wide - change for PVs
		set varabbrev off
		foreach var in `newvar'_all_pvkmgrp {
			qui gen `var' = .
		}
		forvalues i=1/`bin' {
			foreach var in `newvar'_all_pvkmgrp {
				qui replace `var' = `var'`i'[1] in `i'
				qui drop `var'`i'
			}
		}
		set varabbrev on
*/
		
		/* generate PV estimates for 1 - KM using stpsurv: Must be done using mixture KM approach due to dependent censoring */
		forvalues k = 1/`K' {
			
			qui stset `_dta[st_bt]', failure(`eventvar' = `cause`k'_num') scale(`_dta[st_bs]') /// 
										id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
										time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')
										
			qui stpsurv if `datavar_name' == `datavar_num', at(`timepoint') gen(`newvar'_c`k'_pvkm) failure
			
			* Generate spline smoother - better for larger EHRs
			tempvar `pred_c`k''_rcs		
			qui rcsgen `newvar'_c`k'_`binpvar', gen(``pred_c`k''_rcs') df(`rcsdf') percentiles(`rcspercentiles') knots(`rcsknots') rmsknots(`rcsrmsknots')
			qui glm `newvar'_c`k'_pvkmgrp ``pred_c`k''_rcs'?, vce(robust) link(logit) 
			if `e(converged)' == 0 {
				qui glm `newvar'_c`k'_pvkmgrp ``pred_c`k''_rcs'?, vce(robust) link(logit) difficult 
			}
			qui predict `newvar'_c`k'_splnpd
			if "`cicurve'" != "" {
				qui predict `newvar'_c`k'npd_error, stdp
				qui gen `newvar'_c`k'_splnpd_lci = invlogit(logit(`newvar'_c`k'_splnpd) - invnormal(0.975)*`newvar'_c`k'npd_error)
				qui gen `newvar'_c`k'_splnpd_uci = invlogit(logit(`newvar'_c`k'_splnpd) + invnormal(0.975)*`newvar'_c`k'npd_error)
			}
		}
		
/*
		* do for all-causes
		qui stset `_dta[st_bt]', failure(`eventvar' = "`events'") scale(`_dta[st_bs]') /// 
												id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
												time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')
		di in green "Generating pseudovalues for all events"
		qui stpsurv if `datavar_name' == `datavar_num', at(`timepoint') gen(`newvar'_all_pv) failure
		tempvar `pred_all'_rcs													 
		qui rcsgen `newvar'_all_`binpvar', gen(``pred_all'_rcs')  df(`rcsdf') percentiles(`rcspercentiles') knots(`rcsknots') rmsknots(`rcsrmsknots')
		qui glm  `newvar'_all_pvkmgrp ``pred_all'_rcs'?, vce(robust) link(logit)
		if `e(converged)' == 0 {
			qui glm `newvar'_all_pvkmgrp ``pred_all'_rcs'?, vce(robust) link(logit) difficult 
		}
		qui predict `newvar'_all_splnpd
		if "`cicurve'" != "" {
			qui predict `newvar'_allnpd_error, stdp
			qui gen `newvar'_all_splnpd_lci = invlogit(logit(`newvar'_all_splnpd) - invnormal(0.975)*`newvar'_allnpd_error)
			qui gen `newvar'_all_splnpd_uci = invlogit(logit(`newvar'_all_splnpd) + invnormal(0.975)*`newvar'_allnpd_error)
		}
*/
		
		/* prepare spike plot */
		if "`spikeplot'" == "dual" {
			local ifd1 if _d == 1
		}
		
		forvalues k = 1/`K' {
			qui stset `_dta[st_bt]' if `datavar_name' == `datavar_num', failure(`eventvar' = `cause`k'_num') scale(`_dta[st_bs]') /// 
												id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
												time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')
			
			if "`spikeplot'" == "dual" {
				qui twoway__histogram_gen `pred_c`k'' if _d == 0, gen(`newvar'_c`k'_cens_bins `newvar'_c`k'_cens) width(`histwidth')
				local `newvar'_c`k'_cens_width = r(width)
				qui summ `newvar'_c`k'_cens_bins
				qui replace `newvar'_c`k'_cens_bins = ((`newvar'_c`k'_cens_bins/r(max))*-`densaxis') - `densaxis' 
			}   
			   
			qui twoway__histogram_gen `pred_c`k'' `ifd1', gen(`newvar'_c`k'_event_bins `newvar'_c`k'_event) width(`histwidth')
			local `newvar'_c`k'_event_width = r(width)
			qui summ `newvar'_c`k'_event_bins
			qui replace `newvar'_c`k'_event_bins = ((`newvar'_c`k'_event_bins/r(max))*`densaxis') - `densaxis'
			
			if "`spikeplot'" == "single" {
				local spikeplotcode_c`k' (bar `newvar'_c`k'_event_bins `newvar'_c`k'_event, barw(``newvar'_c`k'_event_width') base(-`densaxis') lcolor(white%0) color("214 39 40" ))
			}
			
			if "`spikeplot'" == "dual" {
				local spikeplotcode_c`k' (bar `newvar'_c`k'_event_bins `newvar'_c`k'_event, barw(``newvar'_c`k'_event_width') base(-`densaxis') lcolor(white%0) color("214 39 40" )) ///
			   (bar `newvar'_c`k'_cens_bins `newvar'_c`k'_cens, barw(``newvar'_c`k'_cens_width') base(-`densaxis') lcolor(white%0) color("44 160 44"))
			}
		}
		
/*
		* do for all events
		qui stset `_dta[st_bt]' if `datavar_name' == `datavar_num', failure(`eventvar' = "`events'") scale(`_dta[st_bs]') /// 
												id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
												time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')		
	
		if "`spikeplot'" == "dual" {
			qui twoway__histogram_gen `pred_all' if _d == 0, gen(`newvar'_all_cens_bins `newvar'_all_cens) width(`histwidth')
			local `newvar'_all_cens_width = r(width)
			qui summ `newvar'_all_cens_bins
			qui replace `newvar'_all_cens_bins = ((`newvar'_all_cens_bins/r(max))*-`densaxis') - `densaxis' 
		}
		   
		qui twoway__histogram_gen `pred_all' `ifd1', gen(`newvar'_all_event_bins `newvar'_all_event) width(`histwidth')
		local `newvar'_all_event_width = r(width)
		qui summ `newvar'_all_event_bins
		qui replace `newvar'_all_event_bins = ((`newvar'_all_event_bins/r(max))*`densaxis') - `densaxis'
		
		if "`spikeplot'" == "single" {
			local spikeplotcode_all (bar `newvar'_all_event_bins `newvar'_all_event, barw(``newvar'_all_event_width') base(-`densaxis') lcolor(white%0) color("214 39 40" ))
		}
		
		if "`spikeplot'" == "dual" {
			local spikeplotcode_all (bar `newvar'_all_event_bins `newvar'_all_event, barw(``newvar'_all_event_width') base(-`densaxis') lcolor(white%0) color("214 39 40" )) ///
		   (bar `newvar'_all_cens_bins `newvar'_all_cens, barw(``newvar'_all_cens_width') base(-`densaxis') lcolor(white%0) color("44 160 44"))
		}
*/
		
		if "`plothide'" == "" {
			/* generate cal plots */
			forvalues k = 1/`K' {
				
				if "`cicurve'" != "" {
				    local plotcif_ci (rarea `newvar'_c`k'_splnpd_lci `newvar'_c`k'_splnpd_uci `newvar'_c`k'_`binpvar', fcolor(magenta%10) lpattern(dash) lcolor(magenta%10) sort)
				}
			
				tw (scatter `newvar'_c`k'_pvkmgrp `newvar'_c`k'_`binpvar', m(Oh) color("255 127 14")) ///
				   (line `newvar'_c`k'_splnpd `newvar'_c`k'_`binpvar', sort lcolor(magenta)) ///
				   `spikeplotcode_c`k'' `plotcif_ci' ///
				   (function y=x, range(0 1) lp(-) lcolor(black)), legend(subtitle("{bf:Legend:}", size(vsmall) placement(east)) order(2 "Splines" 3 "Predictor Distribution") textfirst cols(1) size(vsmall) symxsize(4) symysize(4) region(margin(small)) position(2) ring(0.1) placement(east)) ///
				   xlab(`rangeplot', format(%9.2f)) ylab(`rangeplot', format(%9.2f)) ///
				   ytitle("Pseudo-value Net Risk (Mixture)", size(small)) xtitle("Predicted Net Risk (Mixture)", size(small)) ///
				   title("`titleplot'; Event: `cause`k'_name'", size(small)) name(`newvar'_calplot_c`k', replace)  `briernote_c`k''
			
			}
			
/*
			* for all events
			if "`cicurve'" != "" {
				local plotcif_ci (rarea `newvar'_all_splnpd_lci `newvar'_all_splnpd_uci `newvar'_all_`binpvar', fcolor(magenta%10) lpattern(dash) lcolor(magenta%10) sort)
			}
				
			tw (scatter `newvar'_all_pvkmgrp `newvar'_all_`binpvar', m(Oh) color("255 127 14")) ///
			   (line `newvar'_all_splnpd `newvar'_all_`binpvar', sort lcolor(magenta)) ///
			   `spikeplotcode_all' `plotcif_ci' ///
			   (function y=x, range(0 1) lp(-) lcolor(black)), legend(subtitle("{bf:Legend:}", size(vsmall) placement(east)) order(2 "Splines" 3 "Predictor Distribution") textfirst cols(1) size(vsmall) symxsize(4) symysize(4) region(margin(small)) position(2) ring(0.1) placement(east)) ///
			   xlab(`rangeplot', format(%9.2f)) ylab(`rangeplot', format(%9.2f)) ///
			   ytitle("Pseudo-value Net Risk (Mixture)", size(small)) xtitle("Predicted Net Risk (Mixture)", size(small)) ///
			   title("`titleplot'; Event: All", size(small)) name(`newvar'_calplot_all, replace) `briernote_all'
*/
		}
		
	} // end km loop
	***
	
	if "`aj'" != "" { // Calculate predicted cause-specific CIFs
	
			// Check if we want brier scores
			if "`brier'" != "" {
				stgenbrier `varlist', eventlabel(`models') btime(`timepoint') ipcw(`ipcw') genbrier(`genbrier')  genipa(`genipa') df(`stpm2df') knots(`stpm2knots') scale(`scale') ipa(`ipa')		
				forvalues k = 1/`K' {
					local brier_null_c`k' : di %7.2f r(b_null_c`k')*100
					local brier_model_c`k'  : di %7.2f r(b_model_c`k')*100
					local ipa_c`k'  : di %7.2f r(ipa_c`k')*100
					
					if "`grtext'" != "" {
						local briernote_c`k' text(0.30 0.4 "{bf:Overall Performance (%):}" 0.21 0.4 "IPA (Brier{subscript:scaled}) = `ipa_c`k''" 0.12 0.4 "Brier Score (Model) = `brier_model_c`k''" 0.03 0.4 "Brier Score (Null) = `brier_null_c`k''", size(vsmall) placement(east) justification(left))
					}
					if "`grtext'" == "" {
						local briernote_c`k' note("{bf: Overall Performance (%):}" " " "IPA (Brier{subscript:scaled}) = `ipa_c`k''" " " "Brier Score (Model) = `brier_model_c`k''" " " "Brier Score (Null) = `brier_null_c`k''", size(vsmall) position(4) ring(0.1) linegap(0.7) justification(right) alignment(middle) fcolor(white))
					}
				}
			}	

		/* generate PV estimates for CIF using stpci */
		forvalues k = 1/`K' {
		
			qui stset `_dta[st_bt]', failure(`eventvar' = `cause`k'_num') scale(`_dta[st_bs]') /// 
												id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
												time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')
			di in green "Generating pseudovalues for event `cause`k'_name'"
			qui stpci `eventvar' if `datavar_name' == `datavar_num', at(`timepoint') gen(`newvar'_c`k'_pv) ///
																 competingvalues(`complist_`k'')
			tempvar `pred_c`k''_rcs													 
			qui rcsgen `pred_c`k'', gen(``pred_c`k''_rcs')  df(`rcsdf') percentiles(`rcspercentiles') knots(`rcsknots') rmsknots(`rcsrmsknots')
			qui glm  `newvar'_c`k'_pv ``pred_c`k''_rcs'?, vce(robust) link(logit)
			if `e(converged)' == 0 {
				di in yellow "Now trying with difficult option to help convergence"
				qui glm  `newvar'_c`k'_pv ``pred_c`k''_rcs'?, vce(robust) link(logit) difficult 
				if `e(converged)' != 0 {
					di in yellow "Convergence achieved"
				}
			}
			qui predict `newvar'_c`k'_splcif
			if "`cicurve'" != "" {
				qui predict `newvar'_c`k'cif_error, stdp
				qui gen `newvar'_c`k'_splcif_lci = invlogit(logit(`newvar'_c`k'_splcif) - invnormal(0.975)*`newvar'_c`k'cif_error)
				qui gen `newvar'_c`k'_splcif_uci = invlogit(logit(`newvar'_c`k'_splcif) + invnormal(0.975)*`newvar'_c`k'cif_error)
			}

		}
		
		* do for all-causes
		qui stset `_dta[st_bt]', failure(`eventvar' = "`events'") scale(`_dta[st_bs]') /// 
												id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
												time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')
		di in green "Generating pseudovalues for all events"
		qui stpsurv if `datavar_name' == `datavar_num', at(`timepoint') gen(`newvar'_all_pv) failure
		tempvar `pred_all'_rcs													 
		qui rcsgen `pred_all', gen(``pred_all'_rcs')  df(`rcsdf') percentiles(`rcspercentiles') knots(`rcsknots') rmsknots(`rcsrmsknots')
		qui glm  `newvar'_all_pv ``pred_all'_rcs'?, vce(robust) link(logit)
		if `e(converged)' == 0 {
			qui glm  `newvar'_all_pv ``pred_all'_rcs'?, vce(robust) link(logit) difficult 
		}
		qui predict `newvar'_all_splcif
		if "`cicurve'" != "" {
			qui predict `newvar'_allcif_error, stdp
			
			qui gen `newvar'_all_splcif_lci = invlogit(logit(`newvar'_all_splcif) - invnormal(0.975)*`newvar'_allcif_error)
			qui gen `newvar'_all_splcif_uci = invlogit(logit(`newvar'_all_splcif) + invnormal(0.975)*`newvar'_allcif_error)
		}
	
		/* generate PV estimates within groups using mixture approach */
		forvalues k = 1/`K' {
		di in green "Generating pseudovalues for event `cause`k'_name' within `bin' risk groups"
			forvalues i = 1/`bin' {
				qui stset `_dta[st_bt]', failure(`eventvar' = `cause`k'_num') scale(`_dta[st_bs]') /// 
											id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
											time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')
				qui stpci if `datavar_name' == `datavar_num' & `riskgrp_c`k'' == `i', at(`timepoint') gen(`newvar'_c`k'_pvgrp_temp)  competingvalues(`complist_`k'')
				qui sum `newvar'_c`k'_pvgrp_temp
				qui gen `newvar'_c`k'_pvgrp`i' = r(mean)
				cap drop `newvar'_c`k'_pvgrp_temp
			}
			* store predictions in long rather than wide - change for PVs
			set varabbrev off
			foreach var in `newvar'_c`k'_pvgrp {
				qui gen `var' = .
			}
			forvalues i=1/`bin' {
				foreach var in `newvar'_c`k'_pvgrp {
					qui replace `var' = `var'`i'[1] in `i'
					qui drop `var'`i'
				}
			}
			//set varabbrev on
		}
		
		
		
		* do for all-causes
		di in green "Generating pseudovalues for all events within `bin' risk groups"
		forvalues i = 1/`bin' {
			qui stset `_dta[st_bt]', failure(`eventvar' = "`events'") scale(`_dta[st_bs]') /// 
										id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
										time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')
			//di in green "Generating pseudovalues for all events within risk group `i'"
			qui stpsurv if `datavar_name' == `datavar_num' & `riskgrp_all' == `i', at(`timepoint') gen(`newvar'_all_pvgrp_temp) failure
			qui sum `newvar'_all_pvgrp_temp
			qui gen `newvar'_all_pvgrp`i' = r(mean)
			cap drop `newvar'_all_pvgrp_temp
		}
		* store predictions in long rather than wide - change for PVs
		set varabbrev off
		foreach var in `newvar'_all_pvgrp {
			qui gen `var' = .
		}
		forvalues i=1/`bin' {
			foreach var in `newvar'_all_pvgrp {
				qui replace `var' = `var'`i'[1] in `i'
				qui drop `var'`i'
			}
		}
		set varabbrev on
				
		/* prepare spike plot */
		if "`spikeplot'" == "dual" {
			local ifd1 if _d == 1
		}
		
		
		forvalues k = 1/`K' {
			qui stset `_dta[st_bt]' if `datavar_name' == `datavar_num', failure(`eventvar' = `cause`k'_num') scale(`_dta[st_bs]') /// 
												id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
												time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')
			
			if "`spikeplot'" == "dual" {
				qui twoway__histogram_gen `pred_c`k'' if _d == 0, gen(`newvar'_c`k'_cens_bins `newvar'_c`k'_cens) width(`histwidth')
				local `newvar'_c`k'_cens_width = r(width)
				qui summ `newvar'_c`k'_cens_bins
				qui replace `newvar'_c`k'_cens_bins = ((`newvar'_c`k'_cens_bins/r(max))*-`densaxis') - `densaxis' 
			}
			   
			qui twoway__histogram_gen `pred_c`k'' `ifd1', gen(`newvar'_c`k'_event_bins `newvar'_c`k'_event) width(`histwidth')
			local `newvar'_c`k'_event_width = r(width)
			qui summ `newvar'_c`k'_event_bins
			qui replace `newvar'_c`k'_event_bins = ((`newvar'_c`k'_event_bins/r(max))*`densaxis') - `densaxis'
			
			if "`spikeplot'" == "single" {
				local spikeplotcode_c`k' (bar `newvar'_c`k'_event_bins `newvar'_c`k'_event, barw(``newvar'_c`k'_event_width') base(-`densaxis') lcolor(white%0) color("214 39 40" ))
			}
			
			if "`spikeplot'" == "dual" {
				local spikeplotcode_c`k' (bar `newvar'_c`k'_event_bins `newvar'_c`k'_event, barw(``newvar'_c`k'_event_width') base(-`densaxis') lcolor(white%0) color("214 39 40" )) ///
			   (bar `newvar'_c`k'_cens_bins `newvar'_c`k'_cens, barw(``newvar'_c`k'_cens_width') base(-`densaxis') lcolor(white%0) color("44 160 44"))
			}
		}
		
		* do for all events
		qui stset `_dta[st_bt]' if `datavar_name' == `datavar_num', failure(`eventvar' = "`events'") scale(`_dta[st_bs]') /// 
												id(`_dta[st_id]') `_dta[st_show]' exit(`_dta[st_exit]') ///
												time0(`_dta[st_bt0]') enter(`_dta[st_enter]') origin(`_dta[st_orig]')
		
		if "`spikeplot'" == "dual" {
			qui twoway__histogram_gen `pred_all' if _d == 0, gen(`newvar'_all_cens_bins `newvar'_all_cens) width(`histwidth')
			local `newvar'_all_cens_width = r(width)
			qui summ `newvar'_all_cens_bins
			qui replace `newvar'_all_cens_bins = ((`newvar'_all_cens_bins/r(max))*-`densaxis') - `densaxis' 
		}
		   
		qui twoway__histogram_gen `pred_all' `ifd1', gen(`newvar'_all_event_bins `newvar'_all_event) width(`histwidth')
		local `newvar'_all_event_width = r(width)
		qui summ `newvar'_all_event_bins
		qui replace `newvar'_all_event_bins = ((`newvar'_all_event_bins/r(max))*`densaxis') - `densaxis'
		
		if "`spikeplot'" == "single" {
			local spikeplotcode_all (bar `newvar'_all_event_bins `newvar'_all_event, barw(``newvar'_all_event_width') base(-`densaxis') lcolor(white%0) color("214 39 40" ))
		}
		
		if "`spikeplot'" == "dual" {
			local spikeplotcode_all (bar `newvar'_all_event_bins `newvar'_all_event, barw(``newvar'_all_event_width') base(-`densaxis') lcolor(white%0) color("214 39 40" )) ///
		   (bar `newvar'_all_cens_bins `newvar'_all_cens, barw(``newvar'_all_cens_width') base(-`densaxis') lcolor(white%0) color("44 160 44"))
		}
		
		
		/* generate cal plots */
		if "`plothide'" == "" {
			forvalues k = 1/`K' {
				
				if "`cicurve'" != "" {
				    local plotcif_ci (rarea `newvar'_c`k'_splcif_lci `newvar'_c`k'_splcif_uci `pred_c`k'', fcolor(magenta%10) lpattern(dash) lcolor(magenta%10) sort)
				}
				
				
				
				tw (scatter `newvar'_c`k'_pvgrp `newvar'_c`k'_`binpvar', m(Oh) color("255 127 14")) ///
				   (line `newvar'_c`k'_splcif `pred_c`k'', sort lcolor(magenta)) ///
				   `spikeplotcode_c`k'' `plotcif_ci' ///
				   (function y=x, range(0 1) lp(-) lcolor(black)), legend(subtitle("{bf:Legend:}", size(vsmall) placement(east)) order(2 "Splines" 3 "Predictor Distribution") textfirst cols(1) size(vsmall) symxsize(4) symysize(4) region(margin(small)) position(2) ring(0.1) placement(east)) ///
				   xlab(`rangeplot', format(%9.2f)) ylab(`rangeplot', format(%9.2f)) ///
				   ytitle("Pseudo-value CIF", size(small)) xtitle("Predicted CIF", size(small)) ///
				   title("`titleplot'; Event: `cause`k'_name'", size(small)) name(`newvar'_calplot_c`k', replace) `briernote_c`k''
						
			}
			
			* for all events
			
			if "`cicurve'" != "" {
				local plotcif_ci (rarea `newvar'_all_splcif_lci `newvar'_all_splcif_uci `pred_all', fcolor(magenta%10) lpattern(dash) lcolor(magenta%10) sort)
			}
			
			tw (scatter `newvar'_all_pvgrp `newvar'_all_`binpvar', m(Oh) color("255 127 14")) ///
			   (line `newvar'_all_splcif `pred_all', sort lcolor(magenta)) ///
			   `spikeplotcode_all' `plotcif_ci' ///
			   (function y=x, range(0 1) lp(-) lcolor(black)), legend(subtitle("{bf:Legend:}", size(vsmall) placement(east)) order(2 "Splines" 3 "Predictor Distribution") textfirst cols(1) size(vsmall) symxsize(4) symysize(4) region(margin(small)) position(2) ring(0.1) placement(east)) ///
			   xlab(`rangeplot', format(%9.2f)) ylab(`rangeplot', format(%9.2f)) ///
			   ytitle("Pseudo-value CIF", size(small)) xtitle("Predicted CIF", size(small)) ///
			   title("`titleplot'; Event: All", size(small)) name(`newvar'_calplot_all, replace) `briernote_all'
		}
		
	} // end aj loop
	

 	
end
