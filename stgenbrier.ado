*! 1.00 Sarwar Mozumder 29Sep2023

program define stgenbrier, rclass 

version 15 

syntax [varlist(default=none)] [if] [in] ,   ///
		BTime(real)	                                                            /// desired timepoint for the estimation 
		[ 																		///
			EVENTLABel(string)													/// label for events
			Ipcw(varlist)                                                    	/// varlist for IPCW
			GENBrier(string)                                                    /// generate brier score
			GENIpa(string)                                                   	/// generate IPA score
			df(numlist integer) knots(numlist ascending) scale(string) 							/// stpm2 model options
			ipa(varlist)													    /// calculate ipa option
		]

	marksample touse
/*
	cap drop touse
	qui gen touse = 1 if `touse'
*/
	
	local event "`_dta[st_bd]'"
	local event_vals "`_dta[st_ev]'"
	local event_vals = strtrim(subinstr("`event_vals'", ",", " ", .))

	* how many competing events?
	
	local K = wordcount("`varlist'")
	if `K' < 2 {
		di in yellow "Calculating Brier score for non-competing risks since K = 1. For competing risks, include more than 1 prediction in varlist."	
	} 
	
	
	* split varlist
	local count = wordcount("`varlist'")
	tokenize `varlist'
	forvalues i = 1/`count' {
		local Fi_c`i' ``i''
	}
	
	if `count' == 1 {
		local num_all
		local I = wordcount("`event_vals'")
		tokenize `event_vals'
		forvalues i = 1/`I' {
			local num_all `num_all' ,``i''
		}
		//local num_all substr(`num_all',1,length(`num_all')-1)
	}
	

	* split events
	local count = wordcount("`event_vals'")

	if `count' != wordcount("`varlist'") & wordcount("`varlist'") > 1  {
		di in red "varlist must be equal to the number of events specified in stset"
		exit 198
	}
	
	if `count' > 1 {
		di in green "Events are being assigned as:"
		tokenize `event_vals'
		forvalues i = 1/`count' {
			local cause`i'_num ``i''
			local cause`i'_name = word("`eventlabel'", `cause`i'_num')
			di in green "`cause`i'_name' = `cause`i'_num'"
			* list of competing risks (exc. cause of interest)
			local complist _`i'
			forvalues j = 1/`K' {
				if `i' != `j'{
					local complist_`i' `complist_`i'' ``j'' 
				}
			}
		}
	}
	
	if wordcount("`varlist'") == 1 {
			local cause1_num `event_vals'
			local cause1_name = word("`eventlabel'", 1)
			di in green "Events are being assigned as:"
			di in green "`cause1_name' = `cause1_num'"
	}
	
	
	// Generate Gt = the adjusted censoring probability at btime //
	qui stset _t if `touse', failure(`event'=0) //exit(time `btime') // check this whether to include btime
	if "`scale'"=="" {
		local scale "hazard"
	}
	qui stpm2 `ipcw' if `touse', df(`df') knots(`knots') scale(`scale') eform
	tempvar t1 Gt1 Gti
	qui gen `t1' = `btime'
	qui predict `Gt1', surv timevar(`t1')
	// Generate Gt_i = adjusted censoring probability at time _t //
	qui predict `Gti', surv timevar(_t)

	// reset failure in stset
	qui stset _t if `touse', failure(`event' = `event_vals') exit(`_dta[st_exit]')

/*
	* Number who die before time t=btime
	count if _t<=`btime' & _d==1
	* Number who die or are censored after time t=btime
	count if _t>`btime' 
	* Number who are censored before time t=btime
	count if _t<=`btime' & _d==0 
*/
	 
	if "`genbrier'" == "" {
		forvalues k = 1/`K' {
			tempvar genbrier_c`k'
			local bs_c`k' `genbrier_c`k''
		}
	}
	if "`genbrier'" != "" {
		forvalues k = 1/`K' {
			local bs_c`k' = word("`genbrier'", `k')
			confirm new variable `bs_c`k''
		}
	} 

	if `K' > 1 {
	    forvalues k = 1/`K' {
			qui gen `bs_c`k'' = cond(_t <= `btime' & `event' == `cause`k'_num' & _d == 1 & `touse', (1/`Gti') * ((1 - `Fi_c`k'')^2), 0) ///
								+ cond(_t <= `btime' & `event' == `complist_`k'' & _d == 1 & `touse', (1/`Gti') * ((0 - `Fi_c`k'')^2), 0) ///
								+ cond(_t > `btime' & `touse', (1/`Gt1') * ((0 - `Fi_c`k'')^2), 0) 
		}
	}
	if `K' == 1 {
	
	    forvalues k = 1/`K' {
			qui gen `bs_c`k'' = cond(_t <= `btime' & inlist(`event' `num_all') & _d == 1 & `touse', (1/`Gti') * ((1 - `Fi_c`k'')^2), 0) ///
								+ cond(_t > `btime' & `touse', (1/`Gt1') * ((0 - `Fi_c`k'')^2), 0) 
		}
	}
	
	
	
	if "`ipa'" != "" {
		
    
		local count = wordcount("`ipa'")
		tokenize `ipa'
		forvalues i = 1/`count' {
			local Fi_c`i'_null ``i''
		}


		
		// reset failure in stset
		qui stset _t if `touse', failure(`event' = `event_vals') exit(`_dta[st_exit]')
		
		if "`genipa'" == "" {
			forvalues k = 1/`K' {
				tempvar genipa_c`k'
				local ipa_c`k' `genipa_c`k''
			}
		}
		if "`genipa'" != "" {
			forvalues k = 1/`K' {
				local ipa_c`k' = word("`genipa'", `k')
				confirm new variable `ipa_c`k''
			}
		} 

		forvalues k = 1/`K' {
			tempvar null_bs_c`k'
		}
		
		if `K' > 1 {
			forvalues k = 1/`K' {
				qui gen `null_bs_c`k'' = cond(_t <= `btime' & `event' == `cause`k'_num' & _d == 1 & `touse', (1/`Gti') * ((1 - `Fi_c`k'_null')^2), 0) ///
									+ cond(_t <= `btime' & `event' == `complist_`k'' & _d == 1 & `touse', (1/`Gti') * ((0 - `Fi_c`k'_null')^2), 0) ///
									+ cond(_t > `btime' & `touse', (1/`Gt1') * ((0 - `Fi_c`k'_null')^2), 0) 
			}
		}
		if `K' == 1 {
			forvalues k = 1/`K' {
				qui gen `null_bs_c`k'' = cond(_t <= `btime' & inlist(`event' `num_all') & _d == 1 & `touse', (1/`Gti') * ((1 - `Fi_c`k'_null')^2), 0) ///
									+ cond(_t > `btime' & `touse', (1/`Gt1') * ((0 - `Fi_c`k'_null')^2), 0) 
			}
		}

		forvalues k = 1/`K' {
			
			qui summ `null_bs_c`k'' if `touse'
			local null_brier`k' = r(mean)
			qui summ `bs_c`k'' if `touse'
			local brier`k' = r(mean)
			
			qui gen `ipa_c`k'' = 1 - (`brier`k''/`null_brier`k'') if `touse'
			
		}
		
	}
	
	forvalues k = 1/`K' {
		// get mean brier score and save as scalar
		quietly mean `bs_c`k'' if `touse'
		
		// display the results
		tempname meanout n mean se ll ul level
		local `level' = r(level)
		local `n' = e(N)
		matrix `meanout' = r(table)
		matrix colnames `meanout' = "brier"
		scalar `mean' = `meanout'[1,1]
		scalar `se'   = `meanout'[2,1]
		scalar `ll'   = `meanout'[5,1]
		scalar `ul'   = `meanout'[6,1]
		
		
		// get mean brier score (null) and save as scalar
		quietly mean `null_bs_c`k'' if `touse'
		
		// display the results
		tempname meanout_null n_null mean_null se_null ll_null ul_null level_null
		local `level_null' = r(level)
		local `n_null' = e(N)
		matrix `meanout_null' = r(table)
		matrix colnames `meanout_null' = "brier (null)"
		scalar `mean_null' = `meanout_null'[1,1]
		scalar `se_null'   = `meanout_null'[2,1]
		scalar `ll_null'   = `meanout_null'[5,1]
		scalar `ul_null'   = `meanout_null'[6,1]
		
		// get mean IPA and save as scalar
		qui mean `ipa_c`k'' if `touse'
		
		// display the results
		tempname meanout_ipa n_ipa mean_ipa se_ipa ll_ipa ul_ipa level_ipa
		local `level_ipa' = r(level)
		local `n_ipa' = e(N)
		matrix `meanout_ipa' = r(table)
		matrix colnames `meanout_ipa' = "ipa"
		scalar `mean_ipa' = `meanout_ipa'[1,1]
		scalar `se_ipa'   = .
		scalar `ll_ipa'   = .
		scalar `ul_ipa'   = .

		disp ""
		//disp "Mean estimation for event cause`k'_name                 Number of obs   = "  %9.0g as result  ``n_ipa''
		disp ""
		tempname mytab
		.`mytab' = ._tab.new, col(5) lmargin(0)
		.`mytab'.width    22   |12    12    12    12
		.`mytab'.titlefmt  .     .    . %24s   .
		.`mytab'.pad       .     2     3     3    3
		.`mytab'.numfmt    . %9.0g %9.0g %9.0g %9.0g
		.`mytab'.strcolor result  .  .  . . 
		.`mytab'.strfmt    %21s  .  .  . .
		.`mytab'.strcolor   text  .  .  . .
		.`mytab'.sep, top
		.`mytab'.titles "Event `cause`k'_name'"                                     /// 1
						"Mean"                                                          /// 2
						"   Std. Err."                          /// 3
						"[``level''% Conf. Interval]" ""        //  4 5
		.`mytab'.sep, middle										
		.`mytab'.row    "Brier score (Null)"                           ///
				`mean_null'                                                                  ///
				`se_null'                                                                    ///
				`ll_null'                                                                    ///
				`ul_null'															
		.`mytab'.row    "Brier score (Model)"                           ///
				`mean'                                                                  ///
				`se'                                                                    ///
				`ll'                                                                    ///
				`ul'															
		.`mytab'.row    "IPA"                           ///
				`mean_ipa'                                                                  ///
				`se_ipa'                                                                    ///
				`ll_ipa'                                                                    ///
				`ul_ipa'						
		.`mytab'.sep, bottom
	
	tempname output_c`k'
	matrix `output_c`k'' = (`meanout_null',`meanout',`meanout_ipa')
		
	return scalar b_null_c`k' = `output_c`k''["b","brier(null)"]
	return scalar b_model_c`k' = `output_c`k''["b","brier"]
	return scalar ipa_c`k' = `output_c`k''["b","ipa"]
	
	return matrix table = `output_c`k''
	
	}

end
