version 17

cd "/home/kmdl463/paper-cal-comprisks/stata-do"

set seed ${seed}

/* Simulate covariates for TRAINING DATA */

forvalues s = 1/2 {
	
	di "Running simulation `s' for derivation data"
	
	drop _all

	 qui  {
	
		* OPTIONS
		set obs ${obs}
		gen id = _n 
		
		* continuous confounder
		// draw continuous vars from MVN
		local corr = 0
		matrix C = J(5,5,`corr') + I(5)*(1-`corr')

		mat li C
		
		drawnorm z1 z2 z3, corr(C) means(0 0 0) sds(1 1 0.3)
		summ z?

		* create a non-linear term 
		replace z1 = z1
		gen zsq1 = z1^2
		* categorical confounder
		tempvar cuts
		gen z6 = runiform()
		pctile `cuts' = z6, nquantiles(4)
		xtile z6_grps = z6, cutpoints(`cuts')
		drop z6
		tab z6_grps, gen(z6_g)

		/* Determine true CSH functions */
		global csh1_lambda1 0.123
		global csh1_gamma1 1.39
		global csh1_lambda2 0.01
		global csh1_gamma2 0.819
		global csh1_pmix 0.417
		
		if `s' == 1 {
			global csh2_lambda1 0.01
			global csh2_gamma1 1.499
			global csh2_lambda2 0.01
			global csh2_gamma2 1.429
			global csh2_pmix 0.507
		}
		
		if ${changecr} == 0 & `s' == 2 {
			global csh2_lambda1 0.01
			global csh2_gamma1 1.499
			global csh2_lambda2 0.01
			global csh2_gamma2 1.429
			global csh2_pmix 0.507
		}
		
		if ${changecr} == 1 & `s' == 2 {
			global csh2_lambda1 0.01
			global csh2_gamma1 1.351
			global csh2_lambda2 0.034
			global csh2_gamma2 1.929
			global csh2_pmix 0.559
		}

		*Covariable effects
		include cov_effects.do

		// generate equations
		gen HR1 = .
		gen HR2 = .
		global cov_c1_mat exp(
		global cov_c2_mat exp(
		cap drop cov
		range cov 1 9 9
		local i = 0
		foreach var in z1 z2 z3 z6_g2 z6_g3 z6_g4 {
			local i = `i' + 1
			replace HR1 = ${HR1_`var'} if cov==`i'
			replace HR2 = ${HR2_`var'} if cov==`i'
			if `i' == 1 {
				global cov_c1_mat ${cov_c1_mat} `var':*(${log_chr1_`var'}) :+ zsq1:*(${log_chr1_zsq`i'}) :+
				global cov_c2_mat ${cov_c2_mat} `var':*(${log_chr2_`var'}) :+ zsq1:*(${log_chr2_zsq`i'}) :+
			}
			else {
				global cov_c1_mat ${cov_c1_mat} `var':*(${log_chr1_`var'}) :+
				global cov_c2_mat ${cov_c2_mat} `var':*(${log_chr2_`var'}) :+
			}
		}
		global cov_c1_mat ${cov_c1_mat} 0 ) 
		global cov_c2_mat ${cov_c2_mat} 0 )

		di  "${cov_c1_mat}"
		di  "${cov_c2_mat}"

		global cov_c1 subinstr("${cov_c1_mat}",":","",.)
		global cov_c2 subinstr("${cov_c2_mat}",":","",.)


		global csh_c1 	((${csh1_pmix}:*${csh1_lambda1}:*${csh1_gamma1}:*{t}:^(${csh1_gamma1}:-1):*exp(-1:*${csh1_lambda1}:*{t}:^${csh1_gamma1}) :+ ///
						(1:-${csh1_pmix}):*${csh1_lambda2}:*${csh1_gamma2}:*{t}:^(${csh1_gamma2}:-1):*exp(-1:*${csh1_lambda2}:*{t}:^${csh1_gamma2})):/ ///
						(${csh1_pmix}:*exp(-1:*${csh1_lambda1}:*{t}:^${csh1_gamma1}) :+ (1:-${csh1_pmix}):*exp(-1:*${csh1_lambda2}:*{t}:^${csh1_gamma2}))) :*${cov_c1_mat} 

		global csh_c2 	((${csh2_pmix}:*${csh2_lambda1}:*${csh2_gamma1}:*{t}:^(${csh2_gamma1}:-1):*exp(-1:*${csh2_lambda1}:*{t}:^${csh2_gamma1}) :+ ///
						(1:-${csh2_pmix}):*${csh2_lambda2}:*${csh2_gamma2}:*{t}:^(${csh2_gamma2}:-1):*exp(-1:*${csh2_lambda2}:*{t}:^${csh2_gamma2})):/ ///
						(${csh2_pmix}:*exp(-1:*${csh2_lambda1}:*{t}:^${csh2_gamma1}) :+ (1:-${csh2_pmix}):*exp(-1:*${csh2_lambda2}:*{t}:^${csh2_gamma2})))  :*${cov_c2_mat}  
				
		/* Simulate competing risks training data */
		survsim survtime dead, hazard(${csh_c1} :+ ${csh_c2}) maxt(10.5)
		replace survtime = 0.003 if survtime < 0.003

		// Perform binomial experiment on p = csh1/csh1+csh2
					
		gen csh1 = ((${csh1_pmix}*${csh1_lambda1}*${csh1_gamma1}*survtime^(${csh1_gamma1}-1)*exp(-${csh1_lambda1}*survtime^${csh1_gamma1}) + ///
						(1-${csh1_pmix})*${csh1_lambda2}*${csh1_gamma2}*survtime^(${csh1_gamma2}-1)*exp(-${csh1_lambda2}*survtime^${csh1_gamma2}))/ ///
						(${csh1_pmix}*exp(-${csh1_lambda1}*survtime^${csh1_gamma1}) + (1-${csh1_pmix})*exp(-${csh1_lambda2}*survtime^${csh1_gamma2}))) * `=${cov_c1}'
						
		gen csh2 = ((${csh2_pmix}*${csh2_lambda1}*${csh2_gamma1}*survtime^(${csh2_gamma1}-1)*exp(-${csh2_lambda1}*survtime^${csh2_gamma1}) + ///
						(1-${csh2_pmix})*${csh2_lambda2}*${csh2_gamma2}*survtime^(${csh2_gamma2}-1)*exp(-${csh2_lambda2}*survtime^${csh2_gamma2}))/ ///
						(${csh2_pmix}*exp(-${csh2_lambda1}*survtime^${csh2_gamma1}) + (1-${csh2_pmix})*exp(-${csh2_lambda2}*survtime^${csh2_gamma2}))) * `=${cov_c2}'
						
		gen p1_binom = cond(dead==1,csh1/(csh1+csh2),0)

		gen event_c1 = cond(rbinomial(1,p1_binom) == 1 & dead == 1, 1, 0)
		gen event_c2 = cond(event_c1 == 0 & dead == 1,1,0)

		/* Simulate censoring times and generate final data with censoring */
		survsim censtime cens, dist(exp) lambda(0.05) maxt(10.5)
		** Censoring distribution **
		stset censtime, f(d = 0)

		gen stime = min(survtime, censtime) 

		gen d = 0 if stime >= 10
		replace d = 1 if (survtime < censtime) & event_c1 == 1
		replace d = 2 if (survtime < censtime) & event_c2 == 1
		replace d = 0 if (censtime < survtime)

		replace stime = 0.003 if stime < 0.003
	}

	save data/${savesim}_`s', replace
	
}
