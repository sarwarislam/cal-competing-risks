/* Set WD */
clear all
cd "/home/kmdl463/paper-cal-comprisks/stata-do/"

/*** Options for do file ***/

global seed 4362642	// set seed 
global obs 50000	// set no. of obs for simulation in each dataset
global simdata_program master_sim_data	// specify do file generating simulated data  
global savesim simdata	// file name of saved simulated dataset

global groups = 15	// choose number of groups for visualising PVs and to break dependency of censoring for calc KM
global cause1 = 1	// indicator for cause 1
global cause2 = 2	// indicator for cause 2
global klistN 1 2	// list of causes
global klist_names c1 c2 all // list of models
global datalist 0 1	// specify for which data to generate cal plots for 0 = derivation, 1 = validation 
global timepoints 10	// list of timepoints to generate calibration plots

global simulation CHANGECR // options: DGM, MISSCOV, MISSFUNC, CHANGECR <---- 

if "$simulation" == "CHANGECR" {
    global changecr 1	// indicator for whether we want simulated data with change in baseline for competing event in {savedata}_2 data
}
if "$simulation" != "CHANGECR" {
    global changecr 0	// indicator for whether we want simulated data with change in baseline for competing event in {savedata}_2 data
}

/*** Generate calibration plots from simulation ***/

set seed ${seed}

* simulate data *
include ${simdata_program}.do
use "/home/kmdl463/paper-cal-comprisks/stata-do/data/${savesim}_1", clear
append using "/home/kmdl463/paper-cal-comprisks/stata-do/data/${savesim}_2", gen(data)

cap drop  __0* // fix for tempvars left in data
cap drop id
gen id = _n

/* generate model predictions */

* fit cause-specific models and store *
cap drop lp_*

stset stime, f(d = 1) id(id)
if "$simulation" == "DGM" {
	stpm2 z1 zsq1 z2 z3 z6_g2 z6_g3 z6_g4 if data == 0, df(4) scale(h) eform // DGM
}
if "$simulation" == "MISSFUNC" {
	stpm2 z1 z2 z3 z6_g2 z6_g3 z6_g4 if data == 0, df(4) scale(h) eform // MISSING FUNC FORM
}
if "$simulation" == "MISSCOV" {
	stpm2 z1 zsq1 z2 z3 if data == 0, df(4) scale(h) eform // MISSING KEY COV
}
if "$simulation" == "CHANGECR" {
	stpm2 z1 zsq1 z2 z3 z6_g2 z6_g3 z6_g4 if data == 0, df(4) scale(h) eform // DGM
}
estimates store c1
predict lp_c1, xbnobaseline

stset stime, f(d = 2) id(id)
stpm2 z1 zsq1 z2 z3 z6_g2 z6_g3 z6_g4 if data == 0, df(4) scale(h) eform
estimates store c2
predict lp_c2, xbnobaseline

stset stime, f(d = 1 2) id(id)
stpm2 z1 zsq1 z2 z3 z6_g2 z6_g3 z6_g4 if data == 0, df(4) scale(h) eform
estimates store all
predict lp_all, xbnobaseline

* predict cause-specific CIFs *
cap drop _T*
cap drop ind_cif*
cap drop cif_rgrp*
foreach time in $timepoints {
	
	qui gen _T`time' = `time' in 1
	* generate individual predictions
	di in green "Predicting cause-specific CIF for time = `time'"
	foreach d in $datalist {
		cap drop _at?*
		qui standsurv, at1(., atif(data == `d')) genind(ind_cif`d'_t`time') timevar(_T`time') cif crmodels(c1 c2) 
		cap drop _at?*
	}

	* generate risk groups for viz purposes 
	di in green "Generating groups over dist. of risk"
	foreach k in $klistN {
		foreach d in $datalist {			
			qui tempvar cuts`d'_`k'_t`time' 
			qui pctile `cuts`d'_`k'_t`time'' = ind_cif`d'_t`time'_c`k', nquantiles(${groups})
			qui xtile cif`time'_rgrp`d'_c`k' = ind_cif`d'_t`time'_c`k', cutpoints(`cuts`d'_`k'_t`time'')
		}
	}
}

* predict cause-specific NPDs *
cap drop ind_npd*
cap drop npd_rgrp*
foreach time in $timepoints {
	foreach k in $klist_names {
		estimates restore `k'
		* generate individual predictions
		di in green "Predicting cause-specific NPD for time = `time'"
		foreach d in $datalist {
			cap drop _at?*
			qui standsurv, at1(., atif(data == `d')) genind(ind_npd`d'_t`time'_`k') timevar(_T`time') failure
			cap drop _at?*
		}
	}
	
	* generate risk groups for viz purposes and to break dependent censoring
	foreach k in $klist_names {
		foreach d in $datalist {
			tempvar cuts`d'_`k'_t`time'
			pctile `cuts`d'_`k'_t`time'' = lp_`k', nquantiles(${groups})
			xtile npd`time'_rgrp`d'_`k' = lp_`k', cutpoints(`cuts`d'_`k'_t`time'')
		}
	}
}


/* generate predictions for null models */ //CHANGE MAXT TO 10.1 TO GET PREDS AT T = 10?

* predict cause-specific CIFs using AJ *
stset stime, f(d = 1) id(id) 
cap drop nullind_cif*
foreach d in $datalist {
    stcompet nullind_cif`d' = ci if data == `d', compet1(2) 
    foreach time in $timepoints {
	    
		gen nullind_cif`d'_c1 = nullind_cif`d' if d == 1 
		qui summ nullind_cif`d'_c1 if _t >= `time'
		gen nullind_cif`d'_t`time'_c1 = r(mean)
		
		gen nullind_cif`d'_c2 = nullind_cif`d' if d == 2
		qui summ nullind_cif`d'_c2 if _t >= `time'
		gen nullind_cif`d'_t`time'_c2 = r(mean)
	}
}

* predict cause-specific NPDs using marginal/mixture KM for dep censoring

foreach k in $klistN {
	foreach time in $timepoints {
		foreach d in $datalist {
		* Fit FPM for competing event to account for dependent censoring
		stset stime, f(d = `=itrim(subinword("$klistN", "`k'","", .))') id(id)
		stpm2 z1 zsq1 z2 z3 z6_g2 z6_g3 z6_g4 if data == `d', df(5) tvc(z6_g2 z6_g3 z6_g4) dftvc(3) scale(h) eform // fit a complicated model for competing event
		
		stset stime, f(d = `k') id(id) // set data for cause of interest
		
		preserve

		cap drop split_time
		stsplit split_time, every(1) // time split every year - this changes the structure of data
		
		predict crcens_prob, s timevar(split_time) // predict probability of not being censored due to competing event at split times
		cap drop weight
		gen weight = 1/crcens_prob // calculate inverse of censored (due to competing event) prob weights 

		stset _t [iweight=weight], enter(_t0) failure(_d==1)		
		sts gen cmkm = f if data == `d'

		sum cmkm if _t == `time'
		local null_npd_c`k' = r(mean)

		restore
		
		gen nullind_npd`d'_t`time'_c`k' = `null_npd_c`k'' if data == `d'
		
		}
	}
}

// for all causes
cap drop allkm*
foreach time in $timepoints {
	foreach d in $datalist {
		stset stime, f(d = $klistN) id(id) exit(time `time')
		sts gen allkm`d' = f if data == `d'
		summ allkm`d' if _t == `time'
		gen nullind_npd`d'_t`time'_all = r(mean) if data == `d'
	}
}
	
************************************************************************************

/* Generate calibration plots */

* Load programs
capture program drop calplotcr
qui include "/home/kmdl463/paper-cal-comprisks/stata-do/calplotcr.ado"
capture program drop stgenbrier
qui include "/home/kmdl463/paper-cal-comprisks/stata-do/stgenbrier.ado"

cap drop indcif_brier*
cap drop indcif_ipa*
cap drop indnpd_brier*
cap drop indnpd_ipa*
cap drop _cif*
cap drop _npd*

foreach time in $timepoints {
	foreach d in $datalist {
		
	if "`d'" == "0" {
		local title "Derivation"
	}
	
	if "`d'" == "1" {
		local title "Validation"
	}
	
	*CIFs
	stset stime if data == `d', f(d = $klistN) id(id)
	calplotcr ind_cif`d'_t`time'_c1 ind_cif`d'_t`time'_c2, aj models(c1 c2) events(1 2) eventvar(d) ///
									 timepoint(`time') bin($groups) spikeplot(single) rcsrmsknots(5) ///
									 lp(lp_c1 lp_c2) riskgroups(cif`time'_rgrp`d'_c1 cif`time'_rgrp`d'_c2) ///
									 datavar(data `d') gen(_cif`d'_t`time') title(`title') histwidth(0.005) densaxis(0.2) ///
									 brier ipcw(z1 zsq1 z2 z3 z6_g2 z6_g3 z6_g4) /*stpm2df(2)*/ stpm2knots(9.9) genbrier(indcif_brier`d'_t`time'_c1 indcif_brier`d'_t`time'_c2) ///
   									 genipa(indcif_ipa`d'_t`time'_c1 indcif_ipa`d'_t`time'_c2) ipa(nullind_cif`d'_t`time'_c1 nullind_cif`d'_t`time'_c2) grtext cicurve magnify(4)
	
	*NPDs
	stset stime if data == `d', f(d = $klistN) id(id)
	calplotcr ind_npd`d'_t`time'_c1 ind_npd`d'_t`time'_c2, km models(c1 c2) events(1 2) eventvar(d) allpred(ind_npd`d'_t`time'_all) ///
									 timepoint(`time') bin(${groups}) spikeplot(single) rcsrmsknots(5) ///
									 lp(lp_c1 lp_c2) riskgroups(npd`time'_rgrp0_c1 npd`time'_rgrp0_c2) ///
									 datavar(data `d') gen(_npd`d'_t`time') title(`title') histwidth(0.005) densaxis(0.2) ///
									 brier ipcw(z1 zsq1 z2 z3 z6_g2 z6_g3 z6_g4) /*stpm2df(2)*/ stpm2knots(9.9) genbrier(indnpd_brier`d'_t`time'_c1 indnpd_brier`d'_t`time'_c2) ///
   									 genipa(indnpd_ipa`d'_t`time'_c1 indnpd_ipa`d'_t`time'_c2) ipa(nullind_npd`d'_t`time'_c1 nullind_npd`d'_t`time'_c2) allnullpred(nullind_npd`d'_t`time'_all) grtext cicurve  magnify(4)		 
									 
	}					 					 
}

	
/* Combine plots for report */
if "$simulation" == "DGM" {
	global saveanalysis dgm_data	// file name of saved analysis dataset
	foreach t in $timepoints {
		grc1leg _cif0_t`t'_calplot_c1 _cif0_t`t'_calplot_c2 _cif0_t`t'_calplot_all _npd0_t`t'_calplot_c1 _npd0_t`t'_calplot_c2 , cols(3) title("Derivation; Time = `t'") name(dtime`t'_cifplots, replace) nocopies position(4) ring(0)
		grc1leg _cif1_t`t'_calplot_c1 _cif1_t`t'_calplot_c2 _cif1_t`t'_calplot_all _npd1_t`t'_calplot_c1 _npd1_t`t'_calplot_c2 , cols(3) title("Validation; Time = `t'") name(vtime`t'_cifplots, replace) nocopies position(4) ring(0)
	}
}
if "$simulation" == "MISSFUNC" {
    global saveanalysis missfunc_data	// file name of saved analysis dataset
	foreach t in $timepoints {
		grc1leg _cif0_t`t'_calplot_c1 _cif0_t`t'_calplot_c2 _cif0_t`t'_calplot_all _npd0_t`t'_calplot_c1 _npd0_t`t'_calplot_c2, cols(3) title("Derivation; Time = `t'") name(time`t'_cifplots, replace) nocopies position(4) ring(0)
	}
}
if "$simulation" == "MISSCOV" {
    global saveanalysis misscov_data	// file name of saved analysis dataset
	foreach t in $timepoints {
		grc1leg _cif0_t`t'_calplot_c1 _cif0_t`t'_calplot_c2 _cif0_t`t'_calplot_all _npd0_t`t'_calplot_c1 _npd0_t`t'_calplot_c2, cols(3) title("Derivation; Time = `t'") name(time`t'_cifplots, replace) nocopies position(4) ring(0)
	}
}
if "$simulation" == "CHANGECR" {
	global saveanalysis changecr_data	// file name of saved analysis dataset
		foreach t in $timepoints {
		grc1leg _cif0_t`t'_calplot_c1 _cif0_t`t'_calplot_c2 _cif0_t`t'_calplot_all _npd0_t`t'_calplot_c1 _npd0_t`t'_calplot_c2, cols(3) title("Derivation; Time = `t'") name(dtime`t'_cifplots, replace) nocopies
		grc1leg _cif1_t`t'_calplot_c1 _cif1_t`t'_calplot_c2 _cif1_t`t'_calplot_all _npd1_t`t'_calplot_c1 _npd1_t`t'_calplot_c2, cols(3) title("Validation; Time = `t'") name(vtime`t'_cifplots, replace) nocopies position(4) ring(0)
	}
}

save "/home/kmdl463/paper-cal-comprisks/stata-do/results/${saveanalysis}", replace



		





