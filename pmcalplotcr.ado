*! Version 0.0.1 24 July 2020

* Author: Sarwar Mozumder 

* Version 0.0.1, 24 July 2020 - program created [SM].

// pmcalplotcr xb, models(cancer 1 other 2) groups() attime(5) smooth(spline) net

program define pmcalplotcr, rclass sortpreserve
  version 15.1
  syntax varlist(min=1 max=2 numeric) [if] [in],                                ///
    [                                                                           ///
		CRModels(string)                                                        /// List of stored cause-specific models. E.g. crmodels(cancer 1 other 2)
		GRoups(int 10)                                                          /// Number of risk groups                                                    
		ATTime(int 1)                                                           /// Time at which predictions are made for calibration plot
		Smooth(string splines)                                                  /// Smoothing method. E.g. smooth(splines) - default, or smooth(lowess) - computationally intensive & not recommended for larger data
		NET                                                                     /// Calibrate using net probabilities i.e. 1 - S_k(t). This is for checking cause-specific model calibration
		CIF                                                                     /// Calibrate using cause-specific CIFs i.e. F_k(t). This is for checking how well models predict risk
		BINWidth(real 0.001)                                                    /// Width of bins in (histogram) spike plot - this should be small. Wide bins are not recommended.
		CENSored(int 0)                                                         /// Censoring indicator
		CI                                                                      /// Confidence intervals?
		SPIKESCale(real 0.025)                                                  /// Re-scale spike plot
		BYvars(string)                                                          /// Conditional predictions on Z.     
	]		

marksample touse, novarlist
		
/*** Need to do error checking ***/




/*** End error checks ***/


/*** Store options and set tempvars ***/





/*** End option and tempvar checks ***/