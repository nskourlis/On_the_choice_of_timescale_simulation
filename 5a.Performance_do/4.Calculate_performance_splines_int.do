
/* Do file summarizing the results af the analysis step*/
/* Deriving the bias, coverage, convergence, relative efficiency, relative bias*/
/* Monte carlo errors, empirical standard errors and more*/

// Approach d- Splines/Int

// Summarise results
clear

cd  "$N\5b.Performance_results"


	
	/*****************************************************************/
	/*****Approach d- Splines/Int****************************************/
	/****************************************************************/

	foreach scen in   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24   {
			
		forval female=0(1)1 {
		global female=`female'
					  
	// load stored results
	use  "$N\4b.Analyze_results\Approach_d_splines_int/scenario`scen'_spline_int.dta", replace

	include "$N/1b.Scenarios_folder/Scenario`scen'_new.do"


	
	sort givenage

	//Estimate bias for the CIF from each cause, for each time point we are interested over a series of ages in diagnosis
	
	foreach Z in CIF1_f${female}_at1 CIF1_f${female}_at2 CIF1_f${female}_at3 CIF1_f${female}_at4 CIF1_f${female}_at5 CIF1_f${female}_at6 CIF1_f${female}_at7 CIF1_f${female}_at8 CIF1_f${female}_at9 CIF1_f${female}_at10   ///
                 CIF2_f${female}_at1 CIF2_f${female}_at2 CIF2_f${female}_at3 CIF2_f${female}_at4 CIF2_f${female}_at5 CIF2_f${female}_at6 CIF2_f${female}_at7 CIF2_f${female}_at8 CIF2_f${female}_at9 CIF2_f${female}_at10 {
	
	
	forvalues i = 60(10)90 {
	
	// load truths
	include  "$N/2b.Truths_folder/scenario`scen'_Truth_`i'.do"
	
		//summarize estimations
	gen `Z'_`i'=`Z' if givenage==`i'
	summ `Z'_`i'
    local `Z'_`i'_t1t1= `r(mean)'
	
	//Summarize bias
	gen `Z'_bias_`i'=.
	    
	replace `Z'_bias_`i'= `Z' - ${`Z'_true} if givenage==`i'
	
	summ `Z'_bias_`i'	
	}
	
	}
	
  
	
	//Calculate the coverages
	
	foreach X in CIF1_f${female}_at1 CIF1_f${female}_at2 CIF1_f${female}_at3 CIF1_f${female}_at4 CIF1_f${female}_at5 CIF1_f${female}_at6 CIF1_f${female}_at7 CIF1_f${female}_at8 CIF1_f${female}_at9 CIF1_f${female}_at10   ///
                 CIF2_f${female}_at1 CIF2_f${female}_at2 CIF2_f${female}_at3 CIF2_f${female}_at4 CIF2_f${female}_at5 CIF2_f${female}_at6 CIF2_f${female}_at7 CIF2_f${female}_at8 CIF2_f${female}_at9 CIF2_f${female}_at10 {
	 
	 forvalues i = 60(10)90 {
	 
	 // load truths
	 include  "$N/2b.Truths_folder/scenario`scen'_Truth_`i'.do"
	
	 
	 //Generate an index variable of whether in the specific repetition, the 95% CI
	 //of the estimate includes the true value of the parameter
	 gen cov_index_`X'_`i'= .
	 replace cov_index_`X'_`i'=0 if givenage==`i' & canc_model_converged==1 & other_model_converged==1
	 replace cov_index_`X'_`i'=1 if `X'_LL<=${`X'_true} & `X'_UL>=${`X'_true} & givenage==`i'
	 
	 summ cov_index_`X'_`i' if cov_index_`X'_`i'!= .
	 gen coverage_`X'_`i'= `r(mean)'
	 local coverage_`X'_`i'_t1t1= coverage_`X'_`i'
	}
	}
	
	
	// Derive the bias, convergence, relative bias and empirical standard errors 
	
	
	forvalues i = 60(10)90 {
	
		foreach var in  CIF1_f${female}_at1_bias_`i' CIF1_f${female}_at2_bias_`i' CIF1_f${female}_at3_bias_`i' CIF1_f${female}_at4_bias_`i' CIF1_f${female}_at5_bias_`i' ///
	                CIF1_f${female}_at6_bias_`i' CIF1_f${female}_at7_bias_`i' CIF1_f${female}_at8_bias_`i' CIF1_f${female}_at9_bias_`i' CIF1_f${female}_at10_bias_`i'   ///
                    CIF2_f${female}_at1_bias_`i' CIF2_f${female}_at2_bias_`i' CIF2_f${female}_at3_bias_`i' CIF2_f${female}_at4_bias_`i' CIF2_f${female}_at5_bias_`i' ///
					CIF2_f${female}_at6_bias_`i' CIF2_f${female}_at7_bias_`i' CIF2_f${female}_at8_bias_`i' CIF2_f${female}_at9_bias_`i' CIF2_f${female}_at10_bias_`i' {
	summ `var' 
	local b_`var'_t1t1= `r(mean)'
	local r_b_`var'_t1t1=`b_`var'_t1t1'
	local convergence_t1t1= `r(N)'

	}
	}
	
			//MCerrors and EmpSE
	
	
	foreach Z in CIF1_f${female}_at1 CIF1_f${female}_at2 CIF1_f${female}_at3 CIF1_f${female}_at4 CIF1_f${female}_at5 ///
				 CIF1_f${female}_at6 CIF1_f${female}_at7 CIF1_f${female}_at8 CIF1_f${female}_at9 CIF1_f${female}_at10   ///
                 CIF2_f${female}_at1 CIF2_f${female}_at2 CIF2_f${female}_at3 CIF2_f${female}_at4 CIF2_f${female}_at5  ///
				 CIF2_f${female}_at6 CIF2_f${female}_at7 CIF2_f${female}_at8 CIF2_f${female}_at9 CIF2_f${female}_at10 {
	
		
	forvalues i = 60(10)90 {
	
	include  "$N/2b.Truths_folder/scenario`scen'_Truth_`i'.do"
	
	gen sqdiff_`Z'_`i'=. 
	
	summ `Z' if givenage==`i'
	
	local mean_`Z'_`i'=`r(mean)'
	
	replace sqdiff_`Z'_`i'= (`Z'-`mean_`Z'_`i'')^2  if  givenage==`i'
	
	sum sqdiff_`Z'_`i' 
	
    local sum_sqdiff_`Z'_`i'= `r(mean)'*`r(N)'
	
    local MCer_`Z'_`i'_t1t1= sqrt((1/(`r(N)'*(`r(N)'-1)))*`sum_sqdiff_`Z'_`i'')
	
	local rb_MCer_`Z'_`i'_t1t1= (`MCer_`Z'_`i'_t1t1'/ ${`Z'_true} )*100
	
	local EmpSE_`Z'_`i'_t1t1= sqrt((1/(`r(N)'-1))*`sum_sqdiff_`Z'_`i'')
	
	//di `MCer_`Z'_`i'_t1t1'
	//di `EmpSE_`Z'_`i'_t1t1'
	
	}
	
	}
	
	
	
/********************************/
	//Gender bias 
/********************************/

	foreach Z in hr_gender_at1 hr_gender_at2 hr_gender_at3 hr_gender_at4 hr_gender_at5    ///
                 hr_gender_at6 hr_gender_at7 hr_gender_at8 hr_gender_at9 hr_gender_at10 {
	
	forvalues i = 60(10)90 {
    
	// load truths
	include  "$N/2b.Truths_folder/scenario`scen'_Truth_`i'.do"

	gen `Z'_bias_`i'=.
	    
	replace `Z'_bias_`i'= `Z' - ${`Z'_true} if givenage==`i'
	summ `Z'_bias_`i'	
	}
	
	}
	
		forvalues i = 60(10)90 {
	
		foreach var in  hr_gender_at1_bias_`i' hr_gender_at2_bias_`i' hr_gender_at3_bias_`i' hr_gender_at4_bias_`i' hr_gender_at5_bias_`i' ///
	                    hr_gender_at6_bias_`i' hr_gender_at7_bias_`i' hr_gender_at8_bias_`i' hr_gender_at9_bias_`i' hr_gender_at10_bias_`i'    {
	summ `var' 
	local b_`var'_t1t1= `r(mean)'

	}
	}
	
	
/********************************/
//Gender  true values
/********************************/
	forvalues i = 60(10)90 {
	foreach V in hr_gender_at1 hr_gender_at2 hr_gender_at3 hr_gender_at4 hr_gender_at5 /// 
			 hr_gender_at6 hr_gender_at7 hr_gender_at8 hr_gender_at9 hr_gender_at10 { 

	gen `V'_`i'= .
	replace `V'_`i'= `V' if givenage==`i'
				}
				 }
				
	forvalues i = 60(10)90 {
	 
	foreach var in hr_gender_at1_`i' hr_gender_at2_`i' hr_gender_at3_`i' hr_gender_at4_`i' hr_gender_at5_`i' /// 
			 hr_gender_at6_`i' hr_gender_at7_`i' hr_gender_at8_`i' hr_gender_at9_`i' hr_gender_at10_`i' { 
		summ `var' 
	local `var'_t1t1= `r(mean)'
	}
	}
	
	
	//CIF diffs, CIF ratios, total CIF ratio
		forvalues i = 60(10)90 {
				include  "$N/2b.Truths_folder/scenario`scen'_Truth_`i'.do"

	foreach V in    CIF1_diff_at1 CIF1_diff_at2 CIF1_diff_at3 CIF1_diff_at4 CIF1_diff_at5 /// 
				    CIF1_diff_at6 CIF1_diff_at7 CIF1_diff_at8 CIF1_diff_at9 CIF1_diff_at10 ///
					CIF2_diff_at1 CIF2_diff_at2 CIF2_diff_at3 CIF2_diff_at4 CIF2_diff_at5 ///
					CIF2_diff_at6 CIF2_diff_at7 CIF2_diff_at8 CIF2_diff_at9 CIF2_diff_at10 ///
					CIF1_ratio_at1 CIF1_ratio_at2 CIF1_ratio_at3 CIF1_ratio_at4 CIF1_ratio_at5 ///
					CIF1_ratio_at6 CIF1_ratio_at7 CIF1_ratio_at8 CIF1_ratio_at9 CIF1_ratio_at10 ///
					CIF2_ratio_at1 CIF2_ratio_at2 CIF2_ratio_at3 CIF2_ratio_at4 CIF2_ratio_at5 ///
					CIF2_ratio_at6 CIF2_ratio_at7 CIF2_ratio_at8 CIF2_ratio_at9 CIF2_ratio_at10 ///
					CIF_total_at1 CIF_total_at2 CIF_total_at3 CIF_total_at4 CIF_total_at5 ///
					CIF_total_at6 CIF_total_at7 CIF_total_at8 CIF_total_at9 CIF_total_at10  ///
					  {
 
			 
	gen `V'_`i'= .

	replace `V'_`i'= `V' if givenage==`i'
				}
				 }
				
	forvalues i = 60(10)90 {
				include  "$N/2b.Truths_folder/scenario`scen'_Truth_`i'.do"

	 
	foreach var in 	CIF1_diff_at1_`i' CIF1_diff_at2_`i' CIF1_diff_at3_`i' CIF1_diff_at4_`i' CIF1_diff_at5_`i' /// 
				    CIF1_diff_at6_`i' CIF1_diff_at7_`i' CIF1_diff_at8_`i' CIF1_diff_at9_`i' CIF1_diff_at10_`i' ///
					CIF2_diff_at1_`i' CIF2_diff_at2_`i' CIF2_diff_at3_`i' CIF2_diff_at4_`i' CIF2_diff_at5_`i' ///
					CIF2_diff_at6_`i' CIF2_diff_at7_`i' CIF2_diff_at8_`i' CIF2_diff_at9_`i' CIF2_diff_at10_`i' ///
					CIF1_ratio_at1_`i' CIF1_ratio_at2_`i' CIF1_ratio_at3_`i' CIF1_ratio_at4_`i' CIF1_ratio_at5_`i' CIF1_ratio_at6_`i' ///
				    CIF1_ratio_at7_`i' CIF1_ratio_at8_`i' CIF1_ratio_at9_`i' CIF1_ratio_at10_`i' ///
				    CIF2_ratio_at1_`i' CIF2_ratio_at2_`i' CIF2_ratio_at3_`i' CIF2_ratio_at4_`i' CIF2_ratio_at5_`i' ///
				    CIF2_ratio_at6_`i' CIF2_ratio_at7_`i' CIF2_ratio_at8_`i' CIF2_ratio_at9_`i' CIF2_ratio_at10_`i' ///
				    CIF_total_at1_`i' CIF_total_at2_`i' CIF_total_at3_`i' CIF_total_at4_`i' CIF_total_at5_`i' ///
				    CIF_total_at6_`i' CIF_total_at7_`i' CIF_total_at8_`i' CIF_total_at9_`i' CIF_total_at10_`i'  ///
					  {
 
	summ `var' 
	local `var'_t1t1= `r(mean)'
 
	}
	}
	

	/*We use the Approach a- Attained age performance measures as well in order to later derive
	the relative bias and relative effeiciency measures that require both the curret approach and 
	the approach of reference*/
	
	/*****************************************************************/
	/*******Approach a- Attained age***********************************/
	/****************************************************************/


	// load stored results
	use  "$N\4b.Analyze_results\Approach_a_attained_age/scenario`scen'_attained.dta" , replace
	
	sort givenage

	//Estimate bias for the CIF from each cause, for each time point we are interested over a series of ages in diagnosis
	foreach Z in CIF1_f${female}_at1 CIF1_f${female}_at2 CIF1_f${female}_at3 CIF1_f${female}_at4 CIF1_f${female}_at5 CIF1_f${female}_at6 CIF1_f${female}_at7 CIF1_f${female}_at8 CIF1_f${female}_at9 CIF1_f${female}_at10   ///
                 CIF2_f${female}_at1 CIF2_f${female}_at2 CIF2_f${female}_at3 CIF2_f${female}_at4 CIF2_f${female}_at5 CIF2_f${female}_at6 CIF2_f${female}_at7 CIF2_f${female}_at8 CIF2_f${female}_at9 CIF2_f${female}_at10 {
	
	forvalues i = 60(10)90 {
    
	// load truths
	include  "$N/2b.Truths_folder/scenario`scen'_Truth_`i'.do"
	
		//summarize estimations
	gen `Z'_`i'=`Z' if givenage==`i'
	summ `Z'_`i'
    local `Z'_`i'_t1t2= `r(mean)'

	gen `Z'_bias_`i'=.
	    
	replace `Z'_bias_`i'= `Z' - ${`Z'_true} if givenage==`i'
	summ `Z'_bias_`i'	
	}
	
	}
	
	
	//Calculate the coverages
	
	foreach X in CIF1_f${female}_at1 CIF1_f${female}_at2 CIF1_f${female}_at3 CIF1_f${female}_at4 CIF1_f${female}_at5 CIF1_f${female}_at6 CIF1_f${female}_at7 CIF1_f${female}_at8 CIF1_f${female}_at9 CIF1_f${female}_at10   ///
                 CIF2_f${female}_at1 CIF2_f${female}_at2 CIF2_f${female}_at3 CIF2_f${female}_at4 CIF2_f${female}_at5 CIF2_f${female}_at6 CIF2_f${female}_at7 CIF2_f${female}_at8 CIF2_f${female}_at9 CIF2_f${female}_at10 {
	  
	 forvalues i = 60(10)90 {
	 
	 // load truths

	 include  "$N\2b.Truths_folder/scenario`scen'_Truth_`i'.do"

	 //Generate an index variable of whether in the specific repetition, the 95% CI
	 //of the estimate includes the true value of the parameter
	 gen cov_index_`X'_`i'= .
	 replace cov_index_`X'_`i'=0 if givenage==`i' & canc_model_converged==1 & other_model_converged==1
	 replace cov_index_`X'_`i'=1 if `X'_LL<=${`X'_true} & `X'_UL>=${`X'_true} & givenage==`i'
	 
	 summ cov_index_`X'_`i' if cov_index_`X'_`i' != .
	 gen coverage_`X'_`i'= `r(mean)'
	 local coverage_`X'_`i'_t1t2= coverage_`X'_`i'
	}
	}
	
		// Derive the bias, convergence, relative bias and empirical standard errors 

	forvalues i = 60(10)90 {
	
		foreach var in  CIF1_f${female}_at1_bias_`i' CIF1_f${female}_at2_bias_`i' CIF1_f${female}_at3_bias_`i' CIF1_f${female}_at4_bias_`i' CIF1_f${female}_at5_bias_`i' ///
	                CIF1_f${female}_at6_bias_`i' CIF1_f${female}_at7_bias_`i' CIF1_f${female}_at8_bias_`i' CIF1_f${female}_at9_bias_`i' CIF1_f${female}_at10_bias_`i'   ///
                    CIF2_f${female}_at1_bias_`i' CIF2_f${female}_at2_bias_`i' CIF2_f${female}_at3_bias_`i' CIF2_f${female}_at4_bias_`i' CIF2_f${female}_at5_bias_`i' ///
					CIF2_f${female}_at6_bias_`i' CIF2_f${female}_at7_bias_`i' CIF2_f${female}_at8_bias_`i' CIF2_f${female}_at9_bias_`i' CIF2_f${female}_at10_bias_`i' {
	summ `var' 
	local b_`var'_t1t2= `r(mean)'
	local r_b_`var'_t1t2=`b_`var'_t1t2'
	local convergence_t1t2= `r(N)'

	}
	}
	
			//MCerrors and EmpSE
	
	
	foreach Z in CIF1_f${female}_at1 CIF1_f${female}_at2 CIF1_f${female}_at3 CIF1_f${female}_at4 CIF1_f${female}_at5 ///
				 CIF1_f${female}_at6 CIF1_f${female}_at7 CIF1_f${female}_at8 CIF1_f${female}_at9 CIF1_f${female}_at10   ///
                 CIF2_f${female}_at1 CIF2_f${female}_at2 CIF2_f${female}_at3 CIF2_f${female}_at4 CIF2_f${female}_at5  ///
				 CIF2_f${female}_at6 CIF2_f${female}_at7 CIF2_f${female}_at8 CIF2_f${female}_at9 CIF2_f${female}_at10 {
	
		
	forvalues i = 60(10)90 {
	
	include  "$N/2b.Truths_folder/scenario`scen'_Truth_`i'.do"
	
	gen sqdiff_`Z'_`i'=. 
	
	summ `Z' if givenage==`i'
	
	local mean_`Z'_`i'=`r(mean)'
	
	replace sqdiff_`Z'_`i'= (`Z'-`mean_`Z'_`i'')^2  if  givenage==`i'
	
	sum sqdiff_`Z'_`i' 
	
    local sum_sqdiff_`Z'_`i'= `r(mean)'*`r(N)'
	
    local MCer_`Z'_`i'_t1t2= sqrt((1/(`r(N)'*(`r(N)'-1)))*`sum_sqdiff_`Z'_`i'')
	
	local rb_MCer_`Z'_`i'_t1t2= (`MCer_`Z'_`i'_t1t2'/ ${`Z'_true} )*100
	
	local EmpSE_`Z'_`i'_t1t2= sqrt((1/(`r(N)'-1))*`sum_sqdiff_`Z'_`i'')
	
	//di `MCer_`Z'_`i'_t1t2'
	//di `EmpSE_`var'_t1t2'
	
	}
	
	}
	

/********************************/
	//Gender bias 
/********************************/

	foreach Z in hr_gender_at1 hr_gender_at2 hr_gender_at3 hr_gender_at4 hr_gender_at5    ///
                 hr_gender_at6 hr_gender_at7 hr_gender_at8 hr_gender_at9 hr_gender_at10 {
	
	forvalues i = 60(10)90 {
    
	// load truths
	include  "$N/2b.Truths_folder/scenario`scen'_Truth_`i'.do"

	gen `Z'_bias_`i'=.
	    
	replace `Z'_bias_`i'= `Z' - ${`Z'_true} if givenage==`i'
	summ `Z'_bias_`i'	
	}
	
	}
	
		forvalues i = 60(10)90 {
	
		foreach var in  hr_gender_at1_bias_`i' hr_gender_at2_bias_`i' hr_gender_at3_bias_`i' hr_gender_at4_bias_`i' hr_gender_at5_bias_`i' ///
	                    hr_gender_at6_bias_`i' hr_gender_at7_bias_`i' hr_gender_at8_bias_`i' hr_gender_at9_bias_`i' hr_gender_at10_bias_`i'    {
	summ `var' 
	local b_`var'_t1t2= `r(mean)'

	}
	}
	
	
	
/********************************/
	//Gender  true values
/********************************/
	forvalues i = 60(10)90 {
	foreach V in hr_gender_at1 hr_gender_at2 hr_gender_at3 hr_gender_at4 hr_gender_at5 /// 
			 hr_gender_at6 hr_gender_at7 hr_gender_at8 hr_gender_at9 hr_gender_at10 { 

	gen `V'_`i'= .
	replace `V'_`i'= `V' if givenage==`i'
				}
				 }
				
	forvalues i = 60(10)90 {
	 
	foreach var in hr_gender_at1_`i' hr_gender_at2_`i' hr_gender_at3_`i' hr_gender_at4_`i' hr_gender_at5_`i' /// 
			 hr_gender_at6_`i' hr_gender_at7_`i' hr_gender_at8_`i' hr_gender_at9_`i' hr_gender_at10_`i' { 
		summ `var' 
	local `var'_t1t2= `r(mean)'
	}
	}
	
	
	//CIF diffs, CIF ratios, total CIF ratio
		forvalues i = 60(10)90 {
				include  "$N/2b.Truths_folder/scenario`scen'_Truth_`i'.do"

	foreach V in    CIF1_diff_at1 CIF1_diff_at2 CIF1_diff_at3 CIF1_diff_at4 CIF1_diff_at5 /// 
				    CIF1_diff_at6 CIF1_diff_at7 CIF1_diff_at8 CIF1_diff_at9 CIF1_diff_at10 ///
					CIF2_diff_at1 CIF2_diff_at2 CIF2_diff_at3 CIF2_diff_at4 CIF2_diff_at5 ///
					CIF2_diff_at6 CIF2_diff_at7 CIF2_diff_at8 CIF2_diff_at9 CIF2_diff_at10 ///
					CIF1_ratio_at1 CIF1_ratio_at2 CIF1_ratio_at3 CIF1_ratio_at4 CIF1_ratio_at5 ///
					CIF1_ratio_at6 CIF1_ratio_at7 CIF1_ratio_at8 CIF1_ratio_at9 CIF1_ratio_at10 ///
					CIF2_ratio_at1 CIF2_ratio_at2 CIF2_ratio_at3 CIF2_ratio_at4 CIF2_ratio_at5 ///
					CIF2_ratio_at6 CIF2_ratio_at7 CIF2_ratio_at8 CIF2_ratio_at9 CIF2_ratio_at10 ///
					CIF_total_at1 CIF_total_at2 CIF_total_at3 CIF_total_at4 CIF_total_at5 ///
					CIF_total_at6 CIF_total_at7 CIF_total_at8 CIF_total_at9 CIF_total_at10  ///
					  {
 
			 
	gen `V'_`i'= .

	replace `V'_`i'= `V' if givenage==`i'
				}
				 }
				
	forvalues i = 60(10)90 {
				include  "$N/2b.Truths_folder/scenario`scen'_Truth_`i'.do"

	 
	foreach var in 	CIF1_diff_at1_`i' CIF1_diff_at2_`i' CIF1_diff_at3_`i' CIF1_diff_at4_`i' CIF1_diff_at5_`i' /// 
				    CIF1_diff_at6_`i' CIF1_diff_at7_`i' CIF1_diff_at8_`i' CIF1_diff_at9_`i' CIF1_diff_at10_`i' ///
					CIF2_diff_at1_`i' CIF2_diff_at2_`i' CIF2_diff_at3_`i' CIF2_diff_at4_`i' CIF2_diff_at5_`i' ///
					CIF2_diff_at6_`i' CIF2_diff_at7_`i' CIF2_diff_at8_`i' CIF2_diff_at9_`i' CIF2_diff_at10_`i' ///
					CIF1_ratio_at1_`i' CIF1_ratio_at2_`i' CIF1_ratio_at3_`i' CIF1_ratio_at4_`i' CIF1_ratio_at5_`i' CIF1_ratio_at6_`i' ///
				    CIF1_ratio_at7_`i' CIF1_ratio_at8_`i' CIF1_ratio_at9_`i' CIF1_ratio_at10_`i' ///
				    CIF2_ratio_at1_`i' CIF2_ratio_at2_`i' CIF2_ratio_at3_`i' CIF2_ratio_at4_`i' CIF2_ratio_at5_`i' ///
				    CIF2_ratio_at6_`i' CIF2_ratio_at7_`i' CIF2_ratio_at8_`i' CIF2_ratio_at9_`i' CIF2_ratio_at10_`i' ///
				    CIF_total_at1_`i' CIF_total_at2_`i' CIF_total_at3_`i' CIF_total_at4_`i' CIF_total_at5_`i' ///
				    CIF_total_at6_`i' CIF_total_at7_`i' CIF_total_at8_`i' CIF_total_at9_`i' CIF_total_at10_`i'  ///
					  {
 
	summ `var' 
	local `var'_t1t2= `r(mean)'
 
	}
	}
	

//Create a file with the performance measures of the simulation over scenarios,
// the 2 genders and different ages at diagnosis
postutil clear
postfile comparison     CIF1_f${female}_at1_t1t2 CIF1_f${female}_at2_t1t2 CIF1_f${female}_at3_t1t2 CIF1_f${female}_at4_t1t2 CIF1_f${female}_at5_t1t2 CIF1_f${female}_at6_t1t2 ///
						CIF1_f${female}_at7_t1t2 CIF1_f${female}_at8_t1t2 CIF1_f${female}_at9_t1t2 CIF1_f${female}_at10_t1t2  ///
						///
						CIF2_f${female}_at1_t1t2 CIF2_f${female}_at2_t1t2 CIF2_f${female}_at3_t1t2 CIF2_f${female}_at4_t1t2 CIF2_f${female}_at5_t1t2 CIF2_f${female}_at6_t1t2 ///
						CIF2_f${female}_at7_t1t2 CIF2_f${female}_at8_t1t2 CIF2_f${female}_at9_t1t2 CIF2_f${female}_at10_t1t2  ///
						///
						CIF1_f${female}_at1_t1t1 CIF1_f${female}_at2_t1t1 CIF1_f${female}_at3_t1t1 CIF1_f${female}_at4_t1t1 CIF1_f${female}_at5_t1t1 CIF1_f${female}_at6_t1t1 ///
						CIF1_f${female}_at7_t1t1 CIF1_f${female}_at8_t1t1 CIF1_f${female}_at9_t1t1 CIF1_f${female}_at10_t1t1  ///
						///
						CIF2_f${female}_at1_t1t1 CIF2_f${female}_at2_t1t1 CIF2_f${female}_at3_t1t1 CIF2_f${female}_at4_t1t1 CIF2_f${female}_at5_t1t1 CIF2_f${female}_at6_t1t1 ///
						CIF2_f${female}_at7_t1t1 CIF2_f${female}_at8_t1t1 CIF2_f${female}_at9_t1t1 CIF2_f${female}_at10_t1t1  ///
						///
						bias_CIF1_f${female}_at1_t1t2 bias_CIF1_f${female}_at2_t1t2 bias_CIF1_f${female}_at3_t1t2 bias_CIF1_f${female}_at4_t1t2 bias_CIF1_f${female}_at5_t1t2 bias_CIF1_f${female}_at6_t1t2 ///
						bias_CIF1_f${female}_at7_t1t2 bias_CIF1_f${female}_at8_t1t2 bias_CIF1_f${female}_at9_t1t2 bias_CIF1_f${female}_at10_t1t2  ///
						///
					    bias_CIF2_f${female}_at1_t1t2 bias_CIF2_f${female}_at2_t1t2 bias_CIF2_f${female}_at3_t1t2 bias_CIF2_f${female}_at4_t1t2 bias_CIF2_f${female}_at5_t1t2 bias_CIF2_f${female}_at6_t1t2 ///
						bias_CIF2_f${female}_at7_t1t2 bias_CIF2_f${female}_at8_t1t2 bias_CIF2_f${female}_at9_t1t2 bias_CIF2_f${female}_at10_t1t2 ///
					    ///
						bias_CIF1_f${female}_at1_t1t1 bias_CIF1_f${female}_at2_t1t1 bias_CIF1_f${female}_at3_t1t1 bias_CIF1_f${female}_at4_t1t1 bias_CIF1_f${female}_at5_t1t1 bias_CIF1_f${female}_at6_t1t1 ///
						bias_CIF1_f${female}_at7_t1t1 bias_CIF1_f${female}_at8_t1t1 bias_CIF1_f${female}_at9_t1t1 bias_CIF1_f${female}_at10_t1t1  ///
						///
					    bias_CIF2_f${female}_at1_t1t1 bias_CIF2_f${female}_at2_t1t1 bias_CIF2_f${female}_at3_t1t1 bias_CIF2_f${female}_at4_t1t1 bias_CIF2_f${female}_at5_t1t1 bias_CIF2_f${female}_at6_t1t1 ///
						bias_CIF2_f${female}_at7_t1t1 bias_CIF2_f${female}_at8_t1t1 bias_CIF2_f${female}_at9_t1t1 bias_CIF2_f${female}_at10_t1t1 ///
						///
						MCer_CIF1_f${female}_at1_t1t2 MCer_CIF1_f${female}_at2_t1t2 MCer_CIF1_f${female}_at3_t1t2 MCer_CIF1_f${female}_at4_t1t2 MCer_CIF1_f${female}_at5_t1t2 MCer_CIF1_f${female}_at6_t1t2 ///
						MCer_CIF1_f${female}_at7_t1t2 MCer_CIF1_f${female}_at8_t1t2 MCer_CIF1_f${female}_at9_t1t2 MCer_CIF1_f${female}_at10_t1t2  ///
						///
					    MCer_CIF2_f${female}_at1_t1t2 MCer_CIF2_f${female}_at2_t1t2 MCer_CIF2_f${female}_at3_t1t2 MCer_CIF2_f${female}_at4_t1t2 MCer_CIF2_f${female}_at5_t1t2 MCer_CIF2_f${female}_at6_t1t2 ///
						MCer_CIF2_f${female}_at7_t1t2 MCer_CIF2_f${female}_at8_t1t2 MCer_CIF2_f${female}_at9_t1t2 MCer_CIF2_f${female}_at10_t1t2 ///
					    ///
						MCer_CIF1_f${female}_at1_t1t1 MCer_CIF1_f${female}_at2_t1t1 MCer_CIF1_f${female}_at3_t1t1 MCer_CIF1_f${female}_at4_t1t1 MCer_CIF1_f${female}_at5_t1t1 MCer_CIF1_f${female}_at6_t1t1 ///
						MCer_CIF1_f${female}_at7_t1t1 MCer_CIF1_f${female}_at8_t1t1 MCer_CIF1_f${female}_at9_t1t1 MCer_CIF1_f${female}_at10_t1t1  ///
						///
					    MCer_CIF2_f${female}_at1_t1t1 MCer_CIF2_f${female}_at2_t1t1 MCer_CIF2_f${female}_at3_t1t1 MCer_CIF2_f${female}_at4_t1t1 MCer_CIF2_f${female}_at5_t1t1 MCer_CIF2_f${female}_at6_t1t1 ///
						MCer_CIF2_f${female}_at7_t1t1 MCer_CIF2_f${female}_at8_t1t1 MCer_CIF2_f${female}_at9_t1t1 MCer_CIF2_f${female}_at10_t1t1 ///
						///
						rb_MCer_CIF1_f${female}_at1_t1t2 rb_MCer_CIF1_f${female}_at2_t1t2 rb_MCer_CIF1_f${female}_at3_t1t2 rb_MCer_CIF1_f${female}_at4_t1t2 rb_MCer_CIF1_f${female}_at5_t1t2 rb_MCer_CIF1_f${female}_at6_t1t2 ///
						rb_MCer_CIF1_f${female}_at7_t1t2 rb_MCer_CIF1_f${female}_at8_t1t2 rb_MCer_CIF1_f${female}_at9_t1t2 rb_MCer_CIF1_f${female}_at10_t1t2  ///
						///
					    rb_MCer_CIF2_f${female}_at1_t1t2 rb_MCer_CIF2_f${female}_at2_t1t2 rb_MCer_CIF2_f${female}_at3_t1t2 rb_MCer_CIF2_f${female}_at4_t1t2 rb_MCer_CIF2_f${female}_at5_t1t2 rb_MCer_CIF2_f${female}_at6_t1t2 ///
						rb_MCer_CIF2_f${female}_at7_t1t2 rb_MCer_CIF2_f${female}_at8_t1t2 rb_MCer_CIF2_f${female}_at9_t1t2 rb_MCer_CIF2_f${female}_at10_t1t2 ///
					    ///
						rb_MCer_CIF1_f${female}_at1_t1t1 rb_MCer_CIF1_f${female}_at2_t1t1 rb_MCer_CIF1_f${female}_at3_t1t1 rb_MCer_CIF1_f${female}_at4_t1t1 rb_MCer_CIF1_f${female}_at5_t1t1 rb_MCer_CIF1_f${female}_at6_t1t1 ///
						rb_MCer_CIF1_f${female}_at7_t1t1 rb_MCer_CIF1_f${female}_at8_t1t1 rb_MCer_CIF1_f${female}_at9_t1t1 rb_MCer_CIF1_f${female}_at10_t1t1  ///
						///
					    rb_MCer_CIF2_f${female}_at1_t1t1 rb_MCer_CIF2_f${female}_at2_t1t1 rb_MCer_CIF2_f${female}_at3_t1t1 rb_MCer_CIF2_f${female}_at4_t1t1 rb_MCer_CIF2_f${female}_at5_t1t1 rb_MCer_CIF2_f${female}_at6_t1t1 ///
						rb_MCer_CIF2_f${female}_at7_t1t1 rb_MCer_CIF2_f${female}_at8_t1t1 rb_MCer_CIF2_f${female}_at9_t1t1 rb_MCer_CIF2_f${female}_at10_t1t1 ///
						///
						rel_bias_CIF1_f${female}_at1_t1t1  rel_bias_CIF1_f${female}_at2_t1t1  rel_bias_CIF1_f${female}_at3_t1t1 rel_bias_CIF1_f${female}_at4_t1t1  rel_bias_CIF1_f${female}_at5_t1t1  rel_bias_CIF1_f${female}_at6_t1t1   ///
						rel_bias_CIF1_f${female}_at7_t1t1  rel_bias_CIF1_f${female}_at8_t1t1  rel_bias_CIF1_f${female}_at9_t1t1 rel_bias_CIF1_f${female}_at10_t1t1 ///
						///
						rel_bias_CIF2_f${female}_at1_t1t1  rel_bias_CIF2_f${female}_at2_t1t1  rel_bias_CIF2_f${female}_at3_t1t1 rel_bias_CIF2_f${female}_at4_t1t1  rel_bias_CIF2_f${female}_at5_t1t1  rel_bias_CIF2_f${female}_at6_t1t1   ///
						rel_bias_CIF2_f${female}_at7_t1t1  rel_bias_CIF2_f${female}_at8_t1t1  rel_bias_CIF2_f${female}_at9_t1t1 rel_bias_CIF2_f${female}_at10_t1t1 ///
						///
						rel_bias_CIF1_f${female}_at1_t1t2  rel_bias_CIF1_f${female}_at2_t1t2  rel_bias_CIF1_f${female}_at3_t1t2 rel_bias_CIF1_f${female}_at4_t1t2  rel_bias_CIF1_f${female}_at5_t1t2  rel_bias_CIF1_f${female}_at6_t1t2  ///
						rel_bias_CIF1_f${female}_at7_t1t2  rel_bias_CIF1_f${female}_at8_t1t2  rel_bias_CIF1_f${female}_at9_t1t2 rel_bias_CIF1_f${female}_at10_t1t2  ///
						///
						rel_bias_CIF2_f${female}_at1_t1t2  rel_bias_CIF2_f${female}_at2_t1t2  rel_bias_CIF2_f${female}_at3_t1t2 rel_bias_CIF2_f${female}_at4_t1t2  rel_bias_CIF2_f${female}_at5_t1t2  rel_bias_CIF2_f${female}_at6_t1t2  ///
						rel_bias_CIF2_f${female}_at7_t1t2  rel_bias_CIF2_f${female}_at8_t1t2  rel_bias_CIF2_f${female}_at9_t1t2 rel_bias_CIF2_f${female}_at10_t1t2  ///
						///
						EmpSE_CIF1_f${female}_at1_t1t1  EmpSE_CIF1_f${female}_at2_t1t1  EmpSE_CIF1_f${female}_at3_t1t1  EmpSE_CIF1_f${female}_at4_t1t1  EmpSE_CIF1_f${female}_at5_t1t1  EmpSE_CIF1_f${female}_at6_t1t1 ///
						EmpSE_CIF1_f${female}_at7_t1t1  EmpSE_CIF1_f${female}_at8_t1t1  EmpSE_CIF1_f${female}_at9_t1t1 EmpSE_CIF1_f${female}_at10_t1t1 ///
						///
						EmpSE_CIF2_f${female}_at1_t1t1  EmpSE_CIF2_f${female}_at2_t1t1  EmpSE_CIF2_f${female}_at3_t1t1  EmpSE_CIF2_f${female}_at4_t1t1  EmpSE_CIF2_f${female}_at5_t1t1  EmpSE_CIF2_f${female}_at6_t1t1 ///
						EmpSE_CIF2_f${female}_at7_t1t1  EmpSE_CIF2_f${female}_at8_t1t1  EmpSE_CIF2_f${female}_at9_t1t1 EmpSE_CIF2_f${female}_at10_t1t1 ///
						///
						EmpSE_CIF1_f${female}_at1_t1t2  EmpSE_CIF1_f${female}_at2_t1t2  EmpSE_CIF1_f${female}_at3_t1t2 EmpSE_CIF1_f${female}_at4_t1t2  EmpSE_CIF1_f${female}_at5_t1t2  EmpSE_CIF1_f${female}_at6_t1t2   ///
						EmpSE_CIF1_f${female}_at7_t1t2  EmpSE_CIF1_f${female}_at8_t1t2  EmpSE_CIF1_f${female}_at9_t1t2  EmpSE_CIF1_f${female}_at10_t1t2 ///
						///
						EmpSE_CIF2_f${female}_at1_t1t2  EmpSE_CIF2_f${female}_at2_t1t2  EmpSE_CIF2_f${female}_at3_t1t2 EmpSE_CIF2_f${female}_at4_t1t2  EmpSE_CIF2_f${female}_at5_t1t2  EmpSE_CIF2_f${female}_at6_t1t2   ///
						EmpSE_CIF2_f${female}_at7_t1t2  EmpSE_CIF2_f${female}_at8_t1t2  EmpSE_CIF2_f${female}_at9_t1t2  EmpSE_CIF2_f${female}_at10_t1t2 ///
						///
						Rel_prec_CIF1_f${female}_at1  Rel_prec_CIF1_f${female}_at2  Rel_prec_CIF1_f${female}_at3 Rel_prec_CIF1_f${female}_at4  Rel_prec_CIF1_f${female}_at5  Rel_prec_CIF1_f${female}_at6  ///
						Rel_prec_CIF1_f${female}_at7  Rel_prec_CIF1_f${female}_at8  Rel_prec_CIF1_f${female}_at9 Rel_prec_CIF1_f${female}_at10 ///
						///
						Rel_prec_CIF2_f${female}_at1  Rel_prec_CIF2_f${female}_at2  Rel_prec_CIF2_f${female}_at3 Rel_prec_CIF2_f${female}_at4  Rel_prec_CIF2_f${female}_at5  Rel_prec_CIF2_f${female}_at6  ///
						Rel_prec_CIF2_f${female}_at7  Rel_prec_CIF2_f${female}_at8  Rel_prec_CIF2_f${female}_at9 Rel_prec_CIF2_f${female}_at10 ///
						///
						Coverage_CIF1_f${female}_at1_t1t1  Coverage_CIF1_f${female}_at2_t1t1  Coverage_CIF1_f${female}_at3_t1t1 Coverage_CIF1_f${female}_at4_t1t1  Coverage_CIF1_f${female}_at5_t1t1  Coverage_CIF1_f${female}_at6_t1t1 ///
						Coverage_CIF1_f${female}_at7_t1t1  Coverage_CIF1_f${female}_at8_t1t1  Coverage_CIF1_f${female}_at9_t1t1 Coverage_CIF1_f${female}_at10_t1t1  ///
						///
						Coverage_CIF2_f${female}_at1_t1t1  Coverage_CIF2_f${female}_at2_t1t1  Coverage_CIF2_f${female}_at3_t1t1 Coverage_CIF2_f${female}_at4_t1t1  Coverage_CIF2_f${female}_at5_t1t1  Coverage_CIF2_f${female}_at6_t1t1 ///
						Coverage_CIF2_f${female}_at7_t1t1  Coverage_CIF2_f${female}_at8_t1t1  Coverage_CIF2_f${female}_at9_t1t1 Coverage_CIF2_f${female}_at10_t1t1  ///
						///
						Coverage_CIF1_f${female}_at1_t1t2  Coverage_CIF1_f${female}_at2_t1t2  Coverage_CIF1_f${female}_at3_t1t2  Coverage_CIF1_f${female}_at4_t1t2  Coverage_CIF1_f${female}_at5_t1t2  Coverage_CIF1_f${female}_at6_t1t2 ///
						Coverage_CIF1_f${female}_at7_t1t2  Coverage_CIF1_f${female}_at8_t1t2  Coverage_CIF1_f${female}_at9_t1t2 Coverage_CIF1_f${female}_at10_t1t2 ///
						///
				        Coverage_CIF2_f${female}_at1_t1t2  Coverage_CIF2_f${female}_at2_t1t2  Coverage_CIF2_f${female}_at3_t1t2  Coverage_CIF2_f${female}_at4_t1t2  Coverage_CIF2_f${female}_at5_t1t2  Coverage_CIF2_f${female}_at6_t1t2 ///
						Coverage_CIF2_f${female}_at7_t1t2  Coverage_CIF2_f${female}_at8_t1t2  Coverage_CIF2_f${female}_at9_t1t2 Coverage_CIF2_f${female}_at10_t1t2 ///
						///
						hr_gender_at1_t1t2 hr_gender_at2_t1t2 hr_gender_at3_t1t2 hr_gender_at4_t1t2 hr_gender_at5_t1t2  ///
						hr_gender_at6_t1t2 hr_gender_at7_t1t2 hr_gender_at8_t1t2 hr_gender_at9_t1t2 hr_gender_at10_t1t2 ///
						///
						hr_gender_at1_t1t1 hr_gender_at2_t1t1 hr_gender_at3_t1t1 hr_gender_at4_t1t1 hr_gender_at5_t1t1 ///
						hr_gender_at6_t1t1 hr_gender_at7_t1t1 hr_gender_at8_t1t1 hr_gender_at9_t1t1 hr_gender_at10_t1t1 ///
						///
						bias_hr_gender_at1_t1t2 bias_hr_gender_at2_t1t2 bias_hr_gender_at3_t1t2 bias_hr_gender_at4_t1t2 bias_hr_gender_at5_t1t2  ///
						bias_hr_gender_at6_t1t2 bias_hr_gender_at7_t1t2 bias_hr_gender_at8_t1t2 bias_hr_gender_at9_t1t2 bias_hr_gender_at10_t1t2 ///
						///
						bias_hr_gender_at1_t1t1 bias_hr_gender_at2_t1t1 bias_hr_gender_at3_t1t1 bias_hr_gender_at4_t1t1 bias_hr_gender_at5_t1t1 ///
						bias_hr_gender_at6_t1t1 bias_hr_gender_at7_t1t1 bias_hr_gender_at8_t1t1 bias_hr_gender_at9_t1t1 bias_hr_gender_at10_t1t1 ///
						///
						///
						CIF1_diff_at1_t1t2 CIF1_diff_at2_t1t2 CIF1_diff_at3_t1t2 CIF1_diff_at4_t1t2 CIF1_diff_at5_t1t2 CIF1_diff_at6_t1t2 ///
						CIF1_diff_at7_t1t2 CIF1_diff_at8_t1t2 CIF1_diff_at9_t1t2 CIF1_diff_at10_t1t2  ///
						///
			          	CIF1_diff_at1_t1t1 CIF1_diff_at2_t1t1 CIF1_diff_at3_t1t1 CIF1_diff_at4_t1t1 CIF1_diff_at5_t1t1 /// 
						CIF1_diff_at6_t1t1 CIF1_diff_at7_t1t1 CIF1_diff_at8_t1t1 CIF1_diff_at9_t1t1 CIF1_diff_at10_t1t1 ///
						///
					    CIF2_diff_at1_t1t2 CIF2_diff_at2_t1t2 CIF2_diff_at3_t1t2 CIF2_diff_at4_t1t2 CIF2_diff_at5_t1t2 ///
						CIF2_diff_at6_t1t2 CIF2_diff_at7_t1t2 CIF2_diff_at8_t1t2 CIF2_diff_at9_t1t2 CIF2_diff_at10_t1t2 ///
						///
					    CIF2_diff_at1_t1t1 CIF2_diff_at2_t1t1 CIF2_diff_at3_t1t1 CIF2_diff_at4_t1t1 CIF2_diff_at5_t1t1 ///
						CIF2_diff_at6_t1t1 CIF2_diff_at7_t1t1 CIF2_diff_at8_t1t1 CIF2_diff_at9_t1t1 CIF2_diff_at10_t1t1 ///
						///
						CIF1_ratio_at1_t1t2 CIF1_ratio_at2_t1t2 CIF1_ratio_at3_t1t2 CIF1_ratio_at4_t1t2 CIF1_ratio_at5_t1t2 CIF1_ratio_at6_t1t2 ///
						CIF1_ratio_at7_t1t2 CIF1_ratio_at8_t1t2 CIF1_ratio_at9_t1t2 CIF1_ratio_at10_t1t2  ///
						///
			          	CIF1_ratio_at1_t1t1 CIF1_ratio_at2_t1t1 CIF1_ratio_at3_t1t1 CIF1_ratio_at4_t1t1 CIF1_ratio_at5_t1t1 /// 
						CIF1_ratio_at6_t1t1 CIF1_ratio_at7_t1t1 CIF1_ratio_at8_t1t1 CIF1_ratio_at9_t1t1 CIF1_ratio_at10_t1t1 ///
						///
					    CIF2_ratio_at1_t1t2 CIF2_ratio_at2_t1t2 CIF2_ratio_at3_t1t2 CIF2_ratio_at4_t1t2 CIF2_ratio_at5_t1t2 ///
						CIF2_ratio_at6_t1t2 CIF2_ratio_at7_t1t2 CIF2_ratio_at8_t1t2 CIF2_ratio_at9_t1t2 CIF2_ratio_at10_t1t2 ///
						///
					    CIF2_ratio_at1_t1t1 CIF2_ratio_at2_t1t1 CIF2_ratio_at3_t1t1 CIF2_ratio_at4_t1t1 CIF2_ratio_at5_t1t1 ///
						CIF2_ratio_at6_t1t1 CIF2_ratio_at7_t1t1 CIF2_ratio_at8_t1t1 CIF2_ratio_at9_t1t1 CIF2_ratio_at10_t1t1 ///
						///
						CIF_total_at1_t1t2 CIF_total_at2_t1t2 CIF_total_at3_t1t2 CIF_total_at4_t1t2 CIF_total_at5_t1t2 ///
						CIF_total_at6_t1t2 CIF_total_at7_t1t2 CIF_total_at8_t1t2 CIF_total_at9_t1t2 CIF_total_at10_t1t2  ///
						///
					  	CIF_total_at1_t1t1 CIF_total_at2_t1t1 CIF_total_at3_t1t1 CIF_total_at4_t1t1 CIF_total_at5_t1t1 ///
						CIF_total_at6_t1t1 CIF_total_at7_t1t1 CIF_total_at8_t1t1 CIF_total_at9_t1t1 CIF_total_at10_t1t1  ///
						///
						CIF1_f${female}_at1_true CIF1_f${female}_at2_true CIF1_f${female}_at3_true CIF1_f${female}_at4_true CIF1_f${female}_at5_true ///
						CIF1_f${female}_at6_true CIF1_f${female}_at7_true CIF1_f${female}_at8_true CIF1_f${female}_at9_true CIF1_f${female}_at10_true ///
						///
						CIF2_f${female}_at1_true CIF2_f${female}_at2_true CIF2_f${female}_at3_true CIF2_f${female}_at4_true CIF2_f${female}_at5_true ///
						CIF2_f${female}_at6_true CIF2_f${female}_at7_true CIF2_f${female}_at8_true CIF2_f${female}_at9_true CIF2_f${female}_at10_true ///
						///
						CIF1_diff_at1_true 	CIF1_diff_at2_true CIF1_diff_at3_true CIF1_diff_at4_true CIF1_diff_at5_true ///
						CIF1_diff_at6_true  CIF1_diff_at7_true CIF1_diff_at8_true CIF1_diff_at9_true CIF1_diff_at10_true  ///
						///
						CIF2_diff_at1_true CIF2_diff_at2_true CIF2_diff_at3_true CIF2_diff_at4_true CIF2_diff_at5_true ///
						CIF2_diff_at6_true CIF2_diff_at7_true CIF2_diff_at8_true CIF2_diff_at9_true CIF2_diff_at10_true ///
						///
						CIF1_ratio_at1_true CIF1_ratio_at2_true CIF1_ratio_at3_true CIF1_ratio_at4_true CIF1_ratio_at5_true ///
						CIF1_ratio_at6_true CIF1_ratio_at7_true CIF1_ratio_at8_true CIF1_ratio_at9_true CIF1_ratio_at10_true ///
						///
						CIF2_ratio_at1_true CIF2_ratio_at2_true CIF2_ratio_at3_true CIF2_ratio_at4_true CIF2_ratio_at5_true  ///
						CIF2_ratio_at6_true CIF2_ratio_at7_true CIF2_ratio_at8_true CIF2_ratio_at9_true CIF2_ratio_at10_true ///
						///
						CIF_total_at1_true CIF_total_at2_true CIF_total_at3_true CIF_total_at4_true CIF_total_at5_true ///
						CIF_total_at6_true CIF_total_at7_true CIF_total_at8_true CIF_total_at9_true CIF_total_at10_true ///
						///
						hr_gender_at1_true hr_gender_at2_true hr_gender_at3_true hr_gender_at4_true hr_gender_at5_true ///
						hr_gender_at6_true hr_gender_at7_true hr_gender_at8_true hr_gender_at9_true hr_gender_at10_true ///
						///
						obs givenage scenario model independence female convergence_t1t1 convergence_t1t2 sdage str20 prop_sex_other str20 base_function ///
			            using ///
						"$N\5b.Performance_results\compare_splines_int\scen_`scen'_compare_splines_int_f`female'.dta", replace



   
	
forvalues i = 60(10)90 {



	include  "$N/2b.Truths_folder/scenario`scen'_Truth_`i'.do"
	
	foreach k in CIF1_f${female}_at1 CIF1_f${female}_at2 CIF1_f${female}_at3 CIF1_f${female}_at4 CIF1_f${female}_at5 CIF1_f${female}_at6 CIF1_f${female}_at7 CIF1_f${female}_at8 CIF1_f${female}_at9 CIF1_f${female}_at10   ///
                 CIF2_f${female}_at1 CIF2_f${female}_at2 CIF2_f${female}_at3 CIF2_f${female}_at4 CIF2_f${female}_at5 CIF2_f${female}_at6 CIF2_f${female}_at7 CIF2_f${female}_at8 CIF2_f${female}_at9 CIF2_f${female}_at10 { 
	

local `k'_t1t2= ``k'_`i'_t1t2'	

local `k'_t1t1= ``k'_`i'_t1t1'	

local bias_`k'_t1t2= `r_b_`k'_bias_`i'_t1t2'     //Bias in approach a-Attained

local bias_`k'_t1t1= `r_b_`k'_bias_`i'_t1t1'	   //Bias in approach d-Splines/Int

local MCer_`k'_t1t2= `MCer_`k'_`i'_t1t2'      //Monte carlo errors in approach a-Attained

local MCer_`k'_t1t1= `MCer_`k'_`i'_t1t1'	 	 //Monte carlo errors in approach d-Splines/Int


local rb_MCer_`k'_t1t2= `rb_MCer_`k'_`i'_t1t2'
                                 
local rb_MCer_`k'_t1t1= `rb_MCer_`k'_`i'_t1t1'
	
local bias_dif_`k'= .

local rel_bias_`k'_t1t2= (`b_`k'_bias_`i'_t1t2'/${`k'_true} ) *100  //Relative bias in approach a-Attained

local rel_bias_`k'_t1t1= (`b_`k'_bias_`i'_t1t1'/${`k'_true} ) *100  //Relative bias in approach d-Splines/Int

local EmpSE_`k'_t1t1= `EmpSE_`k'_`i'_t1t1'      //Empirical st.errors in approach a-Attained

local EmpSE_`k'_t1t2= `EmpSE_`k'_`i'_t1t2'       //Empirical st.errors in approach d-Splines/Int

local Rel_prec_`k'=100* ((`EmpSE_`k'_`i'_t1t2'/`EmpSE_`k'_`i'_t1t1')^2-1)    //Relative precision approach d-Splines/Int

local Coverage_`k'_t1t2= `coverage_`k'_`i'_t1t2'     //Coverage of approach a-Attained

local Coverage_`k'_t1t1= `coverage_`k'_`i'_t1t1'      //Coverage of approach d-Splines/Int

local givenage= `i'

local scenario= `scen'

local model =3

local convergence_t1t1=`convergence_t1t1'         //Convergence of approach d-Splines/Int
local convergence_t1t2=`convergence_t1t2'         //Convergence of approach a-Attained
local sdage= ${sdage}
local prop_sex_other= "$prop_sex_other"

local base_function= "$base_function"

local obs= $obs


local independence=$independence

local female=$female


}

//Local values of the estimations of HR of gender for other cause mortality

foreach o in hr_gender_at1 hr_gender_at2 hr_gender_at3 hr_gender_at4 hr_gender_at5 /// 
			 hr_gender_at6 hr_gender_at7 hr_gender_at8 hr_gender_at9 hr_gender_at10 { 
				 
	local  bias_`o'_t1t2=	`b_`o'_bias_`i'_t1t2'
	local  bias_`o'_t1t1=	`b_`o'_bias_`i'_t1t1'
	
	local  `o'_t1t2= 	``o'_`i'_t1t2'
	local  `o'_t1t1= 	``o'_`i'_t1t1'

				 }
				 
//Local values of the estimations of CIF ratios and differences between males and females for other cause mortality

				 
	foreach h in CIF1_diff_at1 CIF1_diff_at2 CIF1_diff_at3 CIF1_diff_at4 CIF1_diff_at5 /// 
			  CIF1_diff_at6 CIF1_diff_at7 CIF1_diff_at8 CIF1_diff_at9 CIF1_diff_at10 ///
			  CIF2_diff_at1 CIF2_diff_at2 CIF2_diff_at3 CIF2_diff_at4 CIF2_diff_at5 ///
		      CIF2_diff_at6 CIF2_diff_at7 CIF2_diff_at8 CIF2_diff_at9 CIF2_diff_at10 ///
			  CIF1_ratio_at1 CIF1_ratio_at2 CIF1_ratio_at3 CIF1_ratio_at4 CIF1_ratio_at5 ///
			  CIF1_ratio_at6 CIF1_ratio_at7 CIF1_ratio_at8 CIF1_ratio_at9 CIF1_ratio_at10 ///
			  CIF2_ratio_at1 CIF2_ratio_at2 CIF2_ratio_at3 CIF2_ratio_at4 CIF2_ratio_at5 ///
			  CIF2_ratio_at6 CIF2_ratio_at7 CIF2_ratio_at8 CIF2_ratio_at9 CIF2_ratio_at10 ///
		      CIF_total_at1 CIF_total_at2 CIF_total_at3 CIF_total_at4 CIF_total_at5 ///
			  CIF_total_at6 CIF_total_at7 CIF_total_at8 CIF_total_at9 CIF_total_at10  {
			  
 		
	local  `h'_t1t2= 	``h'_`i'_t1t2' 
	local  `h'_t1t1= 	``h'_`i'_t1t1'	

	local `h'_true= ${`h'_true}

}

//Local values of the estimations of true CIF values

foreach h in  CIF1_f${female}_at1 CIF1_f${female}_at2 CIF1_f${female}_at3 CIF1_f${female}_at4 CIF1_f${female}_at5  ///
              CIF1_f${female}_at6 CIF1_f${female}_at7 CIF1_f${female}_at8 CIF1_f${female}_at9 CIF1_f${female}_at10  ///
              CIF2_f${female}_at1 CIF2_f${female}_at2 CIF2_f${female}_at3 CIF2_f${female}_at4 CIF2_f${female}_at5	 ///	
              CIF2_f${female}_at6 CIF2_f${female}_at7 CIF2_f${female}_at8 CIF2_f${female}_at9 CIF2_f${female}_at10	  ///
			   {
			  
	local `h'_true= ${`h'_true}

}

//Local values of the estimations of true HR values of gender for other cause mortality

foreach y in hr_gender_at1 hr_gender_at2 hr_gender_at3 hr_gender_at4 hr_gender_at5 /// 
			 hr_gender_at6 hr_gender_at7 hr_gender_at8 hr_gender_at9 hr_gender_at10 { 
			 
	local `y'_true= ${`y'_true}
			 }




	post comparison     (`CIF1_f${female}_at1_t1t2') (`CIF1_f${female}_at2_t1t2') (`CIF1_f${female}_at3_t1t2') (`CIF1_f${female}_at4_t1t2') (`CIF1_f${female}_at5_t1t2') (`CIF1_f${female}_at6_t1t2') ///
						(`CIF1_f${female}_at7_t1t2') (`CIF1_f${female}_at8_t1t2') (`CIF1_f${female}_at9_t1t2') (`CIF1_f${female}_at10_t1t2')  ///
						///
						(`CIF2_f${female}_at1_t1t2') (`CIF2_f${female}_at2_t1t2') (`CIF2_f${female}_at3_t1t2') (`CIF2_f${female}_at4_t1t2') (`CIF2_f${female}_at5_t1t2') (`CIF2_f${female}_at6_t1t2') ///
						(`CIF2_f${female}_at7_t1t2') (`CIF2_f${female}_at8_t1t2') (`CIF2_f${female}_at9_t1t2') (`CIF2_f${female}_at10_t1t2')  ///
						///
						(`CIF1_f${female}_at1_t1t1') (`CIF1_f${female}_at2_t1t1') (`CIF1_f${female}_at3_t1t1') (`CIF1_f${female}_at4_t1t1') (`CIF1_f${female}_at5_t1t1') (`CIF1_f${female}_at6_t1t1') ///
						(`CIF1_f${female}_at7_t1t1') (`CIF1_f${female}_at8_t1t1') (`CIF1_f${female}_at9_t1t1') (`CIF1_f${female}_at10_t1t1')  ///
						///
						(`CIF2_f${female}_at1_t1t1') (`CIF2_f${female}_at2_t1t1') (`CIF2_f${female}_at3_t1t1') (`CIF2_f${female}_at4_t1t1') (`CIF2_f${female}_at5_t1t1') (`CIF2_f${female}_at6_t1t1') ///
						(`CIF2_f${female}_at7_t1t1') (`CIF2_f${female}_at8_t1t1') (`CIF2_f${female}_at9_t1t1') (`CIF2_f${female}_at10_t1t1')  ///
						///
                    	(`bias_CIF1_f${female}_at1_t1t2') (`bias_CIF1_f${female}_at2_t1t2') (`bias_CIF1_f${female}_at3_t1t2') (`bias_CIF1_f${female}_at4_t1t2') (`bias_CIF1_f${female}_at5_t1t2') (`bias_CIF1_f${female}_at6_t1t2') ///
						(`bias_CIF1_f${female}_at7_t1t2') (`bias_CIF1_f${female}_at8_t1t2') (`bias_CIF1_f${female}_at9_t1t2') (`bias_CIF1_f${female}_at10_t1t2')  ///
						///
					    (`bias_CIF2_f${female}_at1_t1t2') (`bias_CIF2_f${female}_at2_t1t2') (`bias_CIF2_f${female}_at3_t1t2') (`bias_CIF2_f${female}_at4_t1t2') (`bias_CIF2_f${female}_at5_t1t2') (`bias_CIF2_f${female}_at6_t1t2') ///
						(`bias_CIF2_f${female}_at7_t1t2') (`bias_CIF2_f${female}_at8_t1t2') (`bias_CIF2_f${female}_at9_t1t2') (`bias_CIF2_f${female}_at10_t1t2') ///
					    ///
						(`bias_CIF1_f${female}_at1_t1t1') (`bias_CIF1_f${female}_at2_t1t1') (`bias_CIF1_f${female}_at3_t1t1') (`bias_CIF1_f${female}_at4_t1t1') (`bias_CIF1_f${female}_at5_t1t1') (`bias_CIF1_f${female}_at6_t1t1') ///
						(`bias_CIF1_f${female}_at7_t1t1') (`bias_CIF1_f${female}_at8_t1t1') (`bias_CIF1_f${female}_at9_t1t1') (`bias_CIF1_f${female}_at10_t1t1')  ///
						///
					    (`bias_CIF2_f${female}_at1_t1t1') (`bias_CIF2_f${female}_at2_t1t1') (`bias_CIF2_f${female}_at3_t1t1') (`bias_CIF2_f${female}_at4_t1t1') (`bias_CIF2_f${female}_at5_t1t1') (`bias_CIF2_f${female}_at6_t1t1') ///
						(`bias_CIF2_f${female}_at7_t1t1') (`bias_CIF2_f${female}_at8_t1t1') (`bias_CIF2_f${female}_at9_t1t1') (`bias_CIF2_f${female}_at10_t1t1') ///
						///
						(`MCer_CIF1_f${female}_at1_t1t2') (`MCer_CIF1_f${female}_at2_t1t2') (`MCer_CIF1_f${female}_at3_t1t2') (`MCer_CIF1_f${female}_at4_t1t2') (`MCer_CIF1_f${female}_at5_t1t2') (`MCer_CIF1_f${female}_at6_t1t2') ///
						(`MCer_CIF1_f${female}_at7_t1t2') (`MCer_CIF1_f${female}_at8_t1t2') (`MCer_CIF1_f${female}_at9_t1t2') (`MCer_CIF1_f${female}_at10_t1t2')  ///
						///
					    (`MCer_CIF2_f${female}_at1_t1t2') (`MCer_CIF2_f${female}_at2_t1t2') (`MCer_CIF2_f${female}_at3_t1t2') (`MCer_CIF2_f${female}_at4_t1t2') (`MCer_CIF2_f${female}_at5_t1t2') (`MCer_CIF2_f${female}_at6_t1t2') ///
						(`MCer_CIF2_f${female}_at7_t1t2') (`MCer_CIF2_f${female}_at8_t1t2') (`MCer_CIF2_f${female}_at9_t1t2') (`MCer_CIF2_f${female}_at10_t1t2') ///
					    ///
						(`MCer_CIF1_f${female}_at1_t1t1') (`MCer_CIF1_f${female}_at2_t1t1') (`MCer_CIF1_f${female}_at3_t1t1') (`MCer_CIF1_f${female}_at4_t1t1') (`MCer_CIF1_f${female}_at5_t1t1') (`MCer_CIF1_f${female}_at6_t1t1') ///
						(`MCer_CIF1_f${female}_at7_t1t1') (`MCer_CIF1_f${female}_at8_t1t1') (`MCer_CIF1_f${female}_at9_t1t1') (`MCer_CIF1_f${female}_at10_t1t1')  ///
						///
					    (`MCer_CIF2_f${female}_at1_t1t1') (`MCer_CIF2_f${female}_at2_t1t1') (`MCer_CIF2_f${female}_at3_t1t1') (`MCer_CIF2_f${female}_at4_t1t1') (`MCer_CIF2_f${female}_at5_t1t1') (`MCer_CIF2_f${female}_at6_t1t1') ///
						(`MCer_CIF2_f${female}_at7_t1t1') (`MCer_CIF2_f${female}_at8_t1t1') (`MCer_CIF2_f${female}_at9_t1t1') (`MCer_CIF2_f${female}_at10_t1t1') ///
						///
					    ///
						(`rb_MCer_CIF1_f${female}_at1_t1t2') (`rb_MCer_CIF1_f${female}_at2_t1t2') (`rb_MCer_CIF1_f${female}_at3_t1t2') (`rb_MCer_CIF1_f${female}_at4_t1t2') (`rb_MCer_CIF1_f${female}_at5_t1t2') (`rb_MCer_CIF1_f${female}_at6_t1t2') ///
						(`rb_MCer_CIF1_f${female}_at7_t1t2') (`rb_MCer_CIF1_f${female}_at8_t1t2') (`rb_MCer_CIF1_f${female}_at9_t1t2') (`rb_MCer_CIF1_f${female}_at10_t1t2')  ///
						///
					    (`rb_MCer_CIF2_f${female}_at1_t1t2') (`rb_MCer_CIF2_f${female}_at2_t1t2') (`rb_MCer_CIF2_f${female}_at3_t1t2') (`rb_MCer_CIF2_f${female}_at4_t1t2') (`rb_MCer_CIF2_f${female}_at5_t1t2') (`rb_MCer_CIF2_f${female}_at6_t1t2') ///
						(`rb_MCer_CIF2_f${female}_at7_t1t2') (`rb_MCer_CIF2_f${female}_at8_t1t2') (`rb_MCer_CIF2_f${female}_at9_t1t2') (`rb_MCer_CIF2_f${female}_at10_t1t2') ///
					    ///
						(`rb_MCer_CIF1_f${female}_at1_t1t1') (`rb_MCer_CIF1_f${female}_at2_t1t1') (`rb_MCer_CIF1_f${female}_at3_t1t1') (`rb_MCer_CIF1_f${female}_at4_t1t1') (`rb_MCer_CIF1_f${female}_at5_t1t1') (`rb_MCer_CIF1_f${female}_at6_t1t1') ///
						(`rb_MCer_CIF1_f${female}_at7_t1t1') (`rb_MCer_CIF1_f${female}_at8_t1t1') (`rb_MCer_CIF1_f${female}_at9_t1t1') (`rb_MCer_CIF1_f${female}_at10_t1t1')  ///
						///
					    (`rb_MCer_CIF2_f${female}_at1_t1t1') (`rb_MCer_CIF2_f${female}_at2_t1t1') (`rb_MCer_CIF2_f${female}_at3_t1t1') (`rb_MCer_CIF2_f${female}_at4_t1t1') (`rb_MCer_CIF2_f${female}_at5_t1t1') (`rb_MCer_CIF2_f${female}_at6_t1t1') ///
						(`rb_MCer_CIF2_f${female}_at7_t1t1') (`rb_MCer_CIF2_f${female}_at8_t1t1') (`rb_MCer_CIF2_f${female}_at9_t1t1') (`rb_MCer_CIF2_f${female}_at10_t1t1') ///
						///
						(`rel_bias_CIF1_f${female}_at1_t1t1')  (`rel_bias_CIF1_f${female}_at2_t1t1')  (`rel_bias_CIF1_f${female}_at3_t1t1') (`rel_bias_CIF1_f${female}_at4_t1t1')  (`rel_bias_CIF1_f${female}_at5_t1t1')  (`rel_bias_CIF1_f${female}_at6_t1t1')   ///
						(`rel_bias_CIF1_f${female}_at7_t1t1')  (`rel_bias_CIF1_f${female}_at8_t1t1')  (`rel_bias_CIF1_f${female}_at9_t1t1') (`rel_bias_CIF1_f${female}_at10_t1t1') ///
						///
						(`rel_bias_CIF2_f${female}_at1_t1t1')  (`rel_bias_CIF2_f${female}_at2_t1t1')  (`rel_bias_CIF2_f${female}_at3_t1t1') (`rel_bias_CIF2_f${female}_at4_t1t1')  (`rel_bias_CIF2_f${female}_at5_t1t1')  (`rel_bias_CIF2_f${female}_at6_t1t1')   ///
						(`rel_bias_CIF2_f${female}_at7_t1t1')  (`rel_bias_CIF2_f${female}_at8_t1t1')  (`rel_bias_CIF2_f${female}_at9_t1t1') (`rel_bias_CIF2_f${female}_at10_t1t1') ///
						///
						(`rel_bias_CIF1_f${female}_at1_t1t2')  (`rel_bias_CIF1_f${female}_at2_t1t2')  (`rel_bias_CIF1_f${female}_at3_t1t2') (`rel_bias_CIF1_f${female}_at4_t1t2')  (`rel_bias_CIF1_f${female}_at5_t1t2')  (`rel_bias_CIF1_f${female}_at6_t1t2')  ///
						(`rel_bias_CIF1_f${female}_at7_t1t2')  (`rel_bias_CIF1_f${female}_at8_t1t2')  (`rel_bias_CIF1_f${female}_at9_t1t2') (`rel_bias_CIF1_f${female}_at10_t1t2')  ///
						///
						(`rel_bias_CIF2_f${female}_at1_t1t2')  (`rel_bias_CIF2_f${female}_at2_t1t2')  (`rel_bias_CIF2_f${female}_at3_t1t2') (`rel_bias_CIF2_f${female}_at4_t1t2')  (`rel_bias_CIF2_f${female}_at5_t1t2')  (`rel_bias_CIF2_f${female}_at6_t1t2')  ///
						(`rel_bias_CIF2_f${female}_at7_t1t2')  (`rel_bias_CIF2_f${female}_at8_t1t2')  (`rel_bias_CIF2_f${female}_at9_t1t2') (`rel_bias_CIF2_f${female}_at10_t1t2')  ///
						///
						(`EmpSE_CIF1_f${female}_at1_t1t1')  (`EmpSE_CIF1_f${female}_at2_t1t1')  (`EmpSE_CIF1_f${female}_at3_t1t1')  (`EmpSE_CIF1_f${female}_at4_t1t1')  (`EmpSE_CIF1_f${female}_at5_t1t1')  (`EmpSE_CIF1_f${female}_at6_t1t1') ///
						(`EmpSE_CIF1_f${female}_at7_t1t1')  (`EmpSE_CIF1_f${female}_at8_t1t1')  (`EmpSE_CIF1_f${female}_at9_t1t1') (`EmpSE_CIF1_f${female}_at10_t1t1') ///
						///
						(`EmpSE_CIF2_f${female}_at1_t1t1')  (`EmpSE_CIF2_f${female}_at2_t1t1')  (`EmpSE_CIF2_f${female}_at3_t1t1')  (`EmpSE_CIF2_f${female}_at4_t1t1')  (`EmpSE_CIF2_f${female}_at5_t1t1')  (`EmpSE_CIF2_f${female}_at6_t1t1') ///
						(`EmpSE_CIF2_f${female}_at7_t1t1')  (`EmpSE_CIF2_f${female}_at8_t1t1')  (`EmpSE_CIF2_f${female}_at9_t1t1') (`EmpSE_CIF2_f${female}_at10_t1t1') ///
						///
						(`EmpSE_CIF1_f${female}_at1_t1t2')  (`EmpSE_CIF1_f${female}_at2_t1t2')  (`EmpSE_CIF1_f${female}_at3_t1t2') (`EmpSE_CIF1_f${female}_at4_t1t2')  (`EmpSE_CIF1_f${female}_at5_t1t2')  (`EmpSE_CIF1_f${female}_at6_t1t2')   ///
						(`EmpSE_CIF1_f${female}_at7_t1t2')  (`EmpSE_CIF1_f${female}_at8_t1t2')  (`EmpSE_CIF1_f${female}_at9_t1t2')  (`EmpSE_CIF1_f${female}_at10_t1t2') ///
						///
						(`EmpSE_CIF2_f${female}_at1_t1t2')  (`EmpSE_CIF2_f${female}_at2_t1t2')  (`EmpSE_CIF2_f${female}_at3_t1t2') (`EmpSE_CIF2_f${female}_at4_t1t2')  (`EmpSE_CIF2_f${female}_at5_t1t2')  (`EmpSE_CIF2_f${female}_at6_t1t2')   ///
						(`EmpSE_CIF2_f${female}_at7_t1t2')  (`EmpSE_CIF2_f${female}_at8_t1t2')  (`EmpSE_CIF2_f${female}_at9_t1t2')  (`EmpSE_CIF2_f${female}_at10_t1t2') ///
						///
						(`Rel_prec_CIF1_f${female}_at1')  (`Rel_prec_CIF1_f${female}_at2')  (`Rel_prec_CIF1_f${female}_at3') (`Rel_prec_CIF1_f${female}_at4')  (`Rel_prec_CIF1_f${female}_at5')  (`Rel_prec_CIF1_f${female}_at6')  ///
						(`Rel_prec_CIF1_f${female}_at7')  (`Rel_prec_CIF1_f${female}_at8')  (`Rel_prec_CIF1_f${female}_at9') (`Rel_prec_CIF1_f${female}_at10') ///
						///
						(`Rel_prec_CIF2_f${female}_at1')  (`Rel_prec_CIF2_f${female}_at2')  (`Rel_prec_CIF2_f${female}_at3') (`Rel_prec_CIF2_f${female}_at4')  (`Rel_prec_CIF2_f${female}_at5')  (`Rel_prec_CIF2_f${female}_at6')  ///
						(`Rel_prec_CIF2_f${female}_at7')  (`Rel_prec_CIF2_f${female}_at8')  (`Rel_prec_CIF2_f${female}_at9') (`Rel_prec_CIF2_f${female}_at10') ///
						///
						(`Coverage_CIF1_f${female}_at1_t1t1')  (`Coverage_CIF1_f${female}_at2_t1t1')  (`Coverage_CIF1_f${female}_at3_t1t1') (`Coverage_CIF1_f${female}_at4_t1t1')  (`Coverage_CIF1_f${female}_at5_t1t1')  (`Coverage_CIF1_f${female}_at6_t1t1') ///
						(`Coverage_CIF1_f${female}_at7_t1t1')  (`Coverage_CIF1_f${female}_at8_t1t1')  (`Coverage_CIF1_f${female}_at9_t1t1') (`Coverage_CIF1_f${female}_at10_t1t1')  ///
						///
						(`Coverage_CIF2_f${female}_at1_t1t1')  (`Coverage_CIF2_f${female}_at2_t1t1')  (`Coverage_CIF2_f${female}_at3_t1t1') (`Coverage_CIF2_f${female}_at4_t1t1')  (`Coverage_CIF2_f${female}_at5_t1t1')  (`Coverage_CIF2_f${female}_at6_t1t1') ///
						(`Coverage_CIF2_f${female}_at7_t1t1')  (`Coverage_CIF2_f${female}_at8_t1t1')  (`Coverage_CIF2_f${female}_at9_t1t1') (`Coverage_CIF2_f${female}_at10_t1t1')  ///
						///
						(`Coverage_CIF1_f${female}_at1_t1t2')  (`Coverage_CIF1_f${female}_at2_t1t2')  (`Coverage_CIF1_f${female}_at3_t1t2')  (`Coverage_CIF1_f${female}_at4_t1t2')  (`Coverage_CIF1_f${female}_at5_t1t2')  (`Coverage_CIF1_f${female}_at6_t1t2') ///
						(`Coverage_CIF1_f${female}_at7_t1t2')  (`Coverage_CIF1_f${female}_at8_t1t2')  (`Coverage_CIF1_f${female}_at9_t1t2') (`Coverage_CIF1_f${female}_at10_t1t2') ///
						///
				        (`Coverage_CIF2_f${female}_at1_t1t2')  (`Coverage_CIF2_f${female}_at2_t1t2')  (`Coverage_CIF2_f${female}_at3_t1t2')  (`Coverage_CIF2_f${female}_at4_t1t2')  (`Coverage_CIF2_f${female}_at5_t1t2')  (`Coverage_CIF2_f${female}_at6_t1t2') ///
						(`Coverage_CIF2_f${female}_at7_t1t2')  (`Coverage_CIF2_f${female}_at8_t1t2')  (`Coverage_CIF2_f${female}_at9_t1t2') (`Coverage_CIF2_f${female}_at10_t1t2') ///
						///
						(`hr_gender_at1_t1t2') (`hr_gender_at2_t1t2') (`hr_gender_at3_t1t2') (`hr_gender_at4_t1t2') (`hr_gender_at5_t1t2')  ///
						(`hr_gender_at6_t1t2') (`hr_gender_at7_t1t2') (`hr_gender_at8_t1t2') (`hr_gender_at9_t1t2') (`hr_gender_at10_t1t2') ///
						///
						(`hr_gender_at1_t1t1') (`hr_gender_at2_t1t1') (`hr_gender_at3_t1t1') (`hr_gender_at4_t1t1') (`hr_gender_at5_t1t1') ///
						(`hr_gender_at6_t1t1') (`hr_gender_at7_t1t1') (`hr_gender_at8_t1t1') (`hr_gender_at9_t1t1') (`hr_gender_at10_t1t1') ///
						///
						(`bias_hr_gender_at1_t1t2') (`bias_hr_gender_at2_t1t2') (`bias_hr_gender_at3_t1t2') (`bias_hr_gender_at4_t1t2') (`bias_hr_gender_at5_t1t2')  ///
						(`bias_hr_gender_at6_t1t2') (`bias_hr_gender_at7_t1t2') (`bias_hr_gender_at8_t1t2') (`bias_hr_gender_at9_t1t2') (`bias_hr_gender_at10_t1t2') ///
						///
						(`bias_hr_gender_at1_t1t1') (`bias_hr_gender_at2_t1t1') (`bias_hr_gender_at3_t1t1') (`bias_hr_gender_at4_t1t1') (`bias_hr_gender_at5_t1t1') ///
						(`bias_hr_gender_at6_t1t1') (`bias_hr_gender_at7_t1t1') (`bias_hr_gender_at8_t1t1') (`bias_hr_gender_at9_t1t1') (`bias_hr_gender_at10_t1t1') ///
						///
						(`CIF1_diff_at1_t1t2') (`CIF1_diff_at2_t1t2') (`CIF1_diff_at3_t1t2') (`CIF1_diff_at4_t1t2') (`CIF1_diff_at5_t1t2') (`CIF1_diff_at6_t1t2') ///
						(`CIF1_diff_at7_t1t2') (`CIF1_diff_at8_t1t2') (`CIF1_diff_at9_t1t2') (`CIF1_diff_at10_t1t2')  ///
						///
			          	(`CIF1_diff_at1_t1t1') (`CIF1_diff_at2_t1t1') (`CIF1_diff_at3_t1t1') (`CIF1_diff_at4_t1t1') (`CIF1_diff_at5_t1t1') /// 
						(`CIF1_diff_at6_t1t1') (`CIF1_diff_at7_t1t1') (`CIF1_diff_at8_t1t1') (`CIF1_diff_at9_t1t1') (`CIF1_diff_at10_t1t1') ///
						///
					    (`CIF2_diff_at1_t1t2') (`CIF2_diff_at2_t1t2') (`CIF2_diff_at3_t1t2') (`CIF2_diff_at4_t1t2') (`CIF2_diff_at5_t1t2') ///
						(`CIF2_diff_at6_t1t2') (`CIF2_diff_at7_t1t2') (`CIF2_diff_at8_t1t2') (`CIF2_diff_at9_t1t2') (`CIF2_diff_at10_t1t2') ///
						///
					    (`CIF2_diff_at1_t1t1') (`CIF2_diff_at2_t1t1') (`CIF2_diff_at3_t1t1') (`CIF2_diff_at4_t1t1') (`CIF2_diff_at5_t1t1') ///
						(`CIF2_diff_at6_t1t1') (`CIF2_diff_at7_t1t1') (`CIF2_diff_at8_t1t1') (`CIF2_diff_at9_t1t1') (`CIF2_diff_at10_t1t1') ///
						///
						(`CIF1_ratio_at1_t1t2') (`CIF1_ratio_at2_t1t2') (`CIF1_ratio_at3_t1t2') (`CIF1_ratio_at4_t1t2') (`CIF1_ratio_at5_t1t2') (`CIF1_ratio_at6_t1t2') ///
						(`CIF1_ratio_at7_t1t2') (`CIF1_ratio_at8_t1t2') (`CIF1_ratio_at9_t1t2') (`CIF1_ratio_at10_t1t2')  ///
						///
			          	(`CIF1_ratio_at1_t1t1') (`CIF1_ratio_at2_t1t1') (`CIF1_ratio_at3_t1t1') (`CIF1_ratio_at4_t1t1') (`CIF1_ratio_at5_t1t1') /// 
						(`CIF1_ratio_at6_t1t1') (`CIF1_ratio_at7_t1t1') (`CIF1_ratio_at8_t1t1') (`CIF1_ratio_at9_t1t1') (`CIF1_ratio_at10_t1t1') ///
						///
					    (`CIF2_ratio_at1_t1t2') (`CIF2_ratio_at2_t1t2') (`CIF2_ratio_at3_t1t2') (`CIF2_ratio_at4_t1t2') (`CIF2_ratio_at5_t1t2') ///
						(`CIF2_ratio_at6_t1t2') (`CIF2_ratio_at7_t1t2') (`CIF2_ratio_at8_t1t2') (`CIF2_ratio_at9_t1t2') (`CIF2_ratio_at10_t1t2') ///
						///
					    (`CIF2_ratio_at1_t1t1') (`CIF2_ratio_at2_t1t1') (`CIF2_ratio_at3_t1t1') (`CIF2_ratio_at4_t1t1') (`CIF2_ratio_at5_t1t1') ///
						(`CIF2_ratio_at6_t1t1') (`CIF2_ratio_at7_t1t1') (`CIF2_ratio_at8_t1t1') (`CIF2_ratio_at9_t1t1') (`CIF2_ratio_at10_t1t1') ///
						///
						(`CIF_total_at1_t1t2') (`CIF_total_at2_t1t2') (`CIF_total_at3_t1t2') (`CIF_total_at4_t1t2') (`CIF_total_at5_t1t2') ///
						(`CIF_total_at6_t1t2') (`CIF_total_at7_t1t2') (`CIF_total_at8_t1t2') (`CIF_total_at9_t1t2') (`CIF_total_at10_t1t2')  ///
						///
					  	(`CIF_total_at1_t1t1') (`CIF_total_at2_t1t1') (`CIF_total_at3_t1t1') (`CIF_total_at4_t1t1') (`CIF_total_at5_t1t1') ///
						(`CIF_total_at6_t1t1') (`CIF_total_at7_t1t1') (`CIF_total_at8_t1t1') (`CIF_total_at9_t1t1') (`CIF_total_at10_t1t1')  ///
						///
						(`CIF1_f${female}_at1_true') (`CIF1_f${female}_at2_true') (`CIF1_f${female}_at3_true') (`CIF1_f${female}_at4_true') (`CIF1_f${female}_at5_true')  ///
						(`CIF1_f${female}_at6_true') (`CIF1_f${female}_at7_true') (`CIF1_f${female}_at8_true') (`CIF1_f${female}_at9_true') (`CIF1_f${female}_at10_true')  ///
						///
						(`CIF2_f${female}_at1_true') (`CIF2_f${female}_at2_true') (`CIF2_f${female}_at3_true') (`CIF2_f${female}_at4_true') (`CIF2_f${female}_at5_true')  ///
						(`CIF2_f${female}_at6_true') (`CIF2_f${female}_at7_true') (`CIF2_f${female}_at8_true') (`CIF2_f${female}_at9_true') (`CIF2_f${female}_at10_true')  ///
						///
						(`CIF1_diff_at1_true') 	(`CIF1_diff_at2_true') (`CIF1_diff_at3_true') (`CIF1_diff_at4_true') (`CIF1_diff_at5_true') ///
						(`CIF1_diff_at6_true')  (`CIF1_diff_at7_true') (`CIF1_diff_at8_true') (`CIF1_diff_at9_true')  (`CIF1_diff_at10_true')  ///
						///
						(`CIF2_diff_at1_true') (`CIF2_diff_at2_true') (`CIF2_diff_at3_true') (`CIF2_diff_at4_true') (`CIF2_diff_at5_true') ///
						(`CIF2_diff_at6_true') (`CIF2_diff_at7_true') (`CIF2_diff_at8_true') (`CIF2_diff_at9_true') (`CIF2_diff_at10_true') ///
						///
						(`CIF1_ratio_at1_true') (`CIF1_ratio_at2_true') (`CIF1_ratio_at3_true') (`CIF1_ratio_at4_true') (`CIF1_ratio_at5_true') ///
						(`CIF1_ratio_at6_true') (`CIF1_ratio_at7_true') (`CIF1_ratio_at8_true') (`CIF1_ratio_at9_true') (`CIF1_ratio_at10_true') ///
						///
						(`CIF2_ratio_at1_true') (`CIF2_ratio_at2_true') (`CIF2_ratio_at3_true') (`CIF2_ratio_at4_true') (`CIF2_ratio_at5_true')  ///
						(`CIF2_ratio_at6_true') (`CIF2_ratio_at7_true') (`CIF2_ratio_at8_true') (`CIF2_ratio_at9_true') (`CIF2_ratio_at10_true') ///
						///
						(`CIF_total_at1_true') (`CIF_total_at2_true') (`CIF_total_at3_true') (`CIF_total_at4_true') (`CIF_total_at5_true') ///
						(`CIF_total_at6_true') (`CIF_total_at7_true') (`CIF_total_at8_true') (`CIF_total_at9_true') (`CIF_total_at10_true') ///
						///
						(`hr_gender_at1_true') (`hr_gender_at2_true') (`hr_gender_at3_true') (`hr_gender_at4_true') (`hr_gender_at5_true') ///
						(`hr_gender_at6_true') (`hr_gender_at7_true') (`hr_gender_at8_true') (`hr_gender_at9_true') (`hr_gender_at10_true') ///
						///
						(`obs')  (`givenage') (`scenario')  (`model') (`independence') (`female') (`convergence_t1t1') (`convergence_t1t2') (`sdage') ("`prop_sex_other'")  ("`base_function'") 
						

}
postclose  comparison

}


use "$N\5b.Performance_results\compare_splines_int\scen_`scen'_compare_splines_int_f1.dta", replace

merge 1:1 _n using "$N\5b.Performance_results\compare_splines_int\scen_`scen'_compare_splines_int_f0.dta"
 
save "$N\5b.Performance_results\compare_splines_int\scen_`scen'_compare_splines_int.dta", replace

}

