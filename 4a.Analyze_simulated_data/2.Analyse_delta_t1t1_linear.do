

//Analysis using FPM for the cause specific hazards and command standsurv to derive the CIF
//point estimates and confidence intervals over scenarios and over the simulated datasets

//Approach b- Linear




clear
	
	global  times "1 2 3 4 5 6 7 8 9 10"
	
	global interval = 10
	
 
	//Analysis over the scenarios

	foreach scen in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24   { //
	
	

	// file to store results
 	postutil clear
  
  	postfile sim_t1t1_linear CIF1_f1_at1 CIF1_f1_at2 CIF1_f1_at3 CIF1_f1_at4 CIF1_f1_at5 CIF1_f1_at6 CIF1_f1_at7 CIF1_f1_at8 CIF1_f1_at9 CIF1_f1_at10 ///
					  CIF1_f1_at1_LL CIF1_f1_at1_UL CIF1_f1_at2_LL CIF1_f1_at2_UL CIF1_f1_at3_LL CIF1_f1_at3_UL CIF1_f1_at4_LL CIF1_f1_at4_UL CIF1_f1_at5_LL CIF1_f1_at5_UL ///
					  CIF1_f1_at6_LL CIF1_f1_at6_UL CIF1_f1_at7_LL CIF1_f1_at7_UL CIF1_f1_at8_LL CIF1_f1_at8_UL CIF1_f1_at9_LL CIF1_f1_at9_UL CIF1_f1_at10_LL CIF1_f1_at10_UL ///
					  CIF2_f1_at1 CIF2_f1_at2 CIF2_f1_at3 CIF2_f1_at4 CIF2_f1_at5 CIF2_f1_at6 CIF2_f1_at7 CIF2_f1_at8 CIF2_f1_at9 CIF2_f1_at10 ///
					  CIF2_f1_at1_LL CIF2_f1_at1_UL CIF2_f1_at2_LL CIF2_f1_at2_UL CIF2_f1_at3_LL CIF2_f1_at3_UL CIF2_f1_at4_LL CIF2_f1_at4_UL CIF2_f1_at5_LL CIF2_f1_at5_UL ///
					  CIF2_f1_at6_LL CIF2_f1_at6_UL CIF2_f1_at7_LL CIF2_f1_at7_UL CIF2_f1_at8_LL CIF2_f1_at8_UL CIF2_f1_at9_LL CIF2_f1_at9_UL CIF2_f1_at10_LL CIF2_f1_at10_UL ///
					  ///
					  CIF1_f0_at1 CIF1_f0_at2 CIF1_f0_at3 CIF1_f0_at4 CIF1_f0_at5 CIF1_f0_at6 CIF1_f0_at7 CIF1_f0_at8 CIF1_f0_at9 CIF1_f0_at10 ///
					  CIF1_f0_at1_LL CIF1_f0_at1_UL CIF1_f0_at2_LL CIF1_f0_at2_UL CIF1_f0_at3_LL CIF1_f0_at3_UL CIF1_f0_at4_LL CIF1_f0_at4_UL CIF1_f0_at5_LL CIF1_f0_at5_UL ///
					  CIF1_f0_at6_LL CIF1_f0_at6_UL CIF1_f0_at7_LL CIF1_f0_at7_UL CIF1_f0_at8_LL CIF1_f0_at8_UL CIF1_f0_at9_LL CIF1_f0_at9_UL CIF1_f0_at10_LL CIF1_f0_at10_UL ///
					  CIF2_f0_at1 CIF2_f0_at2 CIF2_f0_at3 CIF2_f0_at4 CIF2_f0_at5 CIF2_f0_at6 CIF2_f0_at7 CIF2_f0_at8 CIF2_f0_at9 CIF2_f0_at10 ///
					  CIF2_f0_at1_LL CIF2_f0_at1_UL CIF2_f0_at2_LL CIF2_f0_at2_UL CIF2_f0_at3_LL CIF2_f0_at3_UL CIF2_f0_at4_LL CIF2_f0_at4_UL CIF2_f0_at5_LL CIF2_f0_at5_UL ///
					  CIF2_f0_at6_LL CIF2_f0_at6_UL CIF2_f0_at7_LL CIF2_f0_at7_UL CIF2_f0_at8_LL CIF2_f0_at8_UL CIF2_f0_at9_LL CIF2_f0_at9_UL CIF2_f0_at10_LL CIF2_f0_at10_UL ///
					  ///
					  CIF1_diff_at1 CIF1_diff_at2 CIF1_diff_at3 CIF1_diff_at4 CIF1_diff_at5 CIF1_diff_at6 CIF1_diff_at7 CIF1_diff_at8 CIF1_diff_at9 CIF1_diff_at10 ///
					  CIF2_diff_at1 CIF2_diff_at2 CIF2_diff_at3 CIF2_diff_at4 CIF2_diff_at5 CIF2_diff_at6 CIF2_diff_at7 CIF2_diff_at8 CIF2_diff_at9 CIF2_diff_at10 ///
					  ///
					  CIF1_ratio_at1 CIF1_ratio_at2 CIF1_ratio_at3 CIF1_ratio_at4 CIF1_ratio_at5 CIF1_ratio_at6 CIF1_ratio_at7 CIF1_ratio_at8 CIF1_ratio_at9 CIF1_ratio_at10 ///
					  CIF2_ratio_at1 CIF2_ratio_at2 CIF2_ratio_at3 CIF2_ratio_at4 CIF2_ratio_at5 CIF2_ratio_at6 CIF2_ratio_at7 CIF2_ratio_at8 CIF2_ratio_at9 CIF2_ratio_at10 ///
					  ///
					  CIF_total_at1 CIF_total_at2 CIF_total_at3 CIF_total_at4 CIF_total_at5 CIF_total_at6 CIF_total_at7 CIF_total_at8 CIF_total_at9 CIF_total_at10  ///
					  ///
					  hr_gender_at1 hr_gender_at2 hr_gender_at3 hr_gender_at4 hr_gender_at5  ///
					  hr_gender_at6 hr_gender_at7 hr_gender_at8 hr_gender_at9 hr_gender_at10 givenage ///
					  canc_model_converged other_model_converged ///
			          using ///		 
					 "$N\4b.Analyze_results\Approach_b_linear/scenario`scen'_linear.dta", replace

					 include "$N\1b.Scenarios_folder/Scenario`scen'_new.do"
	
//Analysis for each dataset

	forvalues i = 1/1000	 {
	
  
	
	capture constraint drop _all
	
	// load data
	clear
	use "$N\3b.Simulated_data_folder\Scenario`scen'\Sim_data`i'.dta", clear
	
	//Make the extremely low or high values of age at diagnosis equal to the 2nd and 98th percentile respectively

	gen id=_n
	_pctile agediag, per(2)
	local age_cutoff_low `r(r1)'
	di "`age_cutoff_low'"
	gen ageadj = cond(agediag<`age_cutoff_low',`age_cutoff_low',agediag)
	_pctile agediag, per(98)
	local age_cutoff_high `r(r1)'
	replace ageadj = cond(ageadj>`age_cutoff_high',`age_cutoff_high',ageadj)

	//Create the spline variables of aage at diagnosis

	rcsgen ageadj, df(4) gen(agercs) 
	global age_knots `r(knots)'
	
	rcsgen, scalar(50) knots(${age_knots}) gen(c50_)
	rcsgen, scalar(60) knots(${age_knots}) gen(c60_)
	rcsgen, scalar(70) knots(${age_knots}) gen(c70_)
	rcsgen, scalar(80) knots(${age_knots}) gen(c80_)
	rcsgen, scalar(90) knots(${age_knots}) gen(c90_)
	
	drop agediagc
	drop agediag
	rename ageadj agediag
	
	sum agediag
	gen agediagc=agediag- `r(mean)' 


    gen attained_age= agediag+time_s_diag
	
	//Analysis for cause I- Cancer mortality- Splines for age at diagnosis

	
			stset time_s_diag, failure(cancerd ==1)  id(id) exit(time $interval) 
	
			capture stpm2 agercs1 agercs2 agercs3 agercs4 female , df(5)  scale(hazard)  eform  failconvlininit  iter(100)	
			
			local canc_model_converged = cond(_rc ,0,1)
			
			di "`canc_model_converged'"

			capture est store cancer
	
	//Analysis for cause II- Other cause mortality- Age at diagnosis not included- Attained age timescale

			stset time_s_diag , failure(otherd ==1)   id(id) exit(time $interval) 
			
	//Analysis allowing only proportional hazards for gender
			
			if "$prop_sex_other"=="PH_sex_other" {
	
				capture stpm2 agediag female,   df(5)  scale(hazard)   eform 
	
				local other_model_converged = cond(_rc ,0,1)
				di "`other_model_converged'"
	
				capture est store other
			
			}
			
	//Analysis allowing non proportional hazards for gender
	
			else if "$prop_sex_other"=="Non_PH_sex_other" {
			
				capture stpm2 agediag female,  tvc(female) dftvc(2)  df(5)  scale(hazard)   eform  failconvlininit  iter(100)
	
				local other_model_converged = cond(_rc ,0,1)
				di "`other_model_converged'"
	
				capture est store other
			
			}
			
			

	
	//should be replaced with time since diag variable but duplicate knots issue

	range tt1 1 10 10
	//replace tt1=1 if tt1==0
	
	/**Attention both the model for cause I and cause II should converge in order to use the standsurv and */
	/*and derive the CIFs. If one of them does not converge then we are just passing null values at the locals*/
	/*for the specific loop for the simulated dataset where convergence was not achieved*/ 	if `canc_model_converged'==1 & `other_model_converged'==1 {
		
			
	forvalues ageloop= 60(10)90 {
	
	global age=`ageloop'
	
	local c50_1_loop `=c50_1'
	local c50_2_loop `=c50_2'
	local c50_3_loop `=c50_3'
	local c50_4_loop `=c50_4'


	local c60_1_loop `=c60_1'
	local c60_2_loop `=c60_2'
	local c60_3_loop `=c60_3'
	local c60_4_loop `=c60_4'

	
	local c70_1_loop `=c70_1'
	local c70_2_loop `=c70_2'
	local c70_3_loop `=c70_3'
	local c70_4_loop `=c70_4'

	
	local c80_1_loop `=c80_1'
	local c80_2_loop `=c80_2'
	local c80_3_loop `=c80_3'
	local c80_4_loop `=c80_4'
	
	
	local c90_1_loop `=c90_1'
	local c90_2_loop `=c90_2'
	local c90_3_loop `=c90_3'
	local c90_4_loop `=c90_4'
		
	gen tt2_$age = tt1+ $age 

	// conditional estimates
	//  age 65 female 0 vs 1
	capture drop F_*
	
	gen age$age =$age in 1
	
	//Use of the standsurv to derive the CIF estimates and the confidence intervals


	standsurv if _n==1         /// restrict to one row as make conditional, not marginal estimates
	, timevar(tt1)          /// times to predict at 
	cif                    /// predict cif
	crmodels(cancer other) /// stored cause-specific models
	nodes(50)              /// number of nodes
	ci                     /// predict cis
	verbose                /// display extra information
	atvar(F_t1t2_${age}_f1 F_t1t2_${age}_f0 ) /// new variables 
	contrast(difference)   /// calculate contrast
	contrastvar(F_diff)    /// name of contrast variable
	toffset(. .)         /// offset for  */ ///
	at1(female 1 agediag ${age} agercs1 `c${age}_1_loop' agercs2 `c${age}_2_loop' agercs3 `c${age}_3_loop' agercs4 `c${age}_4_loop'  ) /// 
	at2(female 0 agediag ${age} agercs1 `c${age}_1_loop' agercs2 `c${age}_2_loop' agercs3 `c${age}_3_loop' agercs4 `c${age}_4_loop') /// 
	ode
	

	/*keep the estimates only for the list of time points specified by global macro times*/
	
	local j=1 
	foreach x in $times  {
	
	local CIF1_f1_at`x' = F_t1t2_${age}_f1_cancer[`j']
	local CIF2_f1_at`x' = F_t1t2_${age}_f1_other[`j']
	local CIF1_f1_at`x'_UL = F_t1t2_${age}_f1_cancer_uci[`j']
	local CIF1_f1_at`x'_LL = F_t1t2_${age}_f1_cancer_lci[`j']
	local CIF2_f1_at`x'_UL = F_t1t2_${age}_f1_other_uci[`j']
	local CIF2_f1_at`x'_LL = F_t1t2_${age}_f1_other_lci[`j']
	
	local CIF1_f0_at`x' = F_t1t2_${age}_f0_cancer[`j']
	local CIF2_f0_at`x' = F_t1t2_${age}_f0_other[`j']
	local CIF1_f0_at`x'_UL = F_t1t2_${age}_f0_cancer_uci[`j']
	local CIF1_f0_at`x'_LL = F_t1t2_${age}_f0_cancer_lci[`j']
	local CIF2_f0_at`x'_UL = F_t1t2_${age}_f0_other_uci[`j']
	local CIF2_f0_at`x'_LL = F_t1t2_${age}_f0_other_lci[`j']
	
	local CIF1_diff_at`x'=  F_t1t2_${age}_f1_cancer[`j']- F_t1t2_${age}_f0_cancer[`j']
	local CIF2_diff_at`x'=  F_t1t2_${age}_f1_other[`j']-F_t1t2_${age}_f0_other[`j']	
		
	local CIF1_ratio_at`x'=  F_t1t2_${age}_f1_cancer[`j']/ F_t1t2_${age}_f0_cancer[`j']
	local CIF2_ratio_at`x'= F_t1t2_${age}_f1_other[`j']/F_t1t2_${age}_f0_other[`j']
	
	local CIF_total_at`x' =(F_t1t2_${age}_f1_cancer[`j']+F_t1t2_${age}_f1_other[`j'])/(F_t1t2_${age}_f0_cancer[`j']+F_t1t2_${age}_f0_other[`j'])
	

	
local j=`j'+1

}
	est restore other
	

    predict hr_gender${age},   hrnum(female 1 agediag ${age}) ///
					         hrdenom(female 0 agediag ${age} ) timevar(tt1) 
	local p=1 
	foreach x in $times {
	local hr_gender_at`x'= hr_gender${age}[`p']
	local p=`p'+1
	}


local givenage = `ageloop'





//use the point estimates of CIFs

// store results

		post sim_t1t1_linear (`CIF1_f1_at1') (`CIF1_f1_at2') (`CIF1_f1_at3') (`CIF1_f1_at4') (`CIF1_f1_at5') (`CIF1_f1_at6') (`CIF1_f1_at7') (`CIF1_f1_at8') (`CIF1_f1_at9') (`CIF1_f1_at10') ///
					  (`CIF1_f1_at1_LL') (`CIF1_f1_at1_UL') (`CIF1_f1_at2_LL') (`CIF1_f1_at2_UL') (`CIF1_f1_at3_LL') (`CIF1_f1_at3_UL') (`CIF1_f1_at4_LL') (`CIF1_f1_at4_UL') (`CIF1_f1_at5_LL') (`CIF1_f1_at5_UL') ///
					  (`CIF1_f1_at6_LL') (`CIF1_f1_at6_UL') (`CIF1_f1_at7_LL') (`CIF1_f1_at7_UL') (`CIF1_f1_at8_LL') (`CIF1_f1_at8_UL') (`CIF1_f1_at9_LL') (`CIF1_f1_at9_UL') (`CIF1_f1_at10_LL') (`CIF1_f1_at10_UL') ///
					  (`CIF2_f1_at1') (`CIF2_f1_at2') (`CIF2_f1_at3') (`CIF2_f1_at4') (`CIF2_f1_at5') (`CIF2_f1_at6') (`CIF2_f1_at7') (`CIF2_f1_at8') (`CIF2_f1_at9') (`CIF2_f1_at10') ///
					  (`CIF2_f1_at1_LL') (`CIF2_f1_at1_UL') (`CIF2_f1_at2_LL') (`CIF2_f1_at2_UL') (`CIF2_f1_at3_LL') (`CIF2_f1_at3_UL') (`CIF2_f1_at4_LL') (`CIF2_f1_at4_UL') (`CIF2_f1_at5_LL') (`CIF2_f1_at5_UL') ///
					  (`CIF2_f1_at6_LL') (`CIF2_f1_at6_UL') (`CIF2_f1_at7_LL') (`CIF2_f1_at7_UL') (`CIF2_f1_at8_LL') (`CIF2_f1_at8_UL') (`CIF2_f1_at9_LL') (`CIF2_f1_at9_UL') (`CIF2_f1_at10_LL') (`CIF2_f1_at10_UL') ///
					  ///
					  (`CIF1_f0_at1') (`CIF1_f0_at2') (`CIF1_f0_at3') (`CIF1_f0_at4') (`CIF1_f0_at5') (`CIF1_f0_at6') (`CIF1_f0_at7') (`CIF1_f0_at8') (`CIF1_f0_at9') (`CIF1_f0_at10') ///
					  (`CIF1_f0_at1_LL') (`CIF1_f0_at1_UL') (`CIF1_f0_at2_LL') (`CIF1_f0_at2_UL') (`CIF1_f0_at3_LL') (`CIF1_f0_at3_UL') (`CIF1_f0_at4_LL') (`CIF1_f0_at4_UL') (`CIF1_f0_at5_LL') (`CIF1_f0_at5_UL') ///
					  (`CIF1_f0_at6_LL') (`CIF1_f0_at6_UL') (`CIF1_f0_at7_LL') (`CIF1_f0_at7_UL') (`CIF1_f0_at8_LL') (`CIF1_f0_at8_UL') (`CIF1_f0_at9_LL') (`CIF1_f0_at9_UL') (`CIF1_f0_at10_LL') (`CIF1_f0_at10_UL') ///
					  (`CIF2_f0_at1') (`CIF2_f0_at2') (`CIF2_f0_at3') (`CIF2_f0_at4') (`CIF2_f0_at5') (`CIF2_f0_at6') (`CIF2_f0_at7') (`CIF2_f0_at8') (`CIF2_f0_at9') (`CIF2_f0_at10') ///
					  (`CIF2_f0_at1_LL') (`CIF2_f0_at1_UL') (`CIF2_f0_at2_LL') (`CIF2_f0_at2_UL') (`CIF2_f0_at3_LL') (`CIF2_f0_at3_UL') (`CIF2_f0_at4_LL') (`CIF2_f0_at4_UL') (`CIF2_f0_at5_LL') (`CIF2_f0_at5_UL') ///
					  (`CIF2_f0_at6_LL') (`CIF2_f0_at6_UL') (`CIF2_f0_at7_LL') (`CIF2_f0_at7_UL') (`CIF2_f0_at8_LL') (`CIF2_f0_at8_UL') (`CIF2_f0_at9_LL') (`CIF2_f0_at9_UL') (`CIF2_f0_at10_LL') (`CIF2_f0_at10_UL') ///
					  ///
					  (`CIF1_diff_at1') (`CIF1_diff_at2') (`CIF1_diff_at3') (`CIF1_diff_at4') (`CIF1_diff_at5') (`CIF1_diff_at6') (`CIF1_diff_at7') (`CIF1_diff_at8') (`CIF1_diff_at9') (`CIF1_diff_at10') ///
					  (`CIF2_diff_at1') (`CIF2_diff_at2') (`CIF2_diff_at3') (`CIF2_diff_at4') (`CIF2_diff_at5') (`CIF2_diff_at6') (`CIF2_diff_at7') (`CIF2_diff_at8') (`CIF2_diff_at9') (`CIF2_diff_at10') ///
					  ///
					  (`CIF1_ratio_at1') (`CIF1_ratio_at2') (`CIF1_ratio_at3') (`CIF1_ratio_at4') (`CIF1_ratio_at5') (`CIF1_ratio_at6') (`CIF1_ratio_at7') (`CIF1_ratio_at8') (`CIF1_ratio_at9') (`CIF1_ratio_at10') ///
					  (`CIF2_ratio_at1') (`CIF2_ratio_at2') (`CIF2_ratio_at3') (`CIF2_ratio_at4') (`CIF2_ratio_at5') (`CIF2_ratio_at6') (`CIF2_ratio_at7') (`CIF2_ratio_at8') (`CIF2_ratio_at9') (`CIF2_ratio_at10') ///
					  ///
					  (`CIF_total_at1') (`CIF_total_at2') (`CIF_total_at3') (`CIF_total_at4') (`CIF_total_at5') (`CIF_total_at6') (`CIF_total_at7') (`CIF_total_at8') (`CIF_total_at9') (`CIF_total_at10') /// 
					  ///
					  (`hr_gender_at1') (`hr_gender_at2') (`hr_gender_at3') (`hr_gender_at4') (`hr_gender_at5') ///
					  (`hr_gender_at6') (`hr_gender_at7') (`hr_gender_at8') (`hr_gender_at9') (`hr_gender_at10') (`givenage') ///
					  (`canc_model_converged') (`other_model_converged')
				
		
	}
	}
	

		else if `canc_model_converged'==0 &  `other_model_converged'==0 {
	
	forvalues ageloop= 60(10)90 {
	
	
	local j= 1
	foreach x in $times {
	
	
	local CIF1_f1_at`x'    = .
	local CIF2_f1_at`x'    = .
	local CIF1_f1_at`x'_UL = .
	local CIF1_f1_at`x'_LL = .
	local CIF2_f1_at`x'_UL = .
	local CIF2_f1_at`x'_LL = .
	
	local CIF1_f0_at`x'    = .
	local CIF2_f0_at`x'    = .
	local CIF1_f0_at`x'_UL = .
	local CIF1_f0_at`x'_LL = .
	local CIF2_f0_at`x'_UL = .
	local CIF2_f0_at`x'_LL = .
	
	local CIF1_diff_at`x'= .
	local CIF2_diff_at`x'= . 
	local CIF1_ratio_at`x'=.
	local CIF2_ratio_at`x'=.
		
	local CIF_total_at`x' =.
		
	local j= `j'+1
}


	local p=1 
	foreach x in $times {
	local hr_gender_at`x'=.
	local p=`p'+1
	}


local givenage = `ageloop'

		post sim_t1t1_linear (`CIF1_f1_at1') (`CIF1_f1_at2') (`CIF1_f1_at3') (`CIF1_f1_at4') (`CIF1_f1_at5') (`CIF1_f1_at6') (`CIF1_f1_at7') (`CIF1_f1_at8') (`CIF1_f1_at9') (`CIF1_f1_at10') ///
					  (`CIF1_f1_at1_LL') (`CIF1_f1_at1_UL') (`CIF1_f1_at2_LL') (`CIF1_f1_at2_UL') (`CIF1_f1_at3_LL') (`CIF1_f1_at3_UL') (`CIF1_f1_at4_LL') (`CIF1_f1_at4_UL') (`CIF1_f1_at5_LL') (`CIF1_f1_at5_UL') ///
					  (`CIF1_f1_at6_LL') (`CIF1_f1_at6_UL') (`CIF1_f1_at7_LL') (`CIF1_f1_at7_UL') (`CIF1_f1_at8_LL') (`CIF1_f1_at8_UL') (`CIF1_f1_at9_LL') (`CIF1_f1_at9_UL') (`CIF1_f1_at10_LL') (`CIF1_f1_at10_UL') ///
					  (`CIF2_f1_at1') (`CIF2_f1_at2') (`CIF2_f1_at3') (`CIF2_f1_at4') (`CIF2_f1_at5') (`CIF2_f1_at6') (`CIF2_f1_at7') (`CIF2_f1_at8') (`CIF2_f1_at9') (`CIF2_f1_at10') ///
					  (`CIF2_f1_at1_LL') (`CIF2_f1_at1_UL') (`CIF2_f1_at2_LL') (`CIF2_f1_at2_UL') (`CIF2_f1_at3_LL') (`CIF2_f1_at3_UL') (`CIF2_f1_at4_LL') (`CIF2_f1_at4_UL') (`CIF2_f1_at5_LL') (`CIF2_f1_at5_UL') ///
					  (`CIF2_f1_at6_LL') (`CIF2_f1_at6_UL') (`CIF2_f1_at7_LL') (`CIF2_f1_at7_UL') (`CIF2_f1_at8_LL') (`CIF2_f1_at8_UL') (`CIF2_f1_at9_LL') (`CIF2_f1_at9_UL') (`CIF2_f1_at10_LL') (`CIF2_f1_at10_UL') ///
					  ///
					  (`CIF1_f0_at1') (`CIF1_f0_at2') (`CIF1_f0_at3') (`CIF1_f0_at4') (`CIF1_f0_at5') (`CIF1_f0_at6') (`CIF1_f0_at7') (`CIF1_f0_at8') (`CIF1_f0_at9') (`CIF1_f0_at10') ///
					  (`CIF1_f0_at1_LL') (`CIF1_f0_at1_UL') (`CIF1_f0_at2_LL') (`CIF1_f0_at2_UL') (`CIF1_f0_at3_LL') (`CIF1_f0_at3_UL') (`CIF1_f0_at4_LL') (`CIF1_f0_at4_UL') (`CIF1_f0_at5_LL') (`CIF1_f0_at5_UL') ///
					  (`CIF1_f0_at6_LL') (`CIF1_f0_at6_UL') (`CIF1_f0_at7_LL') (`CIF1_f0_at7_UL') (`CIF1_f0_at8_LL') (`CIF1_f0_at8_UL') (`CIF1_f0_at9_LL') (`CIF1_f0_at9_UL') (`CIF1_f0_at10_LL') (`CIF1_f0_at10_UL') ///
					  (`CIF2_f0_at1') (`CIF2_f0_at2') (`CIF2_f0_at3') (`CIF2_f0_at4') (`CIF2_f0_at5') (`CIF2_f0_at6') (`CIF2_f0_at7') (`CIF2_f0_at8') (`CIF2_f0_at9') (`CIF2_f0_at10') ///
					  (`CIF2_f0_at1_LL') (`CIF2_f0_at1_UL') (`CIF2_f0_at2_LL') (`CIF2_f0_at2_UL') (`CIF2_f0_at3_LL') (`CIF2_f0_at3_UL') (`CIF2_f0_at4_LL') (`CIF2_f0_at4_UL') (`CIF2_f0_at5_LL') (`CIF2_f0_at5_UL') ///
					  (`CIF2_f0_at6_LL') (`CIF2_f0_at6_UL') (`CIF2_f0_at7_LL') (`CIF2_f0_at7_UL') (`CIF2_f0_at8_LL') (`CIF2_f0_at8_UL') (`CIF2_f0_at9_LL') (`CIF2_f0_at9_UL') (`CIF2_f0_at10_LL') (`CIF2_f0_at10_UL') ///
					  ///
					  (`CIF1_diff_at1') (`CIF1_diff_at2') (`CIF1_diff_at3') (`CIF1_diff_at4') (`CIF1_diff_at5') (`CIF1_diff_at6') (`CIF1_diff_at7') (`CIF1_diff_at8') (`CIF1_diff_at9') (`CIF1_diff_at10') ///
					  (`CIF2_diff_at1') (`CIF2_diff_at2') (`CIF2_diff_at3') (`CIF2_diff_at4') (`CIF2_diff_at5') (`CIF2_diff_at6') (`CIF2_diff_at7') (`CIF2_diff_at8') (`CIF2_diff_at9') (`CIF2_diff_at10') ///
					  ///
					  (`CIF1_ratio_at1') (`CIF1_ratio_at2') (`CIF1_ratio_at3') (`CIF1_ratio_at4') (`CIF1_ratio_at5') (`CIF1_ratio_at6') (`CIF1_ratio_at7') (`CIF1_ratio_at8') (`CIF1_ratio_at9') (`CIF1_ratio_at10') ///
					  (`CIF2_ratio_at1') (`CIF2_ratio_at2') (`CIF2_ratio_at3') (`CIF2_ratio_at4') (`CIF2_ratio_at5') (`CIF2_ratio_at6') (`CIF2_ratio_at7') (`CIF2_ratio_at8') (`CIF2_ratio_at9') (`CIF2_ratio_at10') ///
					  ///
					  (`CIF_total_at1') (`CIF_total_at2') (`CIF_total_at3') (`CIF_total_at4') (`CIF_total_at5') (`CIF_total_at6') (`CIF_total_at7') (`CIF_total_at8') (`CIF_total_at9') (`CIF_total_at10') /// 
					  ///
					  (`hr_gender_at1') (`hr_gender_at2') (`hr_gender_at3') (`hr_gender_at4') (`hr_gender_at5') ///
					  (`hr_gender_at6') (`hr_gender_at7') (`hr_gender_at8') (`hr_gender_at9') (`hr_gender_at10') (`givenage')  ///
					  (`canc_model_converged') (`other_model_converged')
					  
		}	
	}
	
	
			else if `canc_model_converged'==1 & `other_model_converged'==0 {
	
	forvalues ageloop= 60(10)90 {
	
	
	local j= 1
	foreach x in $times {
	
	
	local CIF1_f1_at`x'    = .
	local CIF2_f1_at`x'    = .
	local CIF1_f1_at`x'_UL = .
	local CIF1_f1_at`x'_LL = .
	local CIF2_f1_at`x'_UL = .
	local CIF2_f1_at`x'_LL = .
	
	local CIF1_f0_at`x'    = .
	local CIF2_f0_at`x'    = .
	local CIF1_f0_at`x'_UL = .
	local CIF1_f0_at`x'_LL = .
	local CIF2_f0_at`x'_UL = .
	local CIF2_f0_at`x'_LL = .
	
	local CIF1_diff_at`x'= .
	local CIF2_diff_at`x'= . 
	local CIF1_ratio_at`x'=.
	local CIF2_ratio_at`x'=.
		
	local CIF_total_at`x' =.
	
	local j= `j'+1
}



	local p=1 
	foreach x in $times {
	local hr_gender_at`x'=.
	local p=`p'+1
	}


local givenage = `ageloop'

		post sim_t1t1_linear (`CIF1_f1_at1') (`CIF1_f1_at2') (`CIF1_f1_at3') (`CIF1_f1_at4') (`CIF1_f1_at5') (`CIF1_f1_at6') (`CIF1_f1_at7') (`CIF1_f1_at8') (`CIF1_f1_at9') (`CIF1_f1_at10') ///
					  (`CIF1_f1_at1_LL') (`CIF1_f1_at1_UL') (`CIF1_f1_at2_LL') (`CIF1_f1_at2_UL') (`CIF1_f1_at3_LL') (`CIF1_f1_at3_UL') (`CIF1_f1_at4_LL') (`CIF1_f1_at4_UL') (`CIF1_f1_at5_LL') (`CIF1_f1_at5_UL') ///
					  (`CIF1_f1_at6_LL') (`CIF1_f1_at6_UL') (`CIF1_f1_at7_LL') (`CIF1_f1_at7_UL') (`CIF1_f1_at8_LL') (`CIF1_f1_at8_UL') (`CIF1_f1_at9_LL') (`CIF1_f1_at9_UL') (`CIF1_f1_at10_LL') (`CIF1_f1_at10_UL') ///
					  (`CIF2_f1_at1') (`CIF2_f1_at2') (`CIF2_f1_at3') (`CIF2_f1_at4') (`CIF2_f1_at5') (`CIF2_f1_at6') (`CIF2_f1_at7') (`CIF2_f1_at8') (`CIF2_f1_at9') (`CIF2_f1_at10') ///
					  (`CIF2_f1_at1_LL') (`CIF2_f1_at1_UL') (`CIF2_f1_at2_LL') (`CIF2_f1_at2_UL') (`CIF2_f1_at3_LL') (`CIF2_f1_at3_UL') (`CIF2_f1_at4_LL') (`CIF2_f1_at4_UL') (`CIF2_f1_at5_LL') (`CIF2_f1_at5_UL') ///
					  (`CIF2_f1_at6_LL') (`CIF2_f1_at6_UL') (`CIF2_f1_at7_LL') (`CIF2_f1_at7_UL') (`CIF2_f1_at8_LL') (`CIF2_f1_at8_UL') (`CIF2_f1_at9_LL') (`CIF2_f1_at9_UL') (`CIF2_f1_at10_LL') (`CIF2_f1_at10_UL') ///
					  ///
					  (`CIF1_f0_at1') (`CIF1_f0_at2') (`CIF1_f0_at3') (`CIF1_f0_at4') (`CIF1_f0_at5') (`CIF1_f0_at6') (`CIF1_f0_at7') (`CIF1_f0_at8') (`CIF1_f0_at9') (`CIF1_f0_at10') ///
					  (`CIF1_f0_at1_LL') (`CIF1_f0_at1_UL') (`CIF1_f0_at2_LL') (`CIF1_f0_at2_UL') (`CIF1_f0_at3_LL') (`CIF1_f0_at3_UL') (`CIF1_f0_at4_LL') (`CIF1_f0_at4_UL') (`CIF1_f0_at5_LL') (`CIF1_f0_at5_UL') ///
					  (`CIF1_f0_at6_LL') (`CIF1_f0_at6_UL') (`CIF1_f0_at7_LL') (`CIF1_f0_at7_UL') (`CIF1_f0_at8_LL') (`CIF1_f0_at8_UL') (`CIF1_f0_at9_LL') (`CIF1_f0_at9_UL') (`CIF1_f0_at10_LL') (`CIF1_f0_at10_UL') ///
					  (`CIF2_f0_at1') (`CIF2_f0_at2') (`CIF2_f0_at3') (`CIF2_f0_at4') (`CIF2_f0_at5') (`CIF2_f0_at6') (`CIF2_f0_at7') (`CIF2_f0_at8') (`CIF2_f0_at9') (`CIF2_f0_at10') ///
					  (`CIF2_f0_at1_LL') (`CIF2_f0_at1_UL') (`CIF2_f0_at2_LL') (`CIF2_f0_at2_UL') (`CIF2_f0_at3_LL') (`CIF2_f0_at3_UL') (`CIF2_f0_at4_LL') (`CIF2_f0_at4_UL') (`CIF2_f0_at5_LL') (`CIF2_f0_at5_UL') ///
					  (`CIF2_f0_at6_LL') (`CIF2_f0_at6_UL') (`CIF2_f0_at7_LL') (`CIF2_f0_at7_UL') (`CIF2_f0_at8_LL') (`CIF2_f0_at8_UL') (`CIF2_f0_at9_LL') (`CIF2_f0_at9_UL') (`CIF2_f0_at10_LL') (`CIF2_f0_at10_UL') ///
					  ///
					  (`CIF1_diff_at1') (`CIF1_diff_at2') (`CIF1_diff_at3') (`CIF1_diff_at4') (`CIF1_diff_at5') (`CIF1_diff_at6') (`CIF1_diff_at7') (`CIF1_diff_at8') (`CIF1_diff_at9') (`CIF1_diff_at10') ///
					  (`CIF2_diff_at1') (`CIF2_diff_at2') (`CIF2_diff_at3') (`CIF2_diff_at4') (`CIF2_diff_at5') (`CIF2_diff_at6') (`CIF2_diff_at7') (`CIF2_diff_at8') (`CIF2_diff_at9') (`CIF2_diff_at10') ///
					  ///
					  (`CIF1_ratio_at1') (`CIF1_ratio_at2') (`CIF1_ratio_at3') (`CIF1_ratio_at4') (`CIF1_ratio_at5') (`CIF1_ratio_at6') (`CIF1_ratio_at7') (`CIF1_ratio_at8') (`CIF1_ratio_at9') (`CIF1_ratio_at10') ///
					  (`CIF2_ratio_at1') (`CIF2_ratio_at2') (`CIF2_ratio_at3') (`CIF2_ratio_at4') (`CIF2_ratio_at5') (`CIF2_ratio_at6') (`CIF2_ratio_at7') (`CIF2_ratio_at8') (`CIF2_ratio_at9') (`CIF2_ratio_at10') ///
					  ///
					  (`CIF_total_at1') (`CIF_total_at2') (`CIF_total_at3') (`CIF_total_at4') (`CIF_total_at5') (`CIF_total_at6') (`CIF_total_at7') (`CIF_total_at8') (`CIF_total_at9') (`CIF_total_at10') /// 
					  ///
					  (`hr_gender_at1') (`hr_gender_at2') (`hr_gender_at3') (`hr_gender_at4') (`hr_gender_at5') ///
					  (`hr_gender_at6') (`hr_gender_at7') (`hr_gender_at8') (`hr_gender_at9') (`hr_gender_at10') (`givenage')  ///
					  (`canc_model_converged') (`other_model_converged')
	
	}	
	}
	
	
	
	else if `canc_model_converged'==0 |  `other_model_converged'==1 {
	
	
		
			
	forvalues ageloop= 60(10)90 {
	
	global age=`ageloop'
	
	local c50_1_loop `=c50_1'
	local c50_2_loop `=c50_2'
	local c50_3_loop `=c50_3'
	local c50_4_loop `=c50_4'


	local c60_1_loop `=c60_1'
	local c60_2_loop `=c60_2'
	local c60_3_loop `=c60_3'
	local c60_4_loop `=c60_4'

	
	local c70_1_loop `=c70_1'
	local c70_2_loop `=c70_2'
	local c70_3_loop `=c70_3'
	local c70_4_loop `=c70_4'

	
	local c80_1_loop `=c80_1'
	local c80_2_loop `=c80_2'
	local c80_3_loop `=c80_3'
	local c80_4_loop `=c80_4'
	
	
	local c90_1_loop `=c90_1'
	local c90_2_loop `=c90_2'
	local c90_3_loop `=c90_3'
	local c90_4_loop `=c90_4'
	
	gen tt2_$age = tt1+ $age 

	
	local j= 1
	foreach x in $times {
	
	
	local CIF1_f1_at`x'    = .
	local CIF2_f1_at`x'    = .
	local CIF1_f1_at`x'_UL = .
	local CIF1_f1_at`x'_LL = .
	local CIF2_f1_at`x'_UL = .
	local CIF2_f1_at`x'_LL = .
	
	local CIF1_f0_at`x'    = .
	local CIF2_f0_at`x'    = .
	local CIF1_f0_at`x'_UL = .
	local CIF1_f0_at`x'_LL = .
	local CIF2_f0_at`x'_UL = .
	local CIF2_f0_at`x'_LL = .
	
	local CIF1_diff_at`x'= .
	local CIF2_diff_at`x'= . 
	local CIF1_ratio_at`x'=.
	local CIF2_ratio_at`x'=.
		
	local CIF_total_at`x' =.
	
	local j= `j'+1
}
est restore other

    predict hr_gender${age},   hrnum(female 1 agediag ${age}) ///
					         hrdenom(female 0 agediag ${age} ) timevar(tt1) 
	local p=1 
	foreach x in $times {
	local hr_gender_at`x'= hr_gender${age}[`p']
	local p=`p'+1
	}
	
	local givenage = `ageloop'

// store results


		post sim_t1t1_linear (`CIF1_f1_at1') (`CIF1_f1_at2') (`CIF1_f1_at3') (`CIF1_f1_at4') (`CIF1_f1_at5') (`CIF1_f1_at6') (`CIF1_f1_at7') (`CIF1_f1_at8') (`CIF1_f1_at9') (`CIF1_f1_at10') ///
					  (`CIF1_f1_at1_LL') (`CIF1_f1_at1_UL') (`CIF1_f1_at2_LL') (`CIF1_f1_at2_UL') (`CIF1_f1_at3_LL') (`CIF1_f1_at3_UL') (`CIF1_f1_at4_LL') (`CIF1_f1_at4_UL') (`CIF1_f1_at5_LL') (`CIF1_f1_at5_UL') ///
					  (`CIF1_f1_at6_LL') (`CIF1_f1_at6_UL') (`CIF1_f1_at7_LL') (`CIF1_f1_at7_UL') (`CIF1_f1_at8_LL') (`CIF1_f1_at8_UL') (`CIF1_f1_at9_LL') (`CIF1_f1_at9_UL') (`CIF1_f1_at10_LL') (`CIF1_f1_at10_UL') ///
					  (`CIF2_f1_at1') (`CIF2_f1_at2') (`CIF2_f1_at3') (`CIF2_f1_at4') (`CIF2_f1_at5') (`CIF2_f1_at6') (`CIF2_f1_at7') (`CIF2_f1_at8') (`CIF2_f1_at9') (`CIF2_f1_at10') ///
					  (`CIF2_f1_at1_LL') (`CIF2_f1_at1_UL') (`CIF2_f1_at2_LL') (`CIF2_f1_at2_UL') (`CIF2_f1_at3_LL') (`CIF2_f1_at3_UL') (`CIF2_f1_at4_LL') (`CIF2_f1_at4_UL') (`CIF2_f1_at5_LL') (`CIF2_f1_at5_UL') ///
					  (`CIF2_f1_at6_LL') (`CIF2_f1_at6_UL') (`CIF2_f1_at7_LL') (`CIF2_f1_at7_UL') (`CIF2_f1_at8_LL') (`CIF2_f1_at8_UL') (`CIF2_f1_at9_LL') (`CIF2_f1_at9_UL') (`CIF2_f1_at10_LL') (`CIF2_f1_at10_UL') ///
					  ///
					  (`CIF1_f0_at1') (`CIF1_f0_at2') (`CIF1_f0_at3') (`CIF1_f0_at4') (`CIF1_f0_at5') (`CIF1_f0_at6') (`CIF1_f0_at7') (`CIF1_f0_at8') (`CIF1_f0_at9') (`CIF1_f0_at10') ///
					  (`CIF1_f0_at1_LL') (`CIF1_f0_at1_UL') (`CIF1_f0_at2_LL') (`CIF1_f0_at2_UL') (`CIF1_f0_at3_LL') (`CIF1_f0_at3_UL') (`CIF1_f0_at4_LL') (`CIF1_f0_at4_UL') (`CIF1_f0_at5_LL') (`CIF1_f0_at5_UL') ///
					  (`CIF1_f0_at6_LL') (`CIF1_f0_at6_UL') (`CIF1_f0_at7_LL') (`CIF1_f0_at7_UL') (`CIF1_f0_at8_LL') (`CIF1_f0_at8_UL') (`CIF1_f0_at9_LL') (`CIF1_f0_at9_UL') (`CIF1_f0_at10_LL') (`CIF1_f0_at10_UL') ///
					  (`CIF2_f0_at1') (`CIF2_f0_at2') (`CIF2_f0_at3') (`CIF2_f0_at4') (`CIF2_f0_at5') (`CIF2_f0_at6') (`CIF2_f0_at7') (`CIF2_f0_at8') (`CIF2_f0_at9') (`CIF2_f0_at10') ///
					  (`CIF2_f0_at1_LL') (`CIF2_f0_at1_UL') (`CIF2_f0_at2_LL') (`CIF2_f0_at2_UL') (`CIF2_f0_at3_LL') (`CIF2_f0_at3_UL') (`CIF2_f0_at4_LL') (`CIF2_f0_at4_UL') (`CIF2_f0_at5_LL') (`CIF2_f0_at5_UL') ///
					  (`CIF2_f0_at6_LL') (`CIF2_f0_at6_UL') (`CIF2_f0_at7_LL') (`CIF2_f0_at7_UL') (`CIF2_f0_at8_LL') (`CIF2_f0_at8_UL') (`CIF2_f0_at9_LL') (`CIF2_f0_at9_UL') (`CIF2_f0_at10_LL') (`CIF2_f0_at10_UL') ///
					  ///
					  (`CIF1_diff_at1') (`CIF1_diff_at2') (`CIF1_diff_at3') (`CIF1_diff_at4') (`CIF1_diff_at5') (`CIF1_diff_at6') (`CIF1_diff_at7') (`CIF1_diff_at8') (`CIF1_diff_at9') (`CIF1_diff_at10') ///
					  (`CIF2_diff_at1') (`CIF2_diff_at2') (`CIF2_diff_at3') (`CIF2_diff_at4') (`CIF2_diff_at5') (`CIF2_diff_at6') (`CIF2_diff_at7') (`CIF2_diff_at8') (`CIF2_diff_at9') (`CIF2_diff_at10') ///
					  ///
					  (`CIF1_ratio_at1') (`CIF1_ratio_at2') (`CIF1_ratio_at3') (`CIF1_ratio_at4') (`CIF1_ratio_at5') (`CIF1_ratio_at6') (`CIF1_ratio_at7') (`CIF1_ratio_at8') (`CIF1_ratio_at9') (`CIF1_ratio_at10') ///
					  (`CIF2_ratio_at1') (`CIF2_ratio_at2') (`CIF2_ratio_at3') (`CIF2_ratio_at4') (`CIF2_ratio_at5') (`CIF2_ratio_at6') (`CIF2_ratio_at7') (`CIF2_ratio_at8') (`CIF2_ratio_at9') (`CIF2_ratio_at10') ///
					  ///
					  (`CIF_total_at1') (`CIF_total_at2') (`CIF_total_at3') (`CIF_total_at4') (`CIF_total_at5') (`CIF_total_at6') (`CIF_total_at7') (`CIF_total_at8') (`CIF_total_at9') (`CIF_total_at10') /// 
					  ///
					  (`hr_gender_at1') (`hr_gender_at2') (`hr_gender_at3') (`hr_gender_at4') (`hr_gender_at5') ///
					  (`hr_gender_at6') (`hr_gender_at7') (`hr_gender_at8') (`hr_gender_at9') (`hr_gender_at10') (`givenage') ///
					  (`canc_model_converged') (`other_model_converged')
					  
	}	
	}
	
	
	
	}
	// close postfile
	postclose sim_t1t1_linear
}


