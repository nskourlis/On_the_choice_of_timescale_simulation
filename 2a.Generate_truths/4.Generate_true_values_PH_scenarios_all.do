// Generate Truths and store as globals in separate file

 cd "$N\1b.Scenarios_folder"
 
 //agemean is the age upon which age at diagnosis is centered when modelling its effect on cause I
 global agemean= 65
 
/************************************************************************************/
/* loop over scenarios with proportional hazards of gender for other cause mortality*/
/************************************************************************************/

foreach scen in  1 3 5 7 9 11 13 15 17 19 21 23 {
	clear
	
	/* load scenario parameters*/
	include "$N\1b.Scenarios_folder/Scenario`scen'_new.do"
	
	/* At the end we will keep the estimates only for the list of time since diagnosis points specified by global macro times*/
	global  times "1 2 3 4 5 6 7 8 9 10"
	
   /*Loop for deriving the true cause specific hazards, survival*/ 
   /*and CIFs based on the scenario parameters over ages at diagnosis 50,60,70,80,90*/ 
	
	forvalues ageloop= 50(10)90 {
	global age=`ageloop'

	global female=1
    
	//10 years of maximum follow up for each individual
	global interval =10 
	
	//Generate timevars for the different timescales. Time "runs" the same for both.
	//Note: Attained age= age at diagnosis + time since diagnosis 
	
	//Time dependent effect of age at diagnosis to time since diagnosis

	range tt1 0 $interval 3001    									//Time since diagnosis timevar
    range tt2_`ageloop' `ageloop' `ageloop'+$interval 3001  	    //Attained age timevar
	
    //Hazard and survival functions for males and females for cancer mortality

		gen h1_f1=(((${lambda1} *${gamma1}* tt1^( ${gamma1} -1 )*${pmix} * exp(-${lambda1}  *tt1^(${gamma1} ) ) )+ (${lambda2}* ${gamma2}* (tt1^(${gamma2} -1) )* ///
	       (1-${pmix} )*exp(-${lambda2} *tt1^( ${gamma2} ) ))) /(${pmix} *exp(- ${lambda1} *tt1^( ${gamma1} ) ) + ///
		   (1-${pmix} )*exp(-${lambda2} *tt1^( ${gamma2} ) ))) *exp((${b1}*(`ageloop'-${agemean})+ (${b2}*(`ageloop'-${agemean})^2) +  ${bsex_c1}*${female} ) )  
		   
		gen h1_f0=(((${lambda1} *${gamma1}* tt1^( ${gamma1} -1 )*${pmix} * exp(-${lambda1}  *tt1^(${gamma1} ) ) )+ (${lambda2}* ${gamma2}* (tt1^(${gamma2} -1) )* ///
	       (1-${pmix} )*exp(-${lambda2} *tt1^( ${gamma2} ) ))) /(${pmix} *exp(- ${lambda1} *tt1^( ${gamma1} ) ) + ///
		   (1-${pmix} )*exp(-${lambda2} *tt1^( ${gamma2} ) )))*exp((${b1}*(`ageloop'-${agemean})+ (${b2}*(`ageloop'-${agemean})^2)))
		   
		   
		gen s1_f1=(${pmix} *exp(-${lambda1} *tt1^(${gamma1}) )+(1-${pmix} )*exp(-${lambda2} *tt1^(${gamma2}) )) ///
	        ^exp((${b1}*(`ageloop'-${agemean})+ (${b2}*(`ageloop'-${agemean})^2) +  ${bsex_c1}*${female}) ) 
			
		gen s1_f0=(${pmix} *exp(-${lambda1} *tt1^(${gamma1}) )+(1-${pmix} )*exp(-${lambda2} *tt1^(${gamma2}) ))^exp((${b1}*(`ageloop'-${agemean})+ (${b2}*(`ageloop'-${agemean})^2))) 
	        
	
    //Hazard and survival functions for males and females for other cause mortality
	//on the hazard form for other cause mortality
	
	if "$base_function"=="Weibull" {
	
		gen h2_f1=(${alpha}+ ${lambda} *${gamma} *(`ageloop'+tt1)^(${gamma} -1)) *exp(${bsex_c2}*${female}) // + ${bsex_c2_tvc_1}*${female}*(`ageloop'+tt1) + ${bsex_c2_tvc_2}*${female}*(`ageloop'+tt1)^2  )

		gen s2_f1=(exp(-(${alpha}*(`ageloop'+tt1)+${lambda} *(`ageloop'+tt1)^${gamma} )))^exp(${bsex_c2}*${female}) // + ${bsex_c2_tvc_1}*${female}*(`ageloop'+tt1) + ${bsex_c2_tvc_2}*${female}*(`ageloop'+tt1)^2  )
	
		gen s2_age_f1 =exp(-(${alpha}*`ageloop'+${lambda} *`ageloop'^${gamma}))^exp(${bsex_c2}*${female}) // + ${bsex_c2_tvc_1}*${female}*(`ageloop') + ${bsex_c2_tvc_2}*${female}*(`ageloop')^2  )
   
		gen s2c_f1=s2_f1/s2_age_f1
	

		gen h2_f0=(${alpha}+ ${lambda} *${gamma} *(`ageloop'+tt1)^(${gamma} -1)) 
	
		gen s2_f0=(exp(-(${alpha}*(`ageloop'+tt1)+${lambda} *(`ageloop'+tt1)^${gamma} )))
	
		gen s2_age_f0 =exp(-(${alpha}*`ageloop'+${lambda} *`ageloop'^${gamma}))
	
		gen s2c_f0=s2_f0/s2_age_f0
		
		}

	//We can comment out the gender time interaction termas for the scenarios that we know PH for gender hold

	else if  "$base_function"=="Polynomial" {
	
		gen h2_f1= ((2*${c1}*(`ageloop'+tt1)+${c2}) * exp(${c1}* ((`ageloop'+tt1)^2)+ ${c2}* (`ageloop'+tt1)+  ${c3} ))* ///
			 exp(${bsex_c2}*${female}) // + ${bsex_c2_tvc_1}*${female}*(`ageloop'+tt1) + ${bsex_c2_tvc_2}*${female}*(`ageloop'+tt1)^2  )
	
		gen s2_f1= (exp(-(exp(${c1}* ((`ageloop'+tt1)^2)+ ${c2}* (`ageloop'+tt1)+  ${c3} ) ) ))^ ///
			 exp(${bsex_c2}*${female}) // + ${bsex_c2_tvc_1}*${female}*(`ageloop'+tt1) + ${bsex_c2_tvc_2}*${female}*(`ageloop'+tt1)^2  )
    
		gen s2_age_f1 =  (exp(-(exp(${c1}* ((`ageloop')^2)+ ${c2}* (`ageloop')+  ${c3} ) ) ))^ ///
					exp(${bsex_c2}*${female}) // + ${bsex_c2_tvc_1}*${female}*(`ageloop') + ${bsex_c2_tvc_2}*${female}*(`ageloop')^2  )
	
		gen s2c_f1=s2_f1/s2_age_f1
	
	
		gen h2_f0= ((2*${c1}*(`ageloop'+tt1)+${c2}) * exp(${c1}* ((`ageloop'+tt1)^2)+ ${c2}* (`ageloop'+tt1)+  ${c3} ))
	
		gen s2_f0= (exp(-(exp(${c1}* ((`ageloop'+tt1)^2)+ ${c2}* (`ageloop'+tt1)+  ${c3} ) ) ))
  
		gen s2_age_f0 =  (exp(-(exp(${c1}* ((`ageloop')^2)+ ${c2}* (`ageloop')+  ${c3} ) ) ))
	
		gen s2c_f0=s2_f0/s2_age_f0

	      }
		  
    else if  "$base_function"=="Gompertz" {
	
		gen h2_f1=(${alpha} + ${lambda}* exp(${gamma}*(`ageloop'+tt1)))*exp(${bsex_c2}*${female}) //+ ${bsex_k2_tvc}*${female}*(`ageloop'+tt1)) 
	
		gen s2_f1=  (exp( - ((${alpha}*(`ageloop'+tt1))+ (${lambda}/${gamma}) *(exp(${gamma}*(`ageloop'+tt1))-1) )))^exp(${bsex_c2}*${female}) //+ ${bsex_k2_tvc}*${female}*(`ageloop'+tt1))
	
		gen s2_age_f1 = (exp( - ((${alpha}*(`ageloop'))+ (${lambda}/${gamma}) *(exp(${gamma}*(`ageloop'))-1) )))^exp(${bsex_c2}*${female}) //+ ${bsex_k2_tvc}*${female}*`ageloop')
	
		gen s2c_f1=s2_f1/s2_age_f1
	
	
		gen h2_f0=(${alpha} + ${lambda}* exp(${gamma}*(`ageloop'+tt1)))
	
		gen s2_f0=  (exp( - ((${alpha}*(`ageloop'+tt1))+ (${lambda}/${gamma}) *(exp(${gamma}*(`ageloop'+tt1))-1) )))
	
		gen s2_age_f0 = (exp( - ((${alpha}*(`ageloop'))+ (${lambda}/${gamma}) *(exp(${gamma}*(`ageloop'))-1) )))
	
		gen s2c_f0=s2_f0/s2_age_f0
	
	}
	
	gen hr_gender= exp(${bsex_c2}*${female}) 


	/*Estimating the CIFs*/
	
	// Females
	gen integrand1_f1=s1_f1*s2c_f1*h1_f1
	integ integrand1_f1 tt1, gen(CIF1_f1)

	//Change tt2_$age to tt1 to all the expressions as CIFs are always integrated  between 0 and t for each a0
	gen integrand2_f1=s1_f1*s2c_f1*h2_f1
	integ integrand2_f1 tt1, gen(CIF2_f1) 

	//Males
	gen integrand1_f0=s1_f0*s2c_f0*h1_f0
	integ integrand1_f0 tt1, gen(CIF1_f0)

	//Change tt2_$age to tt1 to all the expressions as CIFs are always integrated  between 0 and t for each a0
	gen integrand2_f0=s1_f0*s2c_f0*h2_f0
	integ integrand2_f0 tt1, gen(CIF2_f0) 
	
	keep tt1 tt2_`ageloop'  h1_f1 s1_f1 h2_f1 s2_f1 s2c_f1 CIF1_f1 CIF2_f1 ///
	hr_gender h1_f0 s1_f0 h2_f0 s2_f0 s2c_f0 CIF1_f0 CIF2_f0
	
	
	/*keep the estimates only for the list of time points specified by global macro times*/
	
	foreach x in $times {
	  preserve
		keep if tt1==`x'
		local CIF1_f1_at`x' = CIF1_f1
		local CIF2_f1_at`x' = CIF2_f1
	
		local CIF1_f0_at`x' = CIF1_f0
		local CIF2_f0_at`x' = CIF2_f0
	
		local CIF1_diff_at`x'=CIF1_f1-CIF1_f0
		local CIF2_diff_at`x'=CIF2_f1-CIF2_f0
	
		local CIF1_ratio_at`x'=CIF1_f1/CIF1_f0
		local CIF2_ratio_at`x'=CIF2_f1/CIF2_f0
	
		local CIF_total_at`x'=(CIF1_f1+CIF2_f1)/(CIF1_f0+CIF2_f0)

	  restore
	}
	
	
	foreach x in $times {
	  preserve
		keep if tt1==`x'
		local hr_gender_at`x' = hr_gender

	  restore
	}
	

	   local bsex_c1=$bsex_c1

       local bsex_c2=$bsex_c2

	   local bsex_c2_tvc_1 = $bsex_c2_tvc_1
	   
	   local bsex_c2_tvc_2 = $bsex_c2_tvc_2
	   
	   local prop_sex_other= "$prop_sex_other"
	   
	   local base_function="$base_function"
	   
	   local independence=$independence

	   local obs= $obs
	   
	   local sdage= $sdage


	// write true values to folder 2b.Truths_folder
	
	file open truth using "$N\2b.Truths_folder\scenario`scen'_Truth_`ageloop'.do", write replace
	
	file write truth "// Scenario`scen'_`ageloop'" _newline ////
				     "global CIF1_f1_at1_true 		`CIF1_f1_at1'" _newline ///
				     "global CIF2_f1_at1_true 		`CIF2_f1_at1'" _newline ///
					 "global CIF1_f1_at2_true 		`CIF1_f1_at2'" _newline ///
				     "global CIF2_f1_at2_true 		`CIF2_f1_at2'" _newline ///
					 "global CIF1_f1_at3_true 		`CIF1_f1_at3'" _newline ///
				     "global CIF2_f1_at3_true 		`CIF2_f1_at3'" _newline ///
					 "global CIF1_f1_at4_true 		`CIF1_f1_at4'" _newline ///
				     "global CIF2_f1_at4_true 		`CIF2_f1_at4'" _newline ///
					 "global CIF1_f1_at5_true 		`CIF1_f1_at5'" _newline ///
				     "global CIF2_f1_at5_true 		`CIF2_f1_at5'" _newline ///
					 "global CIF1_f1_at6_true 		`CIF1_f1_at6'" _newline ///
				     "global CIF2_f1_at6_true 		`CIF2_f1_at6'" _newline ///
					 "global CIF1_f1_at7_true 		`CIF1_f1_at7'" _newline ///
				     "global CIF2_f1_at7_true 		`CIF2_f1_at7'" _newline ///
					 "global CIF1_f1_at8_true 		`CIF1_f1_at8'" _newline ///
				     "global CIF2_f1_at8_true 		`CIF2_f1_at8'" _newline ///
					 "global CIF1_f1_at9_true 		`CIF1_f1_at9'" _newline ///
				     "global CIF2_f1_at9_true 		`CIF2_f1_at9'" _newline ///
					 "global CIF1_f1_at10_true	    `CIF1_f1_at10'" _newline ///
				     "global CIF2_f1_at10_true 	    `CIF2_f1_at10'" _newline ///
					 ///
					 "global CIF1_f0_at1_true 		`CIF1_f0_at1'" _newline ///
				     "global CIF2_f0_at1_true 		`CIF2_f0_at1'" _newline ///
					 "global CIF1_f0_at2_true 		`CIF1_f0_at2'" _newline ///
				     "global CIF2_f0_at2_true 		`CIF2_f0_at2'" _newline ///
					 "global CIF1_f0_at3_true 		`CIF1_f0_at3'" _newline ///
				     "global CIF2_f0_at3_true 		`CIF2_f0_at3'" _newline ///
					 "global CIF1_f0_at4_true 		`CIF1_f0_at4'" _newline ///
				     "global CIF2_f0_at4_true 		`CIF2_f0_at4'" _newline ///
					 "global CIF1_f0_at5_true 		`CIF1_f0_at5'" _newline ///
				     "global CIF2_f0_at5_true 		`CIF2_f0_at5'" _newline ///
					 "global CIF1_f0_at6_true 		`CIF1_f0_at6'" _newline ///
				     "global CIF2_f0_at6_true 		`CIF2_f0_at6'" _newline ///
					 "global CIF1_f0_at7_true 		`CIF1_f0_at7'" _newline ///
				     "global CIF2_f0_at7_true 		`CIF2_f0_at7'" _newline ///
					 "global CIF1_f0_at8_true 		`CIF1_f0_at8'" _newline ///
				     "global CIF2_f0_at8_true 		`CIF2_f0_at8'" _newline ///
					 "global CIF1_f0_at9_true 		`CIF1_f0_at9'" _newline ///
				     "global CIF2_f0_at9_true 		`CIF2_f0_at9'" _newline ///
					 "global CIF1_f0_at10_true	    `CIF1_f0_at10'" _newline ///
				     "global CIF2_f0_at10_true 	    `CIF2_f0_at10'" _newline ///
					 ///
					 "global CIF1_diff_at1_true 		`CIF1_diff_at1'" _newline ///
				     "global CIF2_diff_at1_true 		`CIF2_diff_at1'" _newline ///
					 "global CIF1_diff_at2_true 		`CIF1_diff_at2'" _newline ///
				     "global CIF2_diff_at2_true 		`CIF2_diff_at2'" _newline ///
					 "global CIF1_diff_at3_true 		`CIF1_diff_at3'" _newline ///
				     "global CIF2_diff_at3_true 		`CIF2_diff_at3'" _newline ///
					 "global CIF1_diff_at4_true 		`CIF1_diff_at4'" _newline ///
				     "global CIF2_diff_at4_true 		`CIF2_diff_at4'" _newline ///
					 "global CIF1_diff_at5_true 		`CIF1_diff_at5'" _newline ///
				     "global CIF2_diff_at5_true 		`CIF2_diff_at5'" _newline ///
					 "global CIF1_diff_at6_true 		`CIF1_diff_at6'" _newline ///
				     "global CIF2_diff_at6_true 		`CIF2_diff_at6'" _newline ///
					 "global CIF1_diff_at7_true 		`CIF1_diff_at7'" _newline ///
				     "global CIF2_diff_at7_true 		`CIF2_diff_at7'" _newline ///
					 "global CIF1_diff_at8_true 		`CIF1_diff_at8'" _newline ///
				     "global CIF2_diff_at8_true 		`CIF2_diff_at8'" _newline ///
					 "global CIF1_diff_at9_true 		`CIF1_diff_at9'" _newline ///
				     "global CIF2_diff_at9_true 		`CIF2_diff_at9'" _newline ///
					 "global CIF1_diff_at10_true	    `CIF1_diff_at10'" _newline ///
				     "global CIF2_diff_at10_true 	    `CIF2_diff_at10'" _newline ///
					 ///
					 "global CIF1_ratio_at1_true 		`CIF1_ratio_at1'" _newline ///
				     "global CIF2_ratio_at1_true 		`CIF2_ratio_at1'" _newline ///
					 "global CIF1_ratio_at2_true 		`CIF1_ratio_at2'" _newline ///
				     "global CIF2_ratio_at2_true 		`CIF2_ratio_at2'" _newline ///
					 "global CIF1_ratio_at3_true 		`CIF1_ratio_at3'" _newline ///
				     "global CIF2_ratio_at3_true 		`CIF2_ratio_at3'" _newline ///
					 "global CIF1_ratio_at4_true 		`CIF1_ratio_at4'" _newline ///
				     "global CIF2_ratio_at4_true 		`CIF2_ratio_at4'" _newline ///
					 "global CIF1_ratio_at5_true 		`CIF1_ratio_at5'" _newline ///
				     "global CIF2_ratio_at5_true 		`CIF2_ratio_at5'" _newline ///
					 "global CIF1_ratio_at6_true 		`CIF1_ratio_at6'" _newline ///
				     "global CIF2_ratio_at6_true 		`CIF2_ratio_at6'" _newline ///
					 "global CIF1_ratio_at7_true 		`CIF1_ratio_at7'" _newline ///
				     "global CIF2_ratio_at7_true 		`CIF2_ratio_at7'" _newline ///
					 "global CIF1_ratio_at8_true 		`CIF1_ratio_at8'" _newline ///
				     "global CIF2_ratio_at8_true 		`CIF2_ratio_at8'" _newline ///
					 "global CIF1_ratio_at9_true 		`CIF1_ratio_at9'" _newline ///
				     "global CIF2_ratio_at9_true 		`CIF2_ratio_at9'" _newline ///
					 "global CIF1_ratio_at10_true	    `CIF1_ratio_at10'" _newline ///
				     "global CIF2_ratio_at10_true 	    `CIF2_ratio_at10'" _newline ///
					 ///
					 "global CIF_total_at1_true 		`CIF_total_at1'" _newline ///
				     "global CIF_total_at2_true 		`CIF_total_at2'" _newline ///
					 "global CIF_total_at3_true 		`CIF_total_at3'" _newline ///
				     "global CIF_total_at4_true 		`CIF_total_at4'" _newline ///
					 "global CIF_total_at5_true 		`CIF_total_at5'" _newline ///
				     "global CIF_total_at6_true 		`CIF_total_at6'" _newline ///
					 "global CIF_total_at7_true 		`CIF_total_at7'" _newline ///
				     "global CIF_total_at8_true 		`CIF_total_at8'" _newline ///
					 "global CIF_total_at9_true	        `CIF_total_at9'" _newline ///
				     "global CIF_total_at10_true 	    `CIF_total_at10'" _newline ///
					///
				     "global hr_gender_at1_true `hr_gender_at1'" _newline ///
				     "global hr_gender_at2_true `hr_gender_at2'" _newline ///
				     "global hr_gender_at3_true `hr_gender_at3'" _newline ///
				     "global hr_gender_at4_true `hr_gender_at4'" _newline ///
				     "global hr_gender_at5_true `hr_gender_at5'" _newline ///
				     "global hr_gender_at6_true `hr_gender_at6'" _newline ///
				     "global hr_gender_at7_true `hr_gender_at7'" _newline ///
				     "global hr_gender_at8_true `hr_gender_at8'" _newline ///
				     "global hr_gender_at9_true `hr_gender_at9'" _newline ///
				     "global hr_gender_at10_true `hr_gender_at10'" _newline ///
					 ///
                     "global bsex_c1   			`bsex_c1'" _newline ///
                     "global bsex_c2   			`bsex_c2'" _newline ///
					 "global bsex_c2_tvc_1   		`bsex_c2_tvc_1'" _newline ///
					 "global bsex_c2_tvc_2   		`bsex_c2_tvc_2'" _newline ///
					 "global prop_sex_other   	`prop_sex_other'" _newline ///
					 "global base_function   	`base_function'" _newline ///
					 "global independence   	`independence'" _newline ///
					 "global obs                 `obs'"         _newline ///
					 "global sdage                 `sdage'"         _newline 
					 


	file close truth

	// save "$N\Nick_Simulation\Simulation_mod\2b.Truths_folder\Truth_dataset`scen'_`ageloop'.dta",replace
   drop tt1 h1_f1 s1_f1 h2_f1 s2_f1 s2c_f1 CIF1_f1 CIF2_f1  h1_f0 s1_f0 h2_f0 s2_f0 s2c_f0 CIF1_f0 CIF2_f0  tt2_`ageloop' hr_gender
   
   }
   }
   

 /************************************************************************************/
/* loop over scenarios with non proportional hazards of gender for other cause mortality*/
/************************************************************************************/

foreach scen in 2 4 6 8 10 12 14 16 18 20 22 24  {
	clear
	
	//agemean is the age upon which age at diagnosis is centered when modelling its effect on cause I
    global agemean= 65
	
	/* load scenario parameters*/
	include "$N\1b.Scenarios_folder/Scenario`scen'_new.do"
	
	/* At the end we will keep the estimates only for the list of time since diagnosis points specified by global macro times*/
	global  times "1 2 3 4 5 6 7 8 9 10"
	
   /*Loop for deriving the true cause specific hazards, survival*/ 
   /*and CIFs based on the scenario parameters over ages at diagnosis 50,60,70,80*/ 
   
	
	/*********************************************/
	//Estimate true CIFs for men from analytic types
	/***********************************************/
	
	forvalues ageloop= 50(10)90 {
	
	clear

	global age=`ageloop'
	
	//10 years of maximum follow up for each individual
    global interval =10 
	
	global female=1
	
	//Generate timevars for the different timescales. Time "runs" the same for both. Note: Attained age= age at diagnosis + time since diagnosis 
	
	range tt1 0 $interval 3001    									//Time since diagnosis timevar
    range tt2_`ageloop' `ageloop' `ageloop'+$interval 3001  	    //Attained age timevar
	
	//Hazard and survival functions for males and females for cancer mortality
   
	gen h1_f0=(((${lambda1} *${gamma1}* tt1^( ${gamma1} -1 )*${pmix} * exp(-${lambda1}  *tt1^(${gamma1} ) ) )+ (${lambda2}* ${gamma2}* (tt1^(${gamma2} -1) )* ///
	       (1-${pmix} )*exp(-${lambda2} *tt1^( ${gamma2} ) ))) /(${pmix} *exp(- ${lambda1} *tt1^( ${gamma1} ) ) + ///
		   (1-${pmix} )*exp(-${lambda2} *tt1^( ${gamma2} ) )))*exp(${b1}*(`ageloop'-${agemean})+ (${b2}*(`ageloop'-${agemean})^2))
		   
	gen s1_f0=(${pmix} *exp(-${lambda1} *tt1^(${gamma1}) )+(1-${pmix} )*exp(-${lambda2} *tt1^(${gamma2}) ))^exp(${b1}*(`ageloop'-${agemean})+ (${b2}*(`ageloop'-${agemean})^2)) 
	 
	//Hazard and survival functions for males and females for other cause mortality
	//on the hazard form for other cause mortality
	
	if "$base_function"=="Weibull" {
	
		gen h2_f0=(${alpha}+ ${lambda} *${gamma} *(`ageloop'+tt1)^(${gamma} -1)) 
	
		gen s2_f0=(exp(-(${alpha}*(`ageloop'+tt1)+${lambda} *(`ageloop'+tt1)^${gamma} )))
	
		gen s2_age_f0 =exp(-(${alpha}*`ageloop'+${lambda} *`ageloop'^${gamma}))
	
		gen s2c_f0=s2_f0/s2_age_f0
	
		}
		
	else if  "$base_function"=="Polynomial" {

		gen h2_f0= ((2*${c1}*(`ageloop'+tt1)+${c2}) * exp(${c1}* ((`ageloop'+tt1)^2)+ ${c2}* (`ageloop'+tt1)+  ${c3} ))
		
		gen s2_f0= (exp(-(exp(${c1}* ((`ageloop'+tt1)^2)+ ${c2}* (`ageloop'+tt1)+  ${c3} ) ) ))
		    
		gen s2_age_f0 =  (exp(-(exp(${c1}* ((`ageloop')^2)+ ${c2}* (`ageloop')+  ${c3} ) ) ))
	
		gen s2c_f0=s2_f0/s2_age_f0

	  }
	
			  
    else if  "$base_function"=="Gompertz" {
	
		gen h2_f0=(${alpha} + ${lambda}* exp(${gamma}*(`ageloop'+tt1)))
	
		gen s2_f0=  (exp( - ((${alpha}*(`ageloop'+tt1))+ (${lambda}/${gamma}) *(exp(${gamma}*(`ageloop'+tt1))-1) )))
	
		gen s2_age_f0 = (exp( - ((${alpha}*(`ageloop'))+ (${lambda}/${gamma}) *(exp(${gamma}*(`ageloop'))-1) )))
	
		gen s2c_f0=s2_f0/s2_age_f0
	
	  }
	
	
	gen hr_gender= exp(${bsex_c2}*${female} + ${bsex_c2_tvc_1}*${female}*(`ageloop'+tt1) + ${bsex_c2_tvc_2}*${female}*(`ageloop'+tt1)^2  )


	/*Estimating the CIFs*/
	
	//Males
	gen integrand1_f0=s1_f0*s2c_f0*h1_f0
	integ integrand1_f0 tt1, gen(CIF1_f0)

	//Change tt2_$age to tt1 to all the expressions as CIFs are always integrated  between 0 and t for each a0
	gen integrand2_f0=s1_f0*s2c_f0*h2_f0
	integ integrand2_f0 tt1, gen(CIF2_f0) 

	keep tt1 tt2_`ageloop' hr_gender h1_f0 s1_f0 h2_f0 s2_f0 s2c_f0 CIF1_f0 CIF2_f0

	/*keep the estimates only for the list of time points specified by global macro times*/

	foreach x in $times {
	   preserve
		keep if tt1==`x'
		local CIF1_f0_mata_at`x' = CIF1_f0
		local CIF2_f0_mata_at`x' = CIF2_f0
	   restore
	}
	
	foreach x in $times {
	  preserve
		keep if tt1==`x'
		local hr_gender_at`x' = hr_gender
	  restore
	}
	
//Pass on some of the scenarios basic parameters to the "Truth" file
	
	   local bsex_c1=$bsex_c1

       local bsex_c2=$bsex_c2

	   local bsex_c2_tvc_1 = $bsex_c2_tvc_1
	   
	   local bsex_c2_tvc_2 = $bsex_c2_tvc_2
	   
	   local prop_sex_other= "$prop_sex_other"
	   
	   local base_function="$base_function"
	   
	   local independence=$independence
	   
	   local obs= $obs
	   
	   local sdage= $sdage

	   
	  
	 /***********************************************************************/
	//Import the true CIFs for women from double integration mata outputs
	/************************************************************************/
	   
	use "$N\2b.Truths_folder/True_CIFs_scenario`scen'", replace
	keep if givenage==`ageloop'
	
	/*keep the estimates only for the list of time points specified by global macro times*/

	foreach x in $times {
	  preserve
		keep if year==`x'  
		local CIF1_f1_mata_at`x' = CIF1_true_f1_at
		local CIF2_f1_mata_at`x' = CIF2_true_f1_at
	
		local CIF1_diff_at`x'=CIF1_true_f1_at-`CIF1_f0_mata_at`x''
		local CIF2_diff_at`x'=CIF2_true_f1_at-`CIF2_f0_mata_at`x''
	
		local CIF1_ratio_at`x'=CIF1_true_f1_at/`CIF1_f0_mata_at`x''
		local CIF2_ratio_at`x'=CIF2_true_f1_at/`CIF2_f0_mata_at`x''
	
		local CIF_total_at`x'=(CIF1_true_f1_at+CIF2_true_f1_at)/(`CIF1_f0_mata_at`x''+`CIF2_f0_mata_at`x'')
      restore
	}
	

	// write true values to folder 2b.Truths_folder
	file open truth using "$N\2b.Truths_folder\scenario`scen'_Truth_`ageloop'.do", write replace
	
	file write truth "// Scenario`scen'_`ageloop'" _newline ////
				     "global CIF1_f1_at1_true 		`CIF1_f1_mata_at1'" _newline ///
				     "global CIF2_f1_at1_true 		`CIF2_f1_mata_at1'" _newline ///
					 "global CIF1_f1_at2_true 		`CIF1_f1_mata_at2'" _newline ///
				     "global CIF2_f1_at2_true 		`CIF2_f1_mata_at2'" _newline ///
					 "global CIF1_f1_at3_true 		`CIF1_f1_mata_at3'" _newline ///
				     "global CIF2_f1_at3_true 		`CIF2_f1_mata_at3'" _newline ///
					 "global CIF1_f1_at4_true 		`CIF1_f1_mata_at4'" _newline ///
				     "global CIF2_f1_at4_true 		`CIF2_f1_mata_at4'" _newline ///
					 "global CIF1_f1_at5_true 		`CIF1_f1_mata_at5'" _newline ///
				     "global CIF2_f1_at5_true 		`CIF2_f1_mata_at5'" _newline ///
					 "global CIF1_f1_at6_true 		`CIF1_f1_mata_at6'" _newline ///
				     "global CIF2_f1_at6_true 		`CIF2_f1_mata_at6'" _newline ///
					 "global CIF1_f1_at7_true 		`CIF1_f1_mata_at7'" _newline ///
				     "global CIF2_f1_at7_true 		`CIF2_f1_mata_at7'" _newline ///
					 "global CIF1_f1_at8_true 		`CIF1_f1_mata_at8'" _newline ///
				     "global CIF2_f1_at8_true 		`CIF2_f1_mata_at8'" _newline ///
					 "global CIF1_f1_at9_true 		`CIF1_f1_mata_at9'" _newline ///
				     "global CIF2_f1_at9_true 		`CIF2_f1_mata_at9'" _newline ///
					 "global CIF1_f1_at10_true	    `CIF1_f1_mata_at10'" _newline ///
				     "global CIF2_f1_at10_true 	    `CIF2_f1_mata_at10'" _newline ///
					 ///
					 "global CIF1_f0_at1_true 		`CIF1_f0_mata_at1'" _newline ///
				     "global CIF2_f0_at1_true 		`CIF2_f0_mata_at1'" _newline ///
					 "global CIF1_f0_at2_true 		`CIF1_f0_mata_at2'" _newline ///
				     "global CIF2_f0_at2_true 		`CIF2_f0_mata_at2'" _newline ///
					 "global CIF1_f0_at3_true 		`CIF1_f0_mata_at3'" _newline ///
				     "global CIF2_f0_at3_true 		`CIF2_f0_mata_at3'" _newline ///
					 "global CIF1_f0_at4_true 		`CIF1_f0_mata_at4'" _newline ///
				     "global CIF2_f0_at4_true 		`CIF2_f0_mata_at4'" _newline ///
					 "global CIF1_f0_at5_true 		`CIF1_f0_mata_at5'" _newline ///
				     "global CIF2_f0_at5_true 		`CIF2_f0_mata_at5'" _newline ///
					 "global CIF1_f0_at6_true 		`CIF1_f0_mata_at6'" _newline ///
				     "global CIF2_f0_at6_true 		`CIF2_f0_mata_at6'" _newline ///
					 "global CIF1_f0_at7_true 		`CIF1_f0_mata_at7'" _newline ///
				     "global CIF2_f0_at7_true 		`CIF2_f0_mata_at7'" _newline ///
					 "global CIF1_f0_at8_true 		`CIF1_f0_mata_at8'" _newline ///
				     "global CIF2_f0_at8_true 		`CIF2_f0_mata_at8'" _newline ///
					 "global CIF1_f0_at9_true 		`CIF1_f0_mata_at9'" _newline ///
				     "global CIF2_f0_at9_true 		`CIF2_f0_mata_at9'" _newline ///
					 "global CIF1_f0_at10_true	    `CIF1_f0_mata_at10'" _newline ///
				     "global CIF2_f0_at10_true 	    `CIF2_f0_mata_at10'" _newline ///
					 ///
					 "global CIF1_diff_at1_true 		`CIF1_diff_at1'" _newline ///
				     "global CIF2_diff_at1_true 		`CIF2_diff_at1'" _newline ///
					 "global CIF1_diff_at2_true 		`CIF1_diff_at2'" _newline ///
				     "global CIF2_diff_at2_true 		`CIF2_diff_at2'" _newline ///
					 "global CIF1_diff_at3_true 		`CIF1_diff_at3'" _newline ///
				     "global CIF2_diff_at3_true 		`CIF2_diff_at3'" _newline ///
					 "global CIF1_diff_at4_true 		`CIF1_diff_at4'" _newline ///
				     "global CIF2_diff_at4_true 		`CIF2_diff_at4'" _newline ///
					 "global CIF1_diff_at5_true 		`CIF1_diff_at5'" _newline ///
				     "global CIF2_diff_at5_true 		`CIF2_diff_at5'" _newline ///
					 "global CIF1_diff_at6_true 		`CIF1_diff_at6'" _newline ///
				     "global CIF2_diff_at6_true 		`CIF2_diff_at6'" _newline ///
					 "global CIF1_diff_at7_true 		`CIF1_diff_at7'" _newline ///
				     "global CIF2_diff_at7_true 		`CIF2_diff_at7'" _newline ///
					 "global CIF1_diff_at8_true 		`CIF1_diff_at8'" _newline ///
				     "global CIF2_diff_at8_true 		`CIF2_diff_at8'" _newline ///
					 "global CIF1_diff_at9_true 		`CIF1_diff_at9'" _newline ///
				     "global CIF2_diff_at9_true 		`CIF2_diff_at9'" _newline ///
					 "global CIF1_diff_at10_true	    `CIF1_diff_at10'" _newline ///
				     "global CIF2_diff_at10_true 	    `CIF2_diff_at10'" _newline ///
					///
					 "global CIF1_ratio_at1_true 		`CIF1_ratio_at1'" _newline ///
				     "global CIF2_ratio_at1_true 		`CIF2_ratio_at1'" _newline ///
					 "global CIF1_ratio_at2_true 		`CIF1_ratio_at2'" _newline ///
				     "global CIF2_ratio_at2_true 		`CIF2_ratio_at2'" _newline ///
					 "global CIF1_ratio_at3_true 		`CIF1_ratio_at3'" _newline ///
				     "global CIF2_ratio_at3_true 		`CIF2_ratio_at3'" _newline ///
					 "global CIF1_ratio_at4_true 		`CIF1_ratio_at4'" _newline ///
				     "global CIF2_ratio_at4_true 		`CIF2_ratio_at4'" _newline ///
					 "global CIF1_ratio_at5_true 		`CIF1_ratio_at5'" _newline ///
				     "global CIF2_ratio_at5_true 		`CIF2_ratio_at5'" _newline ///
					 "global CIF1_ratio_at6_true 		`CIF1_ratio_at6'" _newline ///
				     "global CIF2_ratio_at6_true 		`CIF2_ratio_at6'" _newline ///
					 "global CIF1_ratio_at7_true 		`CIF1_ratio_at7'" _newline ///
				     "global CIF2_ratio_at7_true 		`CIF2_ratio_at7'" _newline ///
					 "global CIF1_ratio_at8_true 		`CIF1_ratio_at8'" _newline ///
				     "global CIF2_ratio_at8_true 		`CIF2_ratio_at8'" _newline ///
					 "global CIF1_ratio_at9_true 		`CIF1_ratio_at9'" _newline ///
				     "global CIF2_ratio_at9_true 		`CIF2_ratio_at9'" _newline ///
					 "global CIF1_ratio_at10_true	    `CIF1_ratio_at10'" _newline ///
				     "global CIF2_ratio_at10_true 	    `CIF2_ratio_at10'" _newline ///
					 ///
					 "global CIF_total_at1_true 		`CIF_total_at1'" _newline ///
				     "global CIF_total_at2_true 		`CIF_total_at2'" _newline ///
					 "global CIF_total_at3_true 		`CIF_total_at3'" _newline ///
				     "global CIF_total_at4_true 		`CIF_total_at4'" _newline ///
					 "global CIF_total_at5_true 		`CIF_total_at5'" _newline ///
				     "global CIF_total_at6_true 		`CIF_total_at6'" _newline ///
					 "global CIF_total_at7_true 		`CIF_total_at7'" _newline ///
				     "global CIF_total_at8_true 		`CIF_total_at8'" _newline ///
					 "global CIF_total_at9_true	        `CIF_total_at9'" _newline ///
				     "global CIF_total_at10_true 	    `CIF_total_at10'" _newline ///
					///
				     "global hr_gender_at1_true  `hr_gender_at1'" _newline ///
				     "global hr_gender_at2_true  `hr_gender_at2'" _newline ///
				     "global hr_gender_at3_true  `hr_gender_at3'" _newline ///
				     "global hr_gender_at4_true  `hr_gender_at4'" _newline ///
				     "global hr_gender_at5_true  `hr_gender_at5'" _newline ///
				     "global hr_gender_at6_true  `hr_gender_at6'" _newline ///
				     "global hr_gender_at7_true  `hr_gender_at7'" _newline ///
				     "global hr_gender_at8_true  `hr_gender_at8'" _newline ///
				     "global hr_gender_at9_true  `hr_gender_at9'" _newline ///
				     "global hr_gender_at10_true `hr_gender_at10'" _newline ///
					 ///
                     "global bsex_c1   			 `bsex_c1'" _newline ///
                     "global bsex_c2   			 `bsex_c2'" _newline ///
					 "global bsex_c2_tvc_1   	 `bsex_c2_tvc_1'" _newline ///
					 "global bsex_c2_tvc_2       `bsex_c2_tvc_2'" _newline ///
					 "global prop_sex_other   	 `prop_sex_other'" _newline ///
					 "global base_function   	 `base_function'" _newline ///
					 "global independence   	 `independence'" _newline ///
					 "global obs                 `obs'" _newline ///
					 "global sdage               `sdage'"         _newline 

	file close truth

	}
   }
	
	
