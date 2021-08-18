/* This do file generates the simulated data for tvc */

	cd "$N\"
	
	global interval=10
	// loop over scenarios
    foreach scen in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 {
	
	// load scenario parameters
	include "$N\1b.Scenarios_folder\Scenario`scen'_new.do"
	
	set seed 642764
	
	// loop over datasets to generate the surviva times, number of simulated datasets =1000
	
	forvalues i = 1/1000 {
	
	clear
	set obs 2000
				
	// Generate age at diagnosis and gender for each individual
	
	//If age at diagnosis- gender are independent 
	
	if ${independence}==1 {
	
		gen agediag=rnormal(${meanage},${sdage})

		replace agediag=${minage} if  agediag<${minage}

		replace agediag=105 if  agediag>=105

		gen female=rbinomial(1,0.50 )
}

	//If age at diagnosis- gender are dependent 

	else if ${independence}==0 {
	
		gen female=rbinomial(1,0.50 )

		gen agediag=.

		replace agediag=rnormal(${meanage_men},${sdage}) if female==0

		replace agediag=rnormal(${meanage_women},${sdage}) if female==1		

		replace agediag=${minage} if  agediag<${minage} 

		replace agediag=105 if  agediag>=105
	}
	
summ agediag
gen agediagc=agediag-65
gen agediagc2=agediagc^2

	
			
	//Generate survival times for death due to cancer

	survsim t_cancer d_cancer , dist(weib) mixture lambda(${lambda1} ${lambda2}) gamma(${gamma1} ${gamma2}) pmix(${pmix}) ///
								maxt(10) cov(female ${bsex_c1} agediagc ${b1} agediagc2 ${b2})
								
	replace t_cancer = 0.0027 if t_cancer<0.0027	

	//Generate survival times for other cause mortality
	
	// For adapted weibull basline hazard for other cause mortality
		 if "$base_function"=="Weibull"  {	

	
			survsim t_other d_other, hazard(  ///
	( (${alpha}:+  ${lambda}:*${gamma}:*(agediag:+#t):^(${gamma}:-1)):* ///
			 exp(${bsex_c2}:*female :+ ${bsex_c2_tvc_1}:* female:*(agediag:+#t):+ ${bsex_c2_tvc_2}:*female:*(agediag:+#t):^2  )  ) ///
			)  maxt(10) nodes(200)  centol(1e-12)
		qui replace t_other = 0.0027 if t_other<0.0027	
	
	}

	// For Polynomial baseline hazard for other cause mortality
	
	else if  "$base_function"=="Polynomial" {
			
				survsim t_other d_other , hazard( ///
	(  ( ( 2:*${c1}:*(agediag:+#t)):+ ${c2} ):*  exp(  ${c1}:*  ((agediag:+#t):^2)    :+   ${c2}:* (agediag:+#t)  :+  ${c3} ) ) :*  ///
	 exp(${bsex_c2}:*female :+ ${bsex_c2_tvc_1}:*female:*(agediag:+#t) :+ ${bsex_c2_tvc_2}:*female:*(agediag:+#t):^2  ) )  ///
			mintime(${mintime_cancer}) maxtime(${maxtime_cancer}) 	
							
				}
	
	// For Gompertz/Makeham baseline hazard for other cause mortality

	else if  "$base_function"=="Gompertz" {
		
		survsim t_other d_other, hazard( ///
		( (${alpha} :+ ${lambda}:*exp(${gamma}:*(agediag:+#t)))    :* ///
		exp(${bsex_c2}:*female :+ ${bsex_c2_tvc_1}:*female:*(agediag:+#t):+ ${bsex_c2_tvc_2}:*female:*(agediag:+#t):^2  )   ) 	///
						)	 maxt(10) nodes(200)  centol(1e-12)
		
		qui replace t_other = 0.0027 if t_other<0.0027

			}


// Generate whether the event happened and the observed event time and status for each individual

	gen double time_s_diag = min(t_cancer,t_other)
	gen died=cond(time_s_diag == ${maxtime_cancer},0,cond(t_cancer<t_other,1,2))

	
    gen cancerd=0
    replace cancerd=1 if t_cancer<=t_other & d_cancer==1 

    gen otherd=0
    replace otherd=1 if  t_other<t_cancer  & d_other==1 

    gen status=0
    replace status=1 if cancerd==1
    replace status=2 if otherd==1
	
	sum agediagc,  det

	tab cancerd otherd
	tab status
	
   // save data to scenario specific folder
  quietly save "$N\3b.Simulated_data_folder\Scenario`scen'\Sim_data`i'.dta", replace
 }
  
}
