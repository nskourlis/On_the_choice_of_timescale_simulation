/*********************************/
/* Setting the scenario parameters*/
/**********************************/



/********************************************************************/
/* Parameters that remain constant throughout the different scenarios*/
/********************************************************************/

/*Lambdas and gammas for cause 1*/

//sample size
local obs= 2000

//mean and minimum age at diagnosis 
local meanage=70
local meanage_men=67
local meanage_women=73
local minage 18

/*Lambdas, gammas and pmix for cause 1- cancer mortality*/
local lambda1= exp(-0.3381)
local lambda2= exp(-4.8856)
local gamma1= exp(-0.2658)
local gamma2= exp(-0.2333)
local pmix=invlogit(-0.9365) 

/*Parameters for Quardratic effect of age at diagnosis on hazard for cause 1- cancer mortality*/

local b1=log(0.9969324)
local b2= log(1.000136)


/*Beta coefficient for the effect of gender on hazard for cause 1- cancer mortality*/
local bsex_c1=log(1)


/**Setting minimum and maximum for the simulated time since diagnosis survival times  */
local mintime_cancer 0.0027
local maxtime_cancer 10

/*Parameters for the polynomial hazard for other cause mortality logH(a) = (c1a^2+c2a+c3)*/

local c1= -0.001326
local c2= 0.331
local c3= -18.15


/***Initial hazard used for Gompertz/Makeham and Adapted Weibull hazard****/
local alpha .01
	 

/**********************************************************************************/
/* Loop for creating the parameters that change throughout the different scenarios*/
/**********************************************************************************/
local p=0


//o stands for the standard deviation of age at diagnosis
//i  stands for identifying the other cause mortality hazard
//j identifies the age at diagnosis- gender dependence or independence
//k identifies whether gender has proportional or non proportional hazards for other cause mortality 

foreach  o in  10 15 {
	
		forvalues i =1/3 {

			forvalues j =0/1 {
				  
				    forvalues k =1/2 {
	  
	local sdage=`o'
	 
	  /*******************************************************/
	  
	  tknz "Weibull,Polynomial,Gompertz", s(v) p(,) nochar 
	  local base_function="`v`i''"  
	  
	  local independence=`j' 
	  
	  /*Lambda and gamma for cause II- other cause mortalityfor Adapted weibull*/

	  if "`base_function'"== "`v1'" {
		local lambda=  exp(-55.26)
        local gamma=  exp(2.508) 
	  }
	  
	  /*Lambda and gamma for cause II- other cause mortalityfor Gompertz/Makeham*/

	  else  if "`base_function'"== "`v3'" {
	  local lambda= 0.00007 
      local gamma= 0.087
	  }
	   
	 /*Proportional or non proporional hazards for gender on other cause mortality on attained age timescale*/
	 
	  tknz "PH_sex_other,Non_PH_sex_other", s(t) p(,) nochar 
	  
	   local prop_sex_other="`t`k''" 
	  
	   if "`prop_sex_other'"=="`t1'" {
	   
			local bsex_c2=log(0.7)
	   
			local bsex_c2_tvc_1 = 0
	   
			local bsex_c2_tvc_2 = 0
	   
	   }
	   
	    else if "`prop_sex_other'"=="`t2'" {
	   
			local bsex_c2= -0.9966

			local bsex_c2_tvc_1 = 0.0023
	   
			local bsex_c2_tvc_2 =0.000076667

	   }
	  
	 local p= `p'+1
	   
	
	// write true values to folder 1b.Scenarios_folder
	file open truth using "$N/1b.Scenarios_folder/Scenario`p'_new.do", write replace
	
	file write truth "// Scenario`p'_new" _newline ////
					 "global obs    	        `obs'"     _newline ///
				     "global lambda1     	    `lambda1'"     _newline ///
					 "global gamma1      	    `gamma1'"      _newline ///
				     "global lambda2     	    `lambda2'"     _newline ///
					 "global gamma2      	    `gamma2'"      _newline ///
				     "global pmix        	    `pmix'"        _newline ///
                     "global lambda 	        `lambda'" _newline ///
					 "global gamma  	        `gamma'"  _newline ///
					 "global b1       		    `b1'"          _newline ///
				     "global b2     		    `b2'"          _newline ///
					 "global c1                 `c1'"          _newline ///
				     "global c2     		    `c2'"          _newline ///
					 "global c3       		    `c3'"         _newline ///
					 "global alpha              `alpha'"    _newline ///
					 "global bsex_c1       		`bsex_c1'"     _newline ///
				     "global bsex_c2      		`bsex_c2'"     _newline ///
					 "global bsex_c2_tvc_1      `bsex_c2_tvc_1'"  _newline ///
					 "global bsex_c2_tvc_2      `bsex_c2_tvc_2'"  _newline ///
					 "global minage       		`minage'"        _newline ///
					 "global mintime_cancer 	`mintime_cancer'"  _newline ///
					 "global maxtime_cancer 	`maxtime_cancer'"  _newline ///
                     "global independence 		`independence'" _newline ////
					 "global prop_sex_other   	`prop_sex_other'" _newline ///
					 "global base_function 		`base_function'" _newline ///
					 "global meanage 		    `meanage'" _newline ///
					 "global meanage_men 		`meanage_men'" _newline ///
					 "global meanage_women		`meanage_women'" _newline ///
					 "global sdage 		        `sdage'"            _newline 
					 

		 

	file close truth
	       }
	     }
       } 
     }
	 
/********************************************************************************/
	 

