//Do file that appends the performance measures result files from each model and scenario

/****************************************************************************************************/

// Create a final dataset with the performance measure results (bias, relative bias, convergence, */
// coverage, Monte carlo errors, empirical standard errors, relative efficiency) over gender, different
// ages at diagnosis, time since diagnosis and scenarios 

/****************************************************************************************************/
clear 

//Appending performance measure datasets over the scenarios 

foreach scen in   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24  {
					  
use  "$N\5b.Performance_results\compare_attained\scen_`scen'_compare_attained.dta", replace
save "$N\5b.Performance_results\data_merge/scen_`scen'_data0.dta", replace

use  "$N\5b.Performance_results\compare_linear\scen_`scen'_compare_linear.dta", replace
save "$N\5b.Performance_results\data_merge/scen_`scen'_data1.dta", replace

use "$N\5b.Performance_results\compare_splines\scen_`scen'_compare_splines.dta", replace
save "$N\5b.Performance_results\data_merge/scen_`scen'_data2.dta", replace

use "$N\5b.Performance_results\compare_splines_int\scen_`scen'_compare_splines_int.dta", replace
save "$N\5b.Performance_results\data_merge/scen_`scen'_data3.dta", replace

use "$N\5b.Performance_results\data_merge/scen_`scen'_data0.dta", clear

forval k = 1/3 {   
 
     append using "$N\5b.Performance_results\data_merge/scen_`scen'_data`k'.dta"

	}	
	 
save "$N\5b.Performance_results\data_merge/scen_`scen'_final.dta", replace

 }
 
 use "$N\5b.Performance_results\data_merge/scen_1_final.dta", clear
 
foreach scen in   2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 {

		append using "$N\5b.Performance_results\data_merge/scen_`scen'_final.dta"

		}
 
 save "$N\5b.Performance_results\data_merge/scen_final.dta", replace

/********************************************************************************************/

//Erase the files we do not need any more (files created in the above procedure)*/

 foreach scen in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24  {
      
	      erase  "$N\5b.Performance_results\data_merge/scen_`scen'_final.dta"
	  
	  forval k = 0/3 { 
	  
	      erase  "$N\5b.Performance_results\data_merge/scen_`scen'_data`k'.dta"
 
   }

}

drop female*
drop CIF_total_at*

//preparing data for Rshiny and the bias reactive graph - They will be reshaped to fully long format
gen id=_n

		forval female=0(1)1 {
		global female=`female'

 drop EmpSE_CIF1_f${female}_at* EmpSE_CIF2_f${female}_at*  

foreach k in CIF1 CIF2  {
		foreach s in 1 2 3 4 5 6 7 8 9 10 {
			  
			  rename `k'_f${female}_at`s'_t1t1  `k'_f${female}_at`s'
			  drop `k'_f${female}_at`s'_t1t2
		
	          rename bias_`k'_f${female}_at`s'_t1t1  bias_`k'_f${female}_at`s'
			  drop bias_`k'_f${female}_at`s'_t1t2
			  
			  rename rel_bias_`k'_f${female}_at`s'_t1t1  rel_bias_`k'_f${female}_at`s'
			  drop rel_bias_`k'_f${female}_at`s'_t1t2
			  
			  rename MCer_`k'_f${female}_at`s'_t1t1  MCer_`k'_f${female}_at`s'
			  drop MCer_`k'_f${female}_at`s'_t1t2
			  
			  rename rb_MCer_`k'_f${female}_at`s'_t1t1  rb_MCer_`k'_f${female}_at`s'
			  drop rb_MCer_`k'_f${female}_at`s'_t1t2
			  
			  rename Coverage_`k'_f${female}_at`s'_t1t1  Coverage_`k'_f${female}_at`s'
			  drop Coverage_`k'_f${female}_at`s'_t1t2

			  }
			  }


		}

		foreach s in 1 2 3 4 5 6 7 8 9 10 {
	          rename hr_gender_at`s'_t1t1  hr_gender_at`s'
			  replace hr_gender_at`s'= hr_gender_at`s'_t1t2 if model==0
			  drop hr_gender_at`s'_t1t2
			 
			 rename bias_hr_gender_at`s'_t1t1  bias_hr_gender_at`s'
			 replace bias_hr_gender_at`s'= bias_hr_gender_at`s'_t1t2 if model==0
			  drop  bias_hr_gender_at`s'_t1t2
			  }
			   

foreach k in CIF1 CIF2  {
		foreach s in 1 2 3 4 5 6 7 8 9 10 {
	          rename `k'_diff_at`s'_t1t1  `k'_diff_at`s'
			  replace `k'_diff_at`s'= `k'_diff_at`s'_t1t2 if model==0
			  drop `k'_diff_at`s'_t1t2
			  
			  rename `k'_ratio_at`s'_t1t1  `k'_ratio_at`s'
			  replace `k'_ratio_at`s'= `k'_ratio_at`s'_t1t2 if model==0
			  drop `k'_ratio_at`s'_t1t2
			  
			  }
			  }
			  
foreach t in CIF1 CIF2  {
		foreach s in 1 2 3 4 5 6 7 8 9 10 {
	          rename `t'_diff_at`s'_true `t'_diff_true_at`s'
	          rename `t'_ratio_at`s'_true `t'_ratio_true_at`s'
			  }
			  }			  
		foreach s in 1 2 3 4 5 6 7 8 9 10 {
	          rename hr_gender_at`s'_true hr_gender_true_at`s'
			  }
			  
			  foreach t in CIF1 CIF2  {
			foreach f in 0 1 { 
		foreach s in 1 2 3 4 5 6 7 8 9 10 {

	               rename `t'_f`f'_at`s'_true 	`t'_f`f'_true_at`s'	  
				   
				   }
}
}
	//  drop *true
	 
gen model_lab="Attained" if model==0
replace model_lab="Linear" if model==1
replace model_lab="Splines" if model==2
replace model_lab="Splines/Int" if model==3


 save"$N\5b.Performance_results\data_merge/scen_final.dta", replace

reshape long CIF1_f0_at CIF2_f0_at  ///
			 bias_CIF1_f0_at bias_CIF2_f0_at  ///
			 CIF1_f0_true_at CIF2_f0_true_at ///
			 MCer_CIF1_f0_at MCer_CIF2_f0_at ///
			 rb_MCer_CIF1_f0_at rb_MCer_CIF2_f0_at ///
			 rel_bias_CIF1_f0_at rel_bias_CIF2_f0_at  ///
			 Coverage_CIF1_f0_at Coverage_CIF2_f0_at ///
			 Rel_prec_CIF1_f0_at Rel_prec_CIF2_f0_at ///
			 hr_gender_at bias_hr_gender_at, i(id) j(year)
			 
	order givenage scenario model independence  prop_sex_other base_function model_lab 
	order id year CIF1_f0_at CIF2_f0_at bias_CIF1_f0_at bias_CIF2_f0_at CIF1_f0_true_at CIF2_f0_true_at MCer_CIF1_f0_at MCer_CIF2_f0_at rb_MCer_CIF1_f0_at rb_MCer_CIF2_f0_at  rel_bias_CIF1_f0_at rel_bias_CIF2_f0_at Rel_prec_CIF1_f0_at Rel_prec_CIF2_f0_at Coverage_CIF1_f0_at Coverage_CIF2_f0_at hr_gender_at bias_hr_gender_at	 
	keep  id year CIF1_f0_at CIF2_f0_at bias_CIF1_f0_at bias_CIF2_f0_at CIF1_f0_true_at CIF2_f0_true_at MCer_CIF1_f0_at MCer_CIF2_f0_at rb_MCer_CIF1_f0_at rb_MCer_CIF2_f0_at  rel_bias_CIF1_f0_at rel_bias_CIF2_f0_at Rel_prec_CIF1_f0_at Rel_prec_CIF2_f0_at Coverage_CIF1_f0_at Coverage_CIF2_f0_at hr_gender_at bias_hr_gender_at ///
		givenage scenario model independence  prop_sex_other base_function model_lab convergence_t1t1 convergence_t1t2 obs sdage
	
	rename CIF1_f0_at          CIF1_at 
	rename CIF2_f0_at          CIF2_at 	 
	rename bias_CIF1_f0_at     bias_CIF1_at 
	rename bias_CIF2_f0_at     bias_CIF2_at 
	rename CIF1_f0_true_at     CIF1_true_at
	rename CIF2_f0_true_at     CIF2_true_at
	rename MCer_CIF1_f0_at     MCer_CIF1_at 
	rename MCer_CIF2_f0_at     MCer_CIF2_at 
	rename rb_MCer_CIF1_f0_at  rb_MCer_CIF1_at 
	rename rb_MCer_CIF2_f0_at  rb_MCer_CIF2_at 
	rename rel_bias_CIF1_f0_at rel_bias_CIF1_at
	rename rel_bias_CIF2_f0_at rel_bias_CIF2_at
	rename Coverage_CIF1_f0_at Coverage_CIF1_at
	rename Coverage_CIF2_f0_at Coverage_CIF2_at
	rename Rel_prec_CIF1_f0_at Rel_prec_CIF1_at
	rename Rel_prec_CIF2_f0_at Rel_prec_CIF2_at
		 
 save "$N/5b.Performance_results/data_merge/scen_final_reshape0_part1.dta", replace
  
  use "$N/5b.Performance_results/data_merge/scen_final.dta", replace

  reshape long 		CIF1_diff_at CIF2_diff_at CIF1_ratio_at CIF2_ratio_at ///
					CIF1_diff_true_at CIF2_diff_true_at CIF1_ratio_true_at CIF2_ratio_true_at hr_gender_true_at , i(id) j(year)
				
order givenage scenario model independence  prop_sex_other base_function model_lab
order id year CIF1_diff_at CIF2_diff_at CIF1_ratio_at CIF2_ratio_at CIF1_diff_true_at CIF2_diff_true_at CIF1_ratio_true_at CIF2_ratio_true_at hr_gender_true_at
keep id year CIF1_diff_at CIF2_diff_at CIF1_ratio_at CIF2_ratio_at CIF1_diff_true_at CIF2_diff_true_at CIF1_ratio_true_at CIF2_ratio_true_at  hr_gender_true_at ///
	 givenage scenario model independence  prop_sex_other base_function model_lab convergence_t1t1 convergence_t1t2 obs sdage
	 

	  merge 1:1 id year using "$N/5b.Performance_results/data_merge/scen_final_reshape0_part1.dta"
	drop _merge
		gen female=0
save "$N/5b.Performance_results/data_merge/scen_final_reshape0_final.dta", replace
  

  use "$N/5b.Performance_results/data_merge/scen_final.dta", replace

 reshape long CIF1_f1_at CIF2_f1_at  ///
			  bias_CIF1_f1_at bias_CIF2_f1_at  ///
			  CIF1_f1_true_at CIF2_f1_true_at ///
			  MCer_CIF1_f1_at MCer_CIF2_f1_at  /// 
			  rb_MCer_CIF1_f1_at rb_MCer_CIF2_f1_at  /// 
			  rel_bias_CIF1_f1_at rel_bias_CIF2_f1_at  ///
			  Coverage_CIF1_f1_at Coverage_CIF2_f1_at ///
			  Rel_prec_CIF1_f1_at Rel_prec_CIF2_f1_at ///
			  hr_gender_at bias_hr_gender_at, i(id) j(year)

	order givenage scenario model independence  prop_sex_other base_function model_lab 
	order id year CIF1_f1_at CIF2_f1_at bias_CIF1_f1_at bias_CIF2_f1_at CIF1_f1_true_at CIF2_f1_true_at MCer_CIF1_f1_at MCer_CIF2_f1_at rb_MCer_CIF1_f1_at rb_MCer_CIF2_f1_at rel_bias_CIF1_f1_at rel_bias_CIF2_f1_at Rel_prec_CIF1_f1_at Rel_prec_CIF2_f1_at Coverage_CIF1_f1_at Coverage_CIF2_f1_at hr_gender_at bias_hr_gender_at	 
	keep id year  CIF1_f1_at CIF2_f1_at bias_CIF1_f1_at bias_CIF2_f1_at CIF1_f1_true_at CIF2_f1_true_at MCer_CIF1_f1_at MCer_CIF2_f1_at rb_MCer_CIF1_f1_at rb_MCer_CIF2_f1_at rel_bias_CIF1_f1_at rel_bias_CIF2_f1_at Rel_prec_CIF1_f1_at Rel_prec_CIF2_f1_at Coverage_CIF1_f1_at Coverage_CIF2_f1_at hr_gender_at bias_hr_gender_at ///
		givenage scenario model independence  prop_sex_other base_function model_lab convergence_t1t1 convergence_t1t2 obs sdage
	
	rename CIF1_f1_at          CIF1_at 
	rename CIF2_f1_at          CIF2_at 	 
    rename bias_CIF1_f1_at     bias_CIF1_at 
	rename bias_CIF2_f1_at     bias_CIF2_at 
	rename CIF1_f1_true_at     CIF1_true_at
	rename CIF2_f1_true_at     CIF2_true_at
	rename MCer_CIF1_f1_at     MCer_CIF1_at 
	rename MCer_CIF2_f1_at     MCer_CIF2_at 
	rename rb_MCer_CIF1_f1_at  rb_MCer_CIF1_at 
	rename rb_MCer_CIF2_f1_at  rb_MCer_CIF2_at 
	rename rel_bias_CIF1_f1_at rel_bias_CIF1_at
	rename rel_bias_CIF2_f1_at rel_bias_CIF2_at
	rename Coverage_CIF1_f1_at Coverage_CIF1_at
	rename Coverage_CIF2_f1_at Coverage_CIF2_at
	rename Rel_prec_CIF1_f1_at Rel_prec_CIF1_at
	rename Rel_prec_CIF2_f1_at Rel_prec_CIF2_at
	
	 save "$N/5b.Performance_results/data_merge/scen_final_reshape1_part1.dta", replace
			 
    use "$N/5b.Performance_results/data_merge/scen_final.dta", replace

			 
	  reshape long 	CIF1_diff_at CIF2_diff_at CIF1_ratio_at CIF2_ratio_at ///
					CIF1_diff_true_at CIF2_diff_true_at CIF1_ratio_true_at CIF2_ratio_true_at hr_gender_true_at  , i(id) j(year)
				
order givenage scenario model independence  prop_sex_other base_function model_lab 
order id year CIF1_diff_at CIF2_diff_at CIF1_ratio_at CIF2_ratio_at CIF1_diff_true_at CIF2_diff_true_at CIF1_ratio_true_at CIF2_ratio_true_at hr_gender_true_at	
keep id year CIF1_diff_at CIF2_diff_at CIF1_ratio_at CIF2_ratio_at CIF1_diff_true_at CIF2_diff_true_at CIF1_ratio_true_at CIF2_ratio_true_at hr_gender_true_at ///
	 givenage scenario model independence  prop_sex_other base_function model_lab convergence_t1t1 convergence_t1t2 obs sdage
			 
	
merge 1:1 id year using "$N/5b.Performance_results/data_merge/scen_final_reshape1_part1.dta"
drop _merge

gen female=1
save "$N/5b.Performance_results/data_merge/scen_final_reshape1_final.dta", replace
  		 
sort female givenage scenario model		 

append using "$N/5b.Performance_results/data_merge/scen_final_reshape0_final.dta"
 
keep if scenario== 3 |  scenario== 4 | scenario== 7 | scenario== 8 | scenario== 11 | scenario== 12 | ///
        scenario== 15 | scenario== 16 | scenario== 19 | scenario== 20 | scenario== 23 | scenario== 24 
		
		recode scenario (3=1) (4=2) (7=3) (8=4) (11=5) (12=6) (15=7) (16=8) (19=9) (20=10) (23=11) (24=12) 
 
 		recode scenario (1=1) (3=2) (5=3) (7=4) (9=5) (11=6) (2=7) (4=8) (6=9) (8=10) (10=11) (12=12) 

		sort scenario
		
		li scenario prop_sex_other base_function sdage if givenage==60 & female==1 & model==0 & year==1
 
save "$N/5b.Performance_results/data_merge/scen_final_R_n.dta", replace
use "$N/5b.Performance_results/data_merge/scen_final_R_n.dta", replace

 
 preserve
 keep if female==1		
 keep if givenage==60 | givenage==70 | givenage==80
 sort  givenage scenario model year
 
bysort  givenage scenario model: egen mean_bias_CIF1_at= mean(bias_CIF1_at)
bysort  givenage scenario model: egen mean_bias_CIF2_at= mean(bias_CIF2_at)


sum rel_bias_CIF2
	
 save "$N/5b.Performance_results/data_merge/scen_final_R_n_long.dta", replace
 export delimited using "$N\5b.Performance_results/R_long.csv", replace
 restore
sort female  scenario givenage  model year
order female  scenario givenage  model year  hr_gender_at  hr_gender_true_at


		
export delimited using "$N\5b.Performance_results/Rtable_sim_agegender.csv", replace




li CIF2_at bias_CIF2_at CIF2_true_at rel_bias_CIF2_at givenage year model_lab base_function scenario if givenage==60 & year==5 & model==0 & female==1 & prop_sex_other=="PH_sex_other", abbreviate(30)
sum rel_bias_CIF2_at if  givenage==60 & year==5 & model==0 & female==1 & prop_sex_other=="PH_sex_other"

li CIF2_at bias_CIF2_at CIF2_true_at rel_bias_CIF2_at givenage year model_lab base_function scenario if givenage==70 & year==5 & model==3 & female==1 & prop_sex_other=="Non_PH_sex_other",abbreviate(30)
sum rel_bias_CIF2_at if  givenage==70 & year==5 & model==3 & female==1 & prop_sex_other=="Non_PH_sex_other"
