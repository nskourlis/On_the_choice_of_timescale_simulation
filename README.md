# On_the_choice_of_timescale_simulation
Simulation code for manuscript "On the choice of timescale for other cause mortality in a competing risk setting using flexible parametric survival models"

Guide:
The datasets with the results of the analysis on the simulated datasets are uploaded. However only the 1st out of the 1000 simulated dataset is provided per scenario due to space limitations. The same for the true values of the scenarios. Due to space limitations only the true values for ages 70, 80 and 90 are provided. The user can easily recreate all the simulated datasets and true values of the CIFs by running the simulation code in their local PC.

Example: Template for project readme.txt file
Readme.txt
Name: Skourlis Nikolaos

--------------------------------------------------------------
This document provides explanations regarding the folders included in repository "On_the_choice_of_timescale_simulations"
--------------------------------------------------------------
FOLDERS:
0.Set_pathway: The user has to specify their working directory as a global variable
1a.Scenarios_generate: Generation of 24 scenarios (only the 12 scenarios are used in the manuscript). The scenarios are defined by the fully factorial combination of 4 factors. Factor I is standard deviation of age at diagnosis (sd=10 or sd=15). Factor II is baseline hazard shape for other cause mortality. Factor III is age at diagnosis- gender dependence or independence and Factor IV is  proportional or non proportional hazards of gender on other cause mortality on the attaned age timescale.
1b.Scenarios_folder: Folder with the parameters that are stored for each scenario.
2a.Generate_truths: Do files that are using the parameters stored in 1b.Scenarios_folder and through the appropriate equations the true values of CIFs are derived for each scenario.
2b.Truths_folder: The folder where the true CIFs values for each scenario are stored.
3a.Generate_simulated_data: Contains a do file of the same name which generates 1000 datasets for each scenario.
3b.Simulated_data_folder: The folder where the simulated datasets for each scenario are stored.
4a.Analyze_simulated_data: This folder has 4 do files (1.Analyze_total_attained,2.Analyse_delta_t1t1_linear,3.Analyse_delta_t1t1_splines,4.Analyse_delta_t1t1_splines_int) one do file for the different timescale approach a- Attained age and 3 do files for the common timescale approaches). Through this do files, we analyze the simulated datasets and estimate the CIF for cancer mortality and other cause mortality for males and females for ages at diagnosis 50,60,70,80,90.
4b.Analyze_results: This folder contains the results of the estimation analyses of the 4 different approaches over the different scenarios (4 subfolders).
5a.Performance_do: This folder has 4 do files, one for each of the 4 modelling approaches. In each do file the true and the estimated CIF values of each approach (over scenarios, gender and ages at diagnosis) are called. Then the bias, coverage, relative precision (of each common timescale approach compared to approach a-Attained age), relative bias and Monte carlo errors are computed. Then, the do file "Results_aggregation" aggregated the aforementioned performance measures that are stored in folder "5b.Performance_results" into a common dataset and exports it as a csv/excel file.
5b.Performance_results: The folder that contains/stores the performance measures of the simulation over scenarios, gender and ages at diagnosis for each different modelling approach (4 sub folders). A 5th sub folder "data_merge" contains the intermediate and final results of the aggregation of the perforamce measures.
