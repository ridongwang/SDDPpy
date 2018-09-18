#!/bin/bash

scp dduque@crunch.osl.northwestern.edu:/home/dduque/dduque_projects/SDDP/HydroExamples/Output/DW_Primal/Hydro_R10_AR1_T12_N* ./DW_Primal/
scp dduque@crunch.osl.northwestern.edu:/home/dduque/dduque_projects/SDDP/HydroExamples/Output/DW_Dual/Hydro_R10_AR1_T12_N* ./DW_Dual/


### Hydro_R10_AR1_T12_N30_I100000ESS_DUAL_MC_DW_25.000000_ES_OOS.pickle
rootdirPrimal=DW_Primal/
rootdirDual=DW_Dual/
rootname=Hydro_R10_AR1_T12_N
midnamePrimSC=_I100000ESS_PRIMAL_SC_DW_
midnamePrimMC=_I100000ESS_PRIMAL_MC_DW_
midnameDual=_I100000ESS_DUAL_MC_DW_
### Hydro_R10_AR1_T12_N5_I1000ESS_MC_DW_0.500000_ES
endnameDS=_DS_LBS.pickle
endnameES=_ES_LBS.pickle
DSONLY=DSOnly
for n in {3,5,10,20,30}
do
	for r in {0.100000,1.000000,5.000000,10.000000,15.000000,20.000000,25.000000}
	do
		### {0.100000,1.000000,5.000000,10.000000,15.000000,20.000000,25.000000}
		echo $n $r
		python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py $rootdirPrimal$rootname$n$midnamePrimSC$r$endnameDS $rootdirPrimal$rootname$n$midnamePrimMC$r$endnameDS $rootdirDual$rootname$n$midnameDual$r$endnameDS $rootdirPrimal$rootname$n$midnamePrimSC$r$endnameES $rootdirPrimal$rootname$n$midnamePrimMC$r$endnameES $rootdirDual$rootname$n$midnameDual$r$endnameES --exp_file=$rootname$n --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/ --plot_type=LBS
		python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py $rootdirPrimal$rootname$n$midnamePrimSC$r$endnameDS $rootdirPrimal$rootname$n$midnamePrimMC$r$endnameDS $rootdirDual$rootname$n$midnameDual$r$endnameDS --exp_file=$rootname$n$DSONLY --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/ --plot_type=LBS
	done
done

###python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py DisceteWassersteinSingleCut/Hydro_R10_AR1_T12_N5_I500ESS_Primal_MC_DW_10.000000_DS_LBS.pickle DisceteWassersteinSingleCut/Hydro_R10_AR1_T12_N5_I500ESS_Primal_MC_DW_10.000000_ES_LBS.pickle DisceteWassersteinSingleCut/Hydro_R10_AR1_T12_N5_I500ESS_Primal_SC_DW_10.000000_DS_LBS.pickle DisceteWassersteinSingleCut/Hydro_R10_AR1_T12_N5_I500ESS_Primal_SC_DW_10.000000_ES_LBS.pickle DW_Dual/Hydro_R10_AR1_T12_N5_I500ESS_Dual_MC_DW_10.000000_DS_LBS.pickle DW_Dual/Hydro_R10_AR1_T12_N5_I500ESS_Dual_MC_DW_10.000000_ES_LBS.pickle --exp_file=PrimalLBS --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/ --plot_type=LBS
