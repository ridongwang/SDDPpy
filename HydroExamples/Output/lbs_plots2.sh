#!/bin/bash

scp dduque@crunch.osl.northwestern.edu:/home/dduque/dduque_projects/SDDP/HydroExamples/Output/DisceteWassersteinSingleCut/Hydro_R10_AR1_T12_N* ./DisceteWassersteinSingleCut/
scp dduque@crunch.osl.northwestern.edu:/home/dduque/dduque_projects/SDDP/HydroExamples/Output/DW_Dual/Hydro_R10_AR1_T12_N* ./DW_Dual/

rootdirSC=DisceteWassersteinSingleCut/
rootdirDual=DW_Dual/
rootname=Hydro_R10_AR1_T12_N
midnameSC=_I2000ESS_SC_DW_
midnameDual=_I1000ESS_MC_DW_
Hydro_R10_AR1_T12_N5_I1000ESS_MC_DW_0.500000_ES
endnameDS=_DS_LBS.pickle
endnameES=_ES_LBS.pickle
DSONLY=DynamicSamplingOnly
for n in {5,10,20}
do
	for r in {0.500000,1.000000,5.000000,10.000000,15.000000,20.000000,35.000000,50.000000}
	do
		echo $n $r
		python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py $rootdirSC$rootname$n$midnameSC$r$endnameDS $rootdirSC$rootname$n$midnameSC$r$endnameES $rootdirDual$rootname$n$midnameDual$r$endnameDS $rootdirDual$rootname$n$midnameDual$r$endnameES --exp_file=$rootname$n --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/ --plot_type=LBS
		python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py $rootdirSC$rootname$n$midnameSC$r$endnameDS $rootdirDual$rootname$n$midnameDual$r$endnameDS --exp_file=$rootname$n$DSONLY --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/ --plot_type=LBS
	done
done


### python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py DisceteWassersteinSingleCut/Hydro_R10_AR1_T12_N5_I500ESS_Primal_MC_DW_10.000000_DS_LBS.pickle DisceteWassersteinSingleCut/Hydro_R10_AR1_T12_N5_I500ESS_Primal_MC_DW_10.000000_ES_LBS.pickle DisceteWassersteinSingleCut/Hydro_R10_AR1_T12_N5_I500ESS_Primal_SC_DW_10.000000_DS_LBS.pickle DisceteWassersteinSingleCut/Hydro_R10_AR1_T12_N5_I500ESS_Primal_SC_DW_10.000000_ES_LBS.pickle DW_Dual/Hydro_R10_AR1_T12_N5_I500ESS_Dual_MC_DW_10.000000_DS_LBS.pickle DW_Dual/Hydro_R10_AR1_T12_N5_I500ESS_Dual_MC_DW_10.000000_ES_LBS.pickle --exp_file=PrimalLBS --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/ --plot_type=LBS