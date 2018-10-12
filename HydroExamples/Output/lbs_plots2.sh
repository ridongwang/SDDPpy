#!/bin/bash
# Copy file from crunch
## scp dduque@crunch.osl.northwestern.edu:/home/dduque/dduque_projects/SDDP/HydroExamples/Output/DW_Primal/Hydro_R10_AR1_T5_N* ./DW_Primal/
## scp dduque@crunch.osl.northwestern.edu:/home/dduque/dduque_projects/SDDP/HydroExamples/Output/DW_Dual/Hydro_R10_AR1_T5_N* ./DW_Dual/



rootdirPrimal=DW_Primal/
rootdirDual=DW_Dual/
rootname=Hydro_R10_AR1_T12_N
midnamePrimSC=_I100001ESS_PRIMAL_SC_DW_
midnamePrimMC=_I100001ESS_PRIMAL_MC_DW_
midnameDual=_I100001ESS_DUAL_MC_DW_
endnameDS=_DS_LBS.pickle
endnameES=_ES_LBS.pickle
DSONLY=DSOnly

# 0: OOS, 1: LBS 2:OSSGAP
plot_type=0



if [ $plot_type -eq 0 ]
then
iters=_I100001ESS
for n in {3,5,10,20}
do
	python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py  --exp_file=$rootname$n$iters --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DW_Primal/ --plot_type=OOS --sampling=DS --cut_type=MC --N=$n
	python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py  --exp_file=$rootname$n$iters --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DW_Primal/ --plot_type=OOS --sampling=ES --cut_type=MC --N=$n
	python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py  --exp_file=$rootname$n$iters --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DW_Primal/ --plot_type=OOS --sampling=DS --cut_type=SC --N=$n
	python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py  --exp_file=$rootname$n$iters --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DW_Primal/ --plot_type=OOS --sampling=ES --cut_type=SC --N=$n
	python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py  --exp_file=$rootname$n$iters --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DW_Dual/ --plot_type=OOS --sampling=DS --cut_type=MC --N=$n
	python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py  --exp_file=$rootname$n$iters --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DW_Dual/ --plot_type=OOS --sampling=ES --cut_type=MC --N=$n
done
fi

if [ $plot_type -eq 1 ]
then
for n in 10
do
	for r in 10.000000
	do
		### {0.100000,1.000000,5.000000,10.000000,15.000000,20.000000,25.000000}
		echo $n $r
		python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py $rootdirPrimal$rootname$n$midnamePrimSC$r$endnameDS $rootdirPrimal$rootname$n$midnamePrimMC$r$endnameDS $rootdirDual$rootname$n$midnameDual$r$endnameDS $rootdirPrimal$rootname$n$midnamePrimSC$r$endnameES $rootdirPrimal$rootname$n$midnamePrimMC$r$endnameES $rootdirDual$rootname$n$midnameDual$r$endnameES --exp_file=$rootname$n --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/ --plot_type=LBS
		python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py $rootdirPrimal$rootname$n$midnamePrimSC$r$endnameDS $rootdirPrimal$rootname$n$midnamePrimMC$r$endnameDS $rootdirDual$rootname$n$midnameDual$r$endnameDS --exp_file=$rootname$n$DSONLY --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/ --plot_type=LBS
	done
done
fi

if [ $plot_type -eq 2 ]
then
iters=_I100001ESS
for n in {5,30}
do
	echo $n
	python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py  --exp_file=$rootname$n$iters --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/ --plot_type=OOSGAP --N=$n
done
fi
