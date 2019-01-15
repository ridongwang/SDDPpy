#!/bin/bash

scp dduque@crunch.osl.northwestern.edu:/home/dduque/dduque_projects/SDDP/HydroExamples/Output/DisceteWassersteinSingleCut/Hydro_R10_AR1_T12_I10000* ./

## python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DisceteWassersteinSingleCut/ --exp_file=Hydro_R10_AR1_T12_I500_N3ESS --N=3 --plot_type=OOS
## python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DisceteWassersteinSingleCut/ --exp_file=Hydro_R10_AR1_T12_I500_N5ESS --N=5 --plot_type=OOS
python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DisceteWassersteinSingleCut/ --exp_file=Hydro_R10_AR1_T12_I10000_N10ESS --N=10 --plot_type=OOS
## python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DisceteWassersteinSingleCut/ --exp_file=Hydro_R10_AR1_T12_I500_N20ESS --N=20 --plot_type=OOS
## python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DisceteWassersteinSingleCut/ --exp_file=Hydro_R10_AR1_T12_I500_N30ESS --N=30 --plot_type=OOS
## python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DisceteWassersteinSingleCut/ --exp_file=Hydro_R10_AR1_T12_I500_N50ESS --N=50 --plot_type=OOS
## python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DisceteWassersteinSingleCut/ --exp_file=Hydro_R10_AR1_T12_I500_N100ESS --N=100 --plot_type=OOS
