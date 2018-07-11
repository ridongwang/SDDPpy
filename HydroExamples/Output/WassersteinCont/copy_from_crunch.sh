#!/bin/bash

scp dduque@crunch.osl.northwestern.edu:/home/dduque/dduque_projects/SDDP/HydroExamples/Output/WassersteinCont/Hydro_R10_AR1_T12_I300* ./

python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/WassersteinCont/ --exp_file=Hydro_R10_AR1_T12_I300_N3ESS --N=3 --plot_type=OOS
python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/WassersteinCont/ --exp_file=Hydro_R10_AR1_T12_I300_N5ESS --N=5 --plot_type=OOS
python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/WassersteinCont/ --exp_file=Hydro_R10_AR1_T12_I300_N10ESS --N=10 --plot_type=OOS
python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/WassersteinCont/ --exp_file=Hydro_R10_AR1_T12_I300_N20ESS --N=20 --plot_type=OOS
python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/WassersteinCont/ --exp_file=Hydro_R10_AR1_T12_I300_N30ESS --N=30 --plot_type=OOS
python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/WassersteinCont/ --exp_file=Hydro_R10_AR1_T12_I300_N50ESS --N=50 --plot_type=OOS
python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/WassersteinCont/ --exp_file=Hydro_R10_AR1_T12_I300_N100ESS --N=100 --plot_type=OOS
