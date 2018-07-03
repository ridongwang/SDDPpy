#!/bin/bash

scp dduque@crunch.osl.northwestern.edu:/home/dduque/dduque_projects/SDDP/HydroExamples/Output/WassersteinCont/Hydro_* ./


python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/WassersteinCont/ --exp_file=Hydro_R10_AR1_T12_I100_N5ESS --N=5
python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/WassersteinCont/ --exp_file=Hydro_R10_AR1_T12_I100_N10ESS --N=10
python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/WassersteinCont/ --exp_file=Hydro_R10_AR1_T12_I100_N30ESS --N=30
