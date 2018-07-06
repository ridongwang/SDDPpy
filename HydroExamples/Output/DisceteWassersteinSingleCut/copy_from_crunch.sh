#!/bin/bash

## scp dduque@crunch.osl.northwestern.edu:/home/dduque/dduque_projects/SDDP/HydroExamples/Output/DisceteWassersteinSingleCut/Hydro_R10_AR1_T12_I* ./

python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py --path_to_files=/Users/dduque/Dropbox/WORKSPACE/SDDP/HydroExamples/Output/DisceteWassersteinSingleCut/ --exp_file=Hydro_R10_AR1_T12_I100_N10ESS --N=10 --plot_type=OOS