#!/bin/bash

scp dduque@crunch.osl.northwestern.edu:/home/dduque/dduque_projects/SDDP/HydroExamples/Output/DisceteWassersteinSingleCut/Hydro_R10_AR1_T12_N* ./DisceteWassersteinSingleCut/
scp dduque@crunch.osl.northwestern.edu:/home/dduque/dduque_projects/SDDP/HydroExamples/Output/DW_DUal/Hydro_R10_AR1_T12_N* ./DW_Dual/

python /Users/dduque/Dropbox/WORKSPACE/SDDP/OutputAnalysis/SimulationAnalysis.py 
-