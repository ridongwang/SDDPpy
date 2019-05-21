NAME Hydrovalley_0 37 Variables 18 Constraints
 XL reservoir_level[0]  balance[0]
 XL reservoir_level[2]  balance[1]
 XL outflow[0]  balance[2]
 XL outflow[1]  generationCtr
 XL outflow[2]  R4      
 XL generation  outflowCtr[0]
 XL thermal_gen  outflowCtr[1]
 XL thermal_gen_cost  outflowCtr[2]
 XL dispatch[0,0]  dispatchCtr[1]
 XL dispatch[1,1]  dispatchCtr[2]
 XL dispatch[2,1]  R13     
 XL oracle[0][0]  cut[0,3,0]
 UL dispatch[2,2]
ENDATA
