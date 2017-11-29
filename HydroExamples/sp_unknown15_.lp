\ Model Hydrovalley
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  1000 spill[0] + 1000 spill[1] + 1000 spill[2] + 1000 spill[3]
   + 1000 pour[0] + 1000 pour[1] + 1000 pour[2] + 1000 pour[3]
   - 5 generation + oracle[15]
Subject To
 AR1[0]: inflow[0] - 0.894155785206891 inflow0[0]
   - 0.187191166096471 inflow0[1] + 0.128736587791406 inflow0[2]
   - 0.128227942329335 inflow0[3] - innovations[0] = 0
 AR1[1]: inflow[1] + 0.328397404318301 inflow0[0]
   - 0.6951839351031071 inflow0[1] + 0.0913626214178796 inflow0[2]
   - 0.0289323912343019 inflow0[3] - innovations[1] = 0
 AR1[2]: inflow[2] - 0.585025574269405 inflow0[0]
   - 0.531110268093871 inflow0[1] + 0.113118385105778 inflow0[2]
   + 0.16903444953985 inflow0[3] - innovations[2] = 0
 AR1[3]: inflow[3] - 0.510092039498543 inflow0[0]
   - 0.203814073453963 inflow0[1] - 0.0196469487416793 inflow0[2]
   - 0.150465728201917 inflow0[3] - innovations[3] = 0
 balance[0]: reservoir_level[0] - reservoir_level0[0] - inflow[0]
   + outflow[0] + spill[0] - pour[0] = 0
 balance[1]: reservoir_level[1] - reservoir_level0[1] - inflow[1]
   - outflow[0] + outflow[1] - spill[0] + spill[1] - pour[1] = 0
 balance[2]: reservoir_level[2] - reservoir_level0[2] - inflow[2]
   - outflow[1] + outflow[2] - spill[1] + spill[2] - pour[2] = 0
 balance[3]: reservoir_level[3] - reservoir_level0[3] - inflow[3]
   - outflow[2] + outflow[3] - spill[2] + spill[3] - pour[3] = 0
 generationCtr: generation - 55 dispatch[0,0] - 65 dispatch[0,1]
   - 70 dispatch[0,2] - 55 dispatch[1,0] - 65 dispatch[1,1]
   - 70 dispatch[1,2] - 55 dispatch[2,0] - 65 dispatch[2,1]
   - 70 dispatch[2,2] - 55 dispatch[3,0] - 65 dispatch[3,1]
   - 70 dispatch[3,2] = 0
 outflowCtr[0]: outflow[0] - 50 dispatch[0,0] - 60 dispatch[0,1]
   - 70 dispatch[0,2] = 0
 outflowCtr[1]: outflow[1] - 50 dispatch[1,0] - 60 dispatch[1,1]
   - 70 dispatch[1,2] = 0
 outflowCtr[2]: outflow[2] - 50 dispatch[2,0] - 60 dispatch[2,1]
   - 70 dispatch[2,2] = 0
 outflowCtr[3]: outflow[3] - 50 dispatch[3,0] - 60 dispatch[3,1]
   - 70 dispatch[3,2] = 0
 dispatchCtr[0]: dispatch[0,0] + dispatch[0,1] + dispatch[0,2] <= 1
 dispatchCtr[1]: dispatch[1,0] + dispatch[1,1] + dispatch[1,2] <= 1
 dispatchCtr[2]: dispatch[2,0] + dispatch[2,1] + dispatch[2,2] <= 1
 dispatchCtr[3]: dispatch[3,0] + dispatch[3,1] + dispatch[3,2] <= 1
Bounds
 reservoir_level[0] <= 500
 reservoir_level[1] <= 500
 reservoir_level[2] <= 500
 reservoir_level[3] <= 500
 reservoir_level0[0] = 284.3896787739251
 reservoir_level0[1] = 500
 reservoir_level0[2] = 500
 reservoir_level0[3] = 500
 inflow0[0] = 49.44424968855506
 inflow0[1] = 7.44736233297882
 inflow0[2] = 22.01852147052132
 inflow0[3] = 39.5565816885171
 innovations[0] = -12.626063
 innovations[1] = 10.734489
 innovations[2] = -11.599395
 innovations[3] = -8.867397
 dispatch[0,0] <= 1
 dispatch[0,1] <= 1
 dispatch[0,2] <= 1
 dispatch[1,0] <= 1
 dispatch[1,1] <= 1
 dispatch[1,2] <= 1
 dispatch[2,0] <= 1
 dispatch[2,1] <= 1
 dispatch[2,2] <= 1
 dispatch[3,0] <= 1
 dispatch[3,1] <= 1
 dispatch[3,2] <= 1
 oracle[15] >= -100000
End
