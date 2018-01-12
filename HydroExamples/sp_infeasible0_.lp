\ Model Hydrovalley
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  1000 pour[0] + 1000 pour[1] - generation + oracle[0]
Subject To
 AR_model_1[0]: inflow[0,1] - 0.541013086491916 inflow0[0,1]
   + 0.2277425926103031 inflow0[0,2] - innovations[0] = 0
 AR_model_1[1]: inflow[1,1] - 0.3446844960264348 inflow0[0,1]
   + 0.3279925410235554 inflow0[0,2] - 0.2000393020918058 inflow0[1,1]
   + 0.283230263357532 inflow0[1,2] - innovations[1] = 0
 AR_model[2,0]: inflow[0,2] - inflow0[0,1] = 0
 AR_model[2,1]: inflow[1,2] - inflow0[1,1] = 0
 balance[0]: reservoir_level[0] - reservoir_level0[0]
   - 0.541013086491916 inflow0[0,1] + 0.2277425926103031 inflow0[0,2]
   + outflow[0] + spill[0] - innovations[0] - pour[0] = 0
 balance[1]: reservoir_level[1] - reservoir_level0[1]
   - 0.3446844960264348 inflow0[0,1] + 0.3279925410235554 inflow0[0,2]
   - 0.2000393020918058 inflow0[1,1] + 0.283230263357532 inflow0[1,2]
   - outflow[0] + outflow[1] - spill[0] + spill[1] - innovations[1]
   - pour[1] = 0
 generationCtr: generation - 55 dispatch[0,0] - 65 dispatch[0,1]
   - 70 dispatch[0,2] - 55 dispatch[1,0] - 65 dispatch[1,1]
   - 70 dispatch[1,2] = 0
 outflowCtr[0]: outflow[0] - 50 dispatch[0,0] - 60 dispatch[0,1]
   - 70 dispatch[0,2] = 0
 outflowCtr[1]: outflow[1] - 50 dispatch[1,0] - 60 dispatch[1,1]
   - 70 dispatch[1,2] = 0
 dispatchCtr[0]: dispatch[0,0] + dispatch[0,1] + dispatch[0,2] <= 1
 dispatchCtr[1]: dispatch[1,0] + dispatch[1,1] + dispatch[1,2] <= 1
Bounds
 reservoir_level[0] <= 200
 reservoir_level[1] <= 200
 reservoir_level0[0] = 20
 reservoir_level0[1] = 20
 inflow[0,1] = 3.883370149642889
 inflow[0,2] = 15.82192072596131
 inflow[1,1] = 28.06993597040999
 inflow[1,2] = 20.82004094896339
 inflow0[0,1] = 15.82192072596131
 inflow0[0,2] = 6.422243032264402
 inflow0[1,1] = 20.82004094896339
 inflow0[1,2] = 29.91947275167659
 dispatch[0,0] <= 1
 dispatch[0,1] <= 1
 dispatch[0,2] <= 1
 dispatch[1,0] <= 1
 dispatch[1,1] <= 1
 dispatch[1,2] <= 1
 oracle[0] >= -1e+08
End
