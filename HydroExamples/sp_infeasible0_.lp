\ Model Hydrovalley
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  10000 pour[0] + 10000 pour[1] + 10000 pour[2] + 10000 pour[3]
   + 10000 pour[4] + 10000 pour[5] + 10000 pour[6] + 10000 pour[7]
   + 10000 pour[8] + 10000 pour[9] - 10 generation + 0.2 oracle[0][0]
   + 0.2 oracle[0][1] + 0.2 oracle[0][2] + 0.2 oracle[0][3]
   + 0.2 oracle[0][4] + 0.001 lambda[0]
Subject To
 AR_model_1[0]: inflow[0,1] - 1.058431407580598 inflow0[0,1]
   - innovations[0] = 0
 AR_model_1[1]: inflow[1,1] - 0.5872427904634435 inflow0[0,1]
   - 0.240094325020334 inflow0[1,1] - innovations[1] = 0
 AR_model_1[2]: inflow[2,1] + 0.5863667279258465 inflow0[0,1]
   - 1.12053479408998 inflow0[1,1] - 1.344535919520875 inflow0[2,1]
   - innovations[2] = 0
 AR_model_1[3]: inflow[3,1] + 0.0619313110761347 inflow0[1,1]
   + 0.0908143249786187 inflow0[2,1] - 0.5700530505153536 inflow0[3,1]
   - innovations[3] = 0
 AR_model_1[4]: inflow[4,1] - 0.8725641041777851 inflow0[2,1]
   - 0.0864261426965268 inflow0[3,1] - 0.2463591011630234 inflow0[4,1]
   - innovations[4] = 0
 AR_model_1[5]: inflow[5,1] - 0.2663179396472554 inflow0[3,1]
   - 0.073005009895697 inflow0[4,1] - 0.456622635088196 inflow0[5,1]
   - innovations[5] = 0
 AR_model_1[6]: inflow[6,1] + 0.0410316527531602 inflow0[4,1]
   - 0.2988158146315212 inflow0[5,1] - 0.0667348654748534 inflow0[6,1]
   - innovations[6] = 0
 AR_model_1[7]: inflow[7,1] + 0.5105979631668157 inflow0[5,1]
   + 0.170819147860345 inflow0[6,1] - 0.0626135403301803 inflow0[7,1]
   - innovations[7] = 0
 AR_model_1[8]: inflow[8,1] + 0.1484330040812884 inflow0[6,1]
   - 0.1728872397719012 inflow0[7,1] - 0.1307237190880721 inflow0[8,1]
   - innovations[8] = 0
 AR_model_1[9]: inflow[9,1] - 0.00915170346028921 inflow0[7,1]
   + 0.290873134919753 inflow0[8,1] - 0.4539509247975215 inflow0[9,1]
   - innovations[9] = 0
 balance[0]: reservoir_level[0] - reservoir_level0[0] - inflow[0,1]
   + outflow[0] + spill[0] - pour[0] = 0
 balance[1]: reservoir_level[1] - reservoir_level0[1] - inflow[1,1]
   - outflow[0] + outflow[1] - spill[0] + spill[1] - pour[1] = 0
 balance[2]: reservoir_level[2] - reservoir_level0[2] - inflow[2,1]
   - outflow[1] + outflow[2] - spill[1] + spill[2] - pour[2] = 0
 balance[3]: reservoir_level[3] - reservoir_level0[3] - inflow[3,1]
   - outflow[2] + outflow[3] - spill[2] + spill[3] - pour[3] = 0
 balance[4]: reservoir_level[4] - reservoir_level0[4] - inflow[4,1]
   - outflow[3] + outflow[4] - spill[3] + spill[4] - pour[4] = 0
 balance[5]: reservoir_level[5] - reservoir_level0[5] - inflow[5,1]
   - outflow[4] + outflow[5] - spill[4] + spill[5] - pour[5] = 0
 balance[6]: reservoir_level[6] - reservoir_level0[6] - inflow[6,1]
   - outflow[5] + outflow[6] - spill[5] + spill[6] - pour[6] = 0
 balance[7]: reservoir_level[7] - reservoir_level0[7] - inflow[7,1]
   - outflow[6] + outflow[7] - spill[6] + spill[7] - pour[7] = 0
 balance[8]: reservoir_level[8] - reservoir_level0[8] - inflow[8,1]
   - outflow[7] + outflow[8] - spill[7] + spill[8] - pour[8] = 0
 balance[9]: reservoir_level[9] - reservoir_level0[9] - inflow[9,1]
   - outflow[8] + outflow[9] - spill[8] + spill[9] - pour[9] = 0
 generationCtr: generation - 55 dispatch[0,0] - 65 dispatch[0,1]
   - 70 dispatch[0,2] - 55 dispatch[1,0] - 65 dispatch[1,1]
   - 70 dispatch[1,2] - 55 dispatch[2,0] - 65 dispatch[2,1]
   - 70 dispatch[2,2] - 55 dispatch[3,0] - 65 dispatch[3,1]
   - 70 dispatch[3,2] - 55 dispatch[4,0] - 65 dispatch[4,1]
   - 70 dispatch[4,2] - 55 dispatch[5,0] - 65 dispatch[5,1]
   - 70 dispatch[5,2] - 55 dispatch[6,0] - 65 dispatch[6,1]
   - 70 dispatch[6,2] - 55 dispatch[7,0] - 65 dispatch[7,1]
   - 70 dispatch[7,2] - 55 dispatch[8,0] - 65 dispatch[8,1]
   - 70 dispatch[8,2] - 55 dispatch[9,0] - 65 dispatch[9,1]
   - 70 dispatch[9,2] = 0
 outflowCtr[0]: outflow[0] - 50 dispatch[0,0] - 60 dispatch[0,1]
   - 70 dispatch[0,2] = 0
 outflowCtr[1]: outflow[1] - 50 dispatch[1,0] - 60 dispatch[1,1]
   - 70 dispatch[1,2] = 0
 outflowCtr[2]: outflow[2] - 50 dispatch[2,0] - 60 dispatch[2,1]
   - 70 dispatch[2,2] = 0
 outflowCtr[3]: outflow[3] - 50 dispatch[3,0] - 60 dispatch[3,1]
   - 70 dispatch[3,2] = 0
 outflowCtr[4]: outflow[4] - 50 dispatch[4,0] - 60 dispatch[4,1]
   - 70 dispatch[4,2] = 0
 outflowCtr[5]: outflow[5] - 50 dispatch[5,0] - 60 dispatch[5,1]
   - 70 dispatch[5,2] = 0
 outflowCtr[6]: outflow[6] - 50 dispatch[6,0] - 60 dispatch[6,1]
   - 70 dispatch[6,2] = 0
 outflowCtr[7]: outflow[7] - 50 dispatch[7,0] - 60 dispatch[7,1]
   - 70 dispatch[7,2] = 0
 outflowCtr[8]: outflow[8] - 50 dispatch[8,0] - 60 dispatch[8,1]
   - 70 dispatch[8,2] = 0
 outflowCtr[9]: outflow[9] - 50 dispatch[9,0] - 60 dispatch[9,1]
   - 70 dispatch[9,2] = 0
 dispatchCtr[0]: dispatch[0,0] + dispatch[0,1] + dispatch[0,2] <= 1
 dispatchCtr[1]: dispatch[1,0] + dispatch[1,1] + dispatch[1,2] <= 1
 dispatchCtr[2]: dispatch[2,0] + dispatch[2,1] + dispatch[2,2] <= 1
 dispatchCtr[3]: dispatch[3,0] + dispatch[3,1] + dispatch[3,2] <= 1
 dispatchCtr[4]: dispatch[4,0] + dispatch[4,1] + dispatch[4,2] <= 1
 dispatchCtr[5]: dispatch[5,0] + dispatch[5,1] + dispatch[5,2] <= 1
 dispatchCtr[6]: dispatch[6,0] + dispatch[6,1] + dispatch[6,2] <= 1
 dispatchCtr[7]: dispatch[7,0] + dispatch[7,1] + dispatch[7,2] <= 1
 dispatchCtr[8]: dispatch[8,0] + dispatch[8,1] + dispatch[8,2] <= 1
 dispatchCtr[9]: dispatch[9,0] + dispatch[9,1] + dispatch[9,2] <= 1
Bounds
 30 <= reservoir_level[0] <= 200
 30 <= reservoir_level[1] <= 200
 30 <= reservoir_level[2] <= 200
 30 <= reservoir_level[3] <= 200
 30 <= reservoir_level[4] <= 200
 30 <= reservoir_level[5] <= 200
 30 <= reservoir_level[6] <= 200
 30 <= reservoir_level[7] <= 200
 30 <= reservoir_level[8] <= 200
 30 <= reservoir_level[9] <= 200
 reservoir_level0[0] = 50
 reservoir_level0[1] = 50
 reservoir_level0[2] = 50
 reservoir_level0[3] = 50
 reservoir_level0[4] = 50
 reservoir_level0[5] = 50
 reservoir_level0[6] = 50
 reservoir_level0[7] = 50
 reservoir_level0[8] = 50
 reservoir_level0[9] = 50
 inflow[0,1] = 5.367259049368
 inflow[1,1] = 7.51492688285188
 inflow[2,1] = 10.25194402758205
 inflow[3,1] = 10.6905304047435
 inflow[4,1] = 5.553761321855924
 inflow[5,1] = 8.417942790984373
 inflow[6,1] = 7.546247479617146
 inflow[7,1] = 0.2065267438602336
 inflow[8,1] = 11.59239932418561
 inflow[9,1] = 13.23961785954175
 inflow0[0,1] = 2.872791755683384
 inflow0[1,1] = 9.331631565597478
 inflow0[2,1] = 6.565916085106717
 inflow0[3,1] = 11.78037875570654
 inflow0[4,1] = 11.69963712178205
 inflow0[5,1] = 4.088889079239625
 inflow0[6,1] = 4.14696382714645
 inflow0[7,1] = 12.02808266302529
 inflow0[8,1] = 14.37209030525558
 inflow0[9,1] = 13.13898952113142
 innovations[0] free
 innovations[1] free
 innovations[2] free
 innovations[3] free
 innovations[4] free
 innovations[5] free
 innovations[6] free
 innovations[7] free
 innovations[8] free
 innovations[9] free
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
 dispatch[4,0] <= 1
 dispatch[4,1] <= 1
 dispatch[4,2] <= 1
 dispatch[5,0] <= 1
 dispatch[5,1] <= 1
 dispatch[5,2] <= 1
 dispatch[6,0] <= 1
 dispatch[6,1] <= 1
 dispatch[6,2] <= 1
 dispatch[7,0] <= 1
 dispatch[7,1] <= 1
 dispatch[7,2] <= 1
 dispatch[8,0] <= 1
 dispatch[8,1] <= 1
 dispatch[8,2] <= 1
 dispatch[9,0] <= 1
 dispatch[9,1] <= 1
 dispatch[9,2] <= 1
 oracle[0][0] >= -1e+08
 oracle[0][1] >= -1e+08
 oracle[0][2] >= -1e+08
 oracle[0][3] >= -1e+08
 oracle[0][4] >= -1e+08
 lambda[0] free
End
