\ Model MilkPowder_0
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.0001 W + 0.5 b + 0.275 h + 0.5 i + C19 + 100 C20 + 0 oracle[0][0]
   + 0 oracle[0][1] + 0 oracle[0][2] + 0 oracle[0][3] + 0 oracle[0][4]
   + 0 oracle[0][5] + 0 oracle[0][6] + 0 oracle[0][7] + 0 oracle[0][8]
   + 0 oracle[0][9] + 0 oracle[0][10] + 0 oracle[0][11] + 0 oracle[0][12]
   + 0 oracle[0][13] + 0 oracle[0][14] + 0 oracle[0][15] + 0 oracle[0][16]
   + 0 oracle[0][17] + 0 oracle[0][18] + 0 oracle[0][19] + GlobalOracle[0]
Subject To
 Pasture_Cover: P - P0 + h + fp - 7 gr <= 0
 Supplement_Storage: Q - Q0 - 0.9 h + fs <= 0
 Supplement_Storage: Q - Q0 - 0.9 h + fs <= 0
 CowsAmount: C - C0 <= 0
 Soil_Limit: W <= 150
 milkRel: - 0.012750843874272 mlk + milk <= 0
 milk_BALANCE: M - M0 - milk <= 0
 EnergyBalance: 16.87801873849773 C0 + 11 b + 11 fs + 11 fp - mlk
   >= 916.9969019336006
 maximum_milk: - 888.195937634907 C0 + mlk <= 0
 minimum_milk: - 6.3 C0 + milk >= 0
 pasture_growth_constraints1: - 5.328493184769218 ev + gr <= 0
 growth_aprox0.000000: - 0.0742857142857143 P0 + gr
   <= 7.4285714285714287e-04
 growth_aprox1750.000000: gr <= 65
 growth_aprox3500.000000: 0.0742857142857143 P0 + gr <= 259.9992571428571
 max_irrigation: i <= 0
 evapotranspiration_ctr: ev - evapotranspiration <= 0
 soil-mosture-balance: W - W0 - i + ev - rainfall <= 0
 max_rate_d1_1: C19 >= 0
 max_rate_d1_2: - 0.25 b + C19 >= -15.75
 max_rate_d1_3: - 0.75 b + C19 >= -57.75
 max_rate_d1_4: - b + C19 >= -89.25
Bounds
 M = -10000
 P0 = 2500
 Q0 = 0
 W0 = 150
 C0 = 3
 M0 = 0
 evapotranspiration = 3.755
 rainfall = 58.55333333333333
 oracle[0][0] >= -1e+08
 oracle[0][1] >= -1e+08
 oracle[0][2] >= -1e+08
 oracle[0][3] >= -1e+08
 oracle[0][4] >= -1e+08
 oracle[0][5] >= -1e+08
 oracle[0][6] >= -1e+08
 oracle[0][7] >= -1e+08
 oracle[0][8] >= -1e+08
 oracle[0][9] >= -1e+08
 oracle[0][10] >= -1e+08
 oracle[0][11] >= -1e+08
 oracle[0][12] >= -1e+08
 oracle[0][13] >= -1e+08
 oracle[0][14] >= -1e+08
 oracle[0][15] >= -1e+08
 oracle[0][16] >= -1e+08
 oracle[0][17] >= -1e+08
 oracle[0][18] >= -1e+08
 oracle[0][19] >= -1e+08
 GlobalOracle[0] free
End
