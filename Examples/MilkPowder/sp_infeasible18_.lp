\ Model MilkPowder_18
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.0001 W + 0.5 b + 0.275 h + 0.5 i + C19 + 100 C20 + 20 gamma[18]
   + 0.5 nu[18][0] + 0.5 nu[18][1]
Subject To
 Pasture_Cover: P - P0 + h + fp - 7 gr <= 0
 Supplement_Storage: Q - Q0 - 0.9 h + fs <= 0
 Supplement_Storage: Q - Q0 - 0.9 h + fs <= 0
 CowsAmount: C - C0 <= 0
 Soil_Limit: W <= 150
 milkRel: - 0.0122360501055756 mlk + milk <= 0
 milk_BALANCE: M - M0 - milk <= 0
 EnergyBalance: 3.523269419982503 C0 + 11 b + 11 fs + 11 fp - mlk
   >= 1181.029975391504
 maximum_milk: - 1057.631214988625 C0 + mlk <= 0
 minimum_milk: - 6.3 C0 + milk >= 0
 pasture_growth_constraints1: - 1.488939695327892 ev + gr <= 0
 growth_aprox0.000000: - 0.0742857142857143 P0 + gr
   <= 7.4285714285714287e-04
 growth_aprox184.210526: - 0.0664661654135338 P0 + gr
   <= 0.7208862683023334
 growth_aprox368.421053: - 0.0586466165413534 P0 + gr <= 2.881472892758211
 growth_aprox552.631579: - 0.0508270676691729 P0 + gr <= 6.482502730510483
 growth_aprox736.842105: - 0.0430075187969925 P0 + gr <= 11.52397578155916
 growth_aprox921.052632: - 0.035187969924812 P0 + gr <= 18.00589204590423
 growth_aprox1105.263158: - 0.0273684210526316 P0 + gr <= 25.9282515235457
 growth_aprox1289.473684: - 0.0195488721804511 P0 + gr
   <= 35.29105421448357
 growth_aprox1473.684211: - 0.0117293233082707 P0 + gr
   <= 46.09430011871787
 growth_aprox1657.894737: - 0.00390977443609022 P0 + gr
   <= 58.33798923624853
 growth_aprox1842.105263: 0.00390977443609022 P0 + gr <= 72.02212156707559
 growth_aprox2026.315789: 0.0117293233082707 P0 + gr <= 87.14669711119906
 growth_aprox2210.526316: 0.0195488721804511 P0 + gr <= 103.7117158686189
 growth_aprox2394.736842: 0.0273684210526316 P0 + gr <= 121.7171778393352
 growth_aprox2578.947368: 0.035187969924812 P0 + gr <= 141.1630830233478
 growth_aprox2763.157895: 0.0430075187969925 P0 + gr <= 162.0494314206569
 growth_aprox2947.368421: 0.0508270676691729 P0 + gr <= 184.3762230312624
 growth_aprox3131.578947: 0.0586466165413534 P0 + gr <= 208.1434578551643
 growth_aprox3315.789474: 0.0664661654135338 P0 + gr <= 233.3511358923625
 growth_aprox3500.000000: 0.0742857142857143 P0 + gr <= 259.9992571428571
 max_irrigation: i <= 0
 evapotranspiration_ctr: ev - evapotranspiration <= 0
 soil-mosture-balance: W - W0 - i + ev - rainfall <= 0
 max_rate_d1_1: C19 >= 0
 max_rate_d1_2: - 0.25 b + C19 >= -15.75
 max_rate_d1_3: - 0.75 b + C19 >= -57.75
 max_rate_d1_4: - b + C19 >= -89.25
 dro_dual_ctr[00]: - oracle[18][0] + nu[18][0] >= 0
 dro_dual_ctr[01]: - oracle[18][1] + 4.9 gamma[18] + nu[18][0] >= 0
 dro_dual_ctr[10]: - oracle[18][0] + 4.9 gamma[18] + nu[18][1] >= 0
 dro_dual_ctr[11]: - oracle[18][1] + nu[18][1] >= 0
 cut[18,0,0]: 1.0019323569533771e+08 P + 52.29185759411531 C
   + oracle[18][0] >= -687908.8903603047
 cut[18,0,1]: 1.0019323569533771e+08 P + 52.29185759411531 C
   + oracle[18][1] >= -687908.8903603047
Bounds
 M = -10000
 P0 = 0
 Q0 = 0
 W0 = 150
 C0 = 0
 M0 = -10000
 evapotranspiration = 22.3
 rainfall = 26.23333333333333
 oracle[18][0] >= -1e+08
 oracle[18][1] >= -1e+08
 nu[18][0] free
 nu[18][1] free
End
