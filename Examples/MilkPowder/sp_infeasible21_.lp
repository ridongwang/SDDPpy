\ Model MilkPowder_21
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.0001 W + 0.5 b + 0.275 h + 0.5 i + C19 + 100 C20 + 10 gamma[21]
   + 0.05 nu[21][0] + 0.05 nu[21][1] + 0.05 nu[21][2] + 0.05 nu[21][3]
   + 0.05 nu[21][4] + 0.05 nu[21][5] + 0.05 nu[21][6] + 0.05 nu[21][7]
   + 0.05 nu[21][8] + 0.05 nu[21][9] + 0.05 nu[21][10] + 0.05 nu[21][11]
   + 0.05 nu[21][12] + 0.05 nu[21][13] + 0.05 nu[21][14] + 0.05 nu[21][15]
   + 0.05 nu[21][16] + 0.05 nu[21][17] + 0.05 nu[21][18] + 0.05 nu[21][19]
Subject To
 Pasture_Cover: P - P0 + h + fp - 7 gr <= 0
 Supplement_Storage: Q - Q0 - 0.9 h + fs <= 0
 Supplement_Storage: Q - Q0 - 0.9 h + fs <= 0
 CowsAmount: C - C0 <= 0
 Soil_Limit: W <= 150
 milkRel: - 0.0123071189627373 mlk + milk <= 0
 milk_BALANCE: M - M0 - milk <= 0
 EnergyBalance: 5.064699791224847 C0 + 11 b + 11 fs + 11 fp - mlk
   >= 1202.612434189625
 maximum_milk: - 1035.854952978385 C0 + mlk <= 0
 minimum_milk: - 6.3 C0 + milk >= 0
 pasture_growth_constraints1: - 1.17995735322005 ev + gr <= 0
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
 dro_dual_ctr[00]: - oracle[21][0] + nu[21][0] >= 0
 dro_dual_ctr[01]: - oracle[21][1] + 44.06666666666667 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[02]: - oracle[21][2] + 29.46666666666667 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[03]: - oracle[21][3] + 20.23333333333333 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[04]: - oracle[21][4] + 17 gamma[21] + nu[21][0] >= 0
 dro_dual_ctr[05]: - oracle[21][5] + 21.93333333333333 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[06]: - oracle[21][6] + 6.6 gamma[21] + nu[21][0] >= 0
 dro_dual_ctr[07]: - oracle[21][7] + 43.06666666666666 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[08]: - oracle[21][8] + 21.96666666666667 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[09]: - oracle[21][9] + 9.8 gamma[21] + nu[21][0] >= 0
 dro_dual_ctr[010]: - oracle[21][10] + 18.83333333333334 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[011]: - oracle[21][11] + 30.03333333333332 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[012]: - oracle[21][12] + 29.03333333333334 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[013]: - oracle[21][13] + 24.36666666666667 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[014]: - oracle[21][14] + 20.13333333333333 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[015]: - oracle[21][15] + 61.66666666666667 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[016]: - oracle[21][16] + 11.3 gamma[21] + nu[21][0] >= 0
 dro_dual_ctr[017]: - oracle[21][17] + 20.1 gamma[21] + nu[21][0] >= 0
 dro_dual_ctr[018]: - oracle[21][18] + 38.26666666666667 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[019]: - oracle[21][19] + 21.53333333333333 gamma[21]
   + nu[21][0] >= 0
 dro_dual_ctr[10]: - oracle[21][0] + 44.06666666666667 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[11]: - oracle[21][1] + nu[21][1] >= 0
 dro_dual_ctr[12]: - oracle[21][2] + 20.6 gamma[21] + nu[21][1] >= 0
 dro_dual_ctr[13]: - oracle[21][3] + 23.83333333333334 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[14]: - oracle[21][4] + 27.06666666666667 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[15]: - oracle[21][5] + 30.93333333333333 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[16]: - oracle[21][6] + 37.46666666666667 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[17]: - oracle[21][7] + gamma[21] + nu[21][1] >= 0
 dro_dual_ctr[18]: - oracle[21][8] + 42.23333333333334 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[19]: - oracle[21][9] + 34.26666666666667 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[110]: - oracle[21][10] + 25.23333333333333 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[111]: - oracle[21][11] + 74.1 gamma[21] + nu[21][1] >= 0
 dro_dual_ctr[112]: - oracle[21][12] + 19.63333333333333 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[113]: - oracle[21][13] + 19.7 gamma[21] + nu[21][1] >= 0
 dro_dual_ctr[114]: - oracle[21][14] + 26.13333333333333 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[115]: - oracle[21][15] + 87.93333333333334 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[116]: - oracle[21][16] + 51.96666666666666 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[117]: - oracle[21][17] + 35.76666666666667 gamma[21]
   + nu[21][1] >= 0
 dro_dual_ctr[118]: - oracle[21][18] + 5.8 gamma[21] + nu[21][1] >= 0
 dro_dual_ctr[119]: - oracle[21][19] + 43.2 gamma[21] + nu[21][1] >= 0
 dro_dual_ctr[20]: - oracle[21][0] + 29.46666666666667 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[21]: - oracle[21][1] + 20.6 gamma[21] + nu[21][2] >= 0
 dro_dual_ctr[22]: - oracle[21][2] + nu[21][2] >= 0
 dro_dual_ctr[23]: - oracle[21][3] + 25.83333333333334 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[24]: - oracle[21][4] + 12.86666666666666 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[25]: - oracle[21][5] + 10.33333333333333 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[26]: - oracle[21][6] + 36.06666666666667 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[27]: - oracle[21][7] + 19.6 gamma[21] + nu[21][2] >= 0
 dro_dual_ctr[28]: - oracle[21][8] + 51.43333333333334 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[29]: - oracle[21][9] + 34.86666666666667 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[210]: - oracle[21][10] + 25.63333333333333 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[211]: - oracle[21][11] + 58.5 gamma[21] + nu[21][2] >= 0
 dro_dual_ctr[212]: - oracle[21][12] + 40.23333333333333 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[213]: - oracle[21][13] + 8.5 gamma[21] + nu[21][2] >= 0
 dro_dual_ctr[214]: - oracle[21][14] + 9.333333333333339 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[215]: - oracle[21][15] + 91.13333333333334 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[216]: - oracle[21][16] + 40.76666666666667 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[217]: - oracle[21][17] + 15.16666666666667 gamma[21]
   + nu[21][2] >= 0
 dro_dual_ctr[218]: - oracle[21][18] + 24.8 gamma[21] + nu[21][2] >= 0
 dro_dual_ctr[219]: - oracle[21][19] + 51 gamma[21] + nu[21][2] >= 0
 dro_dual_ctr[30]: - oracle[21][0] + 20.23333333333333 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[31]: - oracle[21][1] + 23.83333333333334 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[32]: - oracle[21][2] + 25.83333333333334 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[33]: - oracle[21][3] + nu[21][3] >= 0
 dro_dual_ctr[34]: - oracle[21][4] + 12.96666666666667 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[35]: - oracle[21][5] + 18.3 gamma[21] + nu[21][3] >= 0
 dro_dual_ctr[36]: - oracle[21][6] + 13.63333333333333 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[37]: - oracle[21][7] + 22.83333333333333 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[38]: - oracle[21][8] + 25.6 gamma[21] + nu[21][3] >= 0
 dro_dual_ctr[39]: - oracle[21][9] + 10.43333333333333 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[310]: - oracle[21][10] + 1.4 gamma[21] + nu[21][3] >= 0
 dro_dual_ctr[311]: - oracle[21][11] + 50.26666666666667 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[312]: - oracle[21][12] + 14.4 gamma[21] + nu[21][3] >= 0
 dro_dual_ctr[313]: - oracle[21][13] + 17.33333333333333 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[314]: - oracle[21][14] + 16.5 gamma[21] + nu[21][3] >= 0
 dro_dual_ctr[315]: - oracle[21][15] + 65.3 gamma[21] + nu[21][3] >= 0
 dro_dual_ctr[316]: - oracle[21][16] + 28.13333333333333 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[317]: - oracle[21][17] + 16.46666666666667 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[318]: - oracle[21][18] + 18.03333333333333 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[319]: - oracle[21][19] + 25.16666666666667 gamma[21]
   + nu[21][3] >= 0
 dro_dual_ctr[40]: - oracle[21][0] + 17 gamma[21] + nu[21][4] >= 0
 dro_dual_ctr[41]: - oracle[21][1] + 27.06666666666667 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[42]: - oracle[21][2] + 12.86666666666666 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[43]: - oracle[21][3] + 12.96666666666667 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[44]: - oracle[21][4] + nu[21][4] >= 0
 dro_dual_ctr[45]: - oracle[21][5] + 5.333333333333329 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[46]: - oracle[21][6] + 23.2 gamma[21] + nu[21][4] >= 0
 dro_dual_ctr[47]: - oracle[21][7] + 26.06666666666667 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[48]: - oracle[21][8] + 38.56666666666668 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[49]: - oracle[21][9] + 22 gamma[21] + nu[21][4] >= 0
 dro_dual_ctr[410]: - oracle[21][10] + 12.76666666666667 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[411]: - oracle[21][11] + 47.03333333333332 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[412]: - oracle[21][12] + 27.36666666666667 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[413]: - oracle[21][13] + 7.366666666666675 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[414]: - oracle[21][14] + 3.533333333333325 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[415]: - oracle[21][15] + 78.26666666666667 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[416]: - oracle[21][16] + 27.9 gamma[21] + nu[21][4] >= 0
 dro_dual_ctr[417]: - oracle[21][17] + 8.7 gamma[21] + nu[21][4] >= 0
 dro_dual_ctr[418]: - oracle[21][18] + 21.26666666666667 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[419]: - oracle[21][19] + 38.13333333333334 gamma[21]
   + nu[21][4] >= 0
 dro_dual_ctr[50]: - oracle[21][0] + 21.93333333333333 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[51]: - oracle[21][1] + 30.93333333333333 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[52]: - oracle[21][2] + 10.33333333333333 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[53]: - oracle[21][3] + 18.3 gamma[21] + nu[21][5] >= 0
 dro_dual_ctr[54]: - oracle[21][4] + 5.333333333333329 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[55]: - oracle[21][5] + nu[21][5] >= 0
 dro_dual_ctr[56]: - oracle[21][6] + 28.53333333333333 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[57]: - oracle[21][7] + 29.93333333333333 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[58]: - oracle[21][8] + 43.9 gamma[21] + nu[21][5] >= 0
 dro_dual_ctr[59]: - oracle[21][9] + 27.33333333333333 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[510]: - oracle[21][10] + 18.1 gamma[21] + nu[21][5] >= 0
 dro_dual_ctr[511]: - oracle[21][11] + 50.96666666666667 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[512]: - oracle[21][12] + 32.7 gamma[21] + nu[21][5] >= 0
 dro_dual_ctr[513]: - oracle[21][13] + 11.23333333333333 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[514]: - oracle[21][14] + 4.8 gamma[21] + nu[21][5] >= 0
 dro_dual_ctr[515]: - oracle[21][15] + 83.6 gamma[21] + nu[21][5] >= 0
 dro_dual_ctr[516]: - oracle[21][16] + 33.23333333333333 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[517]: - oracle[21][17] + 4.833333333333334 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[518]: - oracle[21][18] + 25.13333333333333 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[519]: - oracle[21][19] + 43.46666666666667 gamma[21]
   + nu[21][5] >= 0
 dro_dual_ctr[60]: - oracle[21][0] + 6.6 gamma[21] + nu[21][6] >= 0
 dro_dual_ctr[61]: - oracle[21][1] + 37.46666666666667 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[62]: - oracle[21][2] + 36.06666666666667 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[63]: - oracle[21][3] + 13.63333333333333 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[64]: - oracle[21][4] + 23.2 gamma[21] + nu[21][6] >= 0
 dro_dual_ctr[65]: - oracle[21][5] + 28.53333333333333 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[66]: - oracle[21][6] + nu[21][6] >= 0
 dro_dual_ctr[67]: - oracle[21][7] + 36.46666666666667 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[68]: - oracle[21][8] + 15.36666666666667 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[69]: - oracle[21][9] + 3.2 gamma[21] + nu[21][6] >= 0
 dro_dual_ctr[610]: - oracle[21][10] + 12.23333333333333 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[611]: - oracle[21][11] + 36.63333333333333 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[612]: - oracle[21][12] + 22.43333333333334 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[613]: - oracle[21][13] + 27.56666666666667 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[614]: - oracle[21][14] + 26.73333333333333 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[615]: - oracle[21][15] + 55.06666666666666 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[616]: - oracle[21][16] + 14.5 gamma[21] + nu[21][6] >= 0
 dro_dual_ctr[617]: - oracle[21][17] + 26.7 gamma[21] + nu[21][6] >= 0
 dro_dual_ctr[618]: - oracle[21][18] + 31.66666666666667 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[619]: - oracle[21][19] + 14.93333333333333 gamma[21]
   + nu[21][6] >= 0
 dro_dual_ctr[70]: - oracle[21][0] + 43.06666666666666 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[71]: - oracle[21][1] + gamma[21] + nu[21][7] >= 0
 dro_dual_ctr[72]: - oracle[21][2] + 19.6 gamma[21] + nu[21][7] >= 0
 dro_dual_ctr[73]: - oracle[21][3] + 22.83333333333333 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[74]: - oracle[21][4] + 26.06666666666667 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[75]: - oracle[21][5] + 29.93333333333333 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[76]: - oracle[21][6] + 36.46666666666667 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[77]: - oracle[21][7] + nu[21][7] >= 0
 dro_dual_ctr[78]: - oracle[21][8] + 41.23333333333333 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[79]: - oracle[21][9] + 33.26666666666667 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[710]: - oracle[21][10] + 24.23333333333333 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[711]: - oracle[21][11] + 73.1 gamma[21] + nu[21][7] >= 0
 dro_dual_ctr[712]: - oracle[21][12] + 20.63333333333334 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[713]: - oracle[21][13] + 18.7 gamma[21] + nu[21][7] >= 0
 dro_dual_ctr[714]: - oracle[21][14] + 25.13333333333333 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[715]: - oracle[21][15] + 86.93333333333334 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[716]: - oracle[21][16] + 50.96666666666665 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[717]: - oracle[21][17] + 34.76666666666667 gamma[21]
   + nu[21][7] >= 0
 dro_dual_ctr[718]: - oracle[21][18] + 5.2 gamma[21] + nu[21][7] >= 0
 dro_dual_ctr[719]: - oracle[21][19] + 42.2 gamma[21] + nu[21][7] >= 0
 dro_dual_ctr[80]: - oracle[21][0] + 21.96666666666667 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[81]: - oracle[21][1] + 42.23333333333334 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[82]: - oracle[21][2] + 51.43333333333334 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[83]: - oracle[21][3] + 25.6 gamma[21] + nu[21][8] >= 0
 dro_dual_ctr[84]: - oracle[21][4] + 38.56666666666668 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[85]: - oracle[21][5] + 43.9 gamma[21] + nu[21][8] >= 0
 dro_dual_ctr[86]: - oracle[21][6] + 15.36666666666667 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[87]: - oracle[21][7] + 41.23333333333333 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[88]: - oracle[21][8] + nu[21][8] >= 0
 dro_dual_ctr[89]: - oracle[21][9] + 16.56666666666667 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[810]: - oracle[21][10] + 25.8 gamma[21] + nu[21][8] >= 0
 dro_dual_ctr[811]: - oracle[21][11] + 31.86666666666665 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[812]: - oracle[21][12] + 27.2 gamma[21] + nu[21][8] >= 0
 dro_dual_ctr[813]: - oracle[21][13] + 42.93333333333334 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[814]: - oracle[21][14] + 42.1 gamma[21] + nu[21][8] >= 0
 dro_dual_ctr[815]: - oracle[21][15] + 45.7 gamma[21] + nu[21][8] >= 0
 dro_dual_ctr[816]: - oracle[21][16] + 10.66666666666667 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[817]: - oracle[21][17] + 42.06666666666667 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[818]: - oracle[21][18] + 36.43333333333334 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[819]: - oracle[21][19] + 0.9666666666666544 gamma[21]
   + nu[21][8] >= 0
 dro_dual_ctr[90]: - oracle[21][0] + 9.8 gamma[21] + nu[21][9] >= 0
 dro_dual_ctr[91]: - oracle[21][1] + 34.26666666666667 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[92]: - oracle[21][2] + 34.86666666666667 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[93]: - oracle[21][3] + 10.43333333333333 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[94]: - oracle[21][4] + 22 gamma[21] + nu[21][9] >= 0
 dro_dual_ctr[95]: - oracle[21][5] + 27.33333333333333 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[96]: - oracle[21][6] + 3.2 gamma[21] + nu[21][9] >= 0
 dro_dual_ctr[97]: - oracle[21][7] + 33.26666666666667 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[98]: - oracle[21][8] + 16.56666666666667 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[99]: - oracle[21][9] + nu[21][9] >= 0
 dro_dual_ctr[910]: - oracle[21][10] + 9.233333333333334 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[911]: - oracle[21][11] + 39.83333333333333 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[912]: - oracle[21][12] + 19.23333333333333 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[913]: - oracle[21][13] + 26.36666666666666 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[914]: - oracle[21][14] + 25.53333333333333 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[915]: - oracle[21][15] + 56.26666666666667 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[916]: - oracle[21][16] + 17.7 gamma[21] + nu[21][9] >= 0
 dro_dual_ctr[917]: - oracle[21][17] + 25.5 gamma[21] + nu[21][9] >= 0
 dro_dual_ctr[918]: - oracle[21][18] + 28.46666666666667 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[919]: - oracle[21][19] + 16.13333333333333 gamma[21]
   + nu[21][9] >= 0
 dro_dual_ctr[100]: - oracle[21][0] + 18.83333333333334 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[101]: - oracle[21][1] + 25.23333333333333 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[102]: - oracle[21][2] + 25.63333333333333 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[103]: - oracle[21][3] + 1.4 gamma[21] + nu[21][10] >= 0
 dro_dual_ctr[104]: - oracle[21][4] + 12.76666666666667 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[105]: - oracle[21][5] + 18.1 gamma[21] + nu[21][10] >= 0
 dro_dual_ctr[106]: - oracle[21][6] + 12.23333333333333 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[107]: - oracle[21][7] + 24.23333333333333 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[108]: - oracle[21][8] + 25.8 gamma[21] + nu[21][10] >= 0
 dro_dual_ctr[109]: - oracle[21][9] + 9.233333333333334 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[1010]: - oracle[21][10] + nu[21][10] >= 0
 dro_dual_ctr[1011]: - oracle[21][11] + 48.86666666666666 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[1012]: - oracle[21][12] + 14.6 gamma[21] + nu[21][10] >= 0
 dro_dual_ctr[1013]: - oracle[21][13] + 17.13333333333333 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[1014]: - oracle[21][14] + 16.3 gamma[21] + nu[21][10] >= 0
 dro_dual_ctr[1015]: - oracle[21][15] + 65.5 gamma[21] + nu[21][10] >= 0
 dro_dual_ctr[1016]: - oracle[21][16] + 26.73333333333333 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[1017]: - oracle[21][17] + 16.26666666666667 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[1018]: - oracle[21][18] + 19.43333333333333 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[1019]: - oracle[21][19] + 25.36666666666667 gamma[21]
   + nu[21][10] >= 0
 dro_dual_ctr[110]: - oracle[21][0] + 30.03333333333332 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[111]: - oracle[21][1] + 74.1 gamma[21] + nu[21][11] >= 0
 dro_dual_ctr[112]: - oracle[21][2] + 58.5 gamma[21] + nu[21][11] >= 0
 dro_dual_ctr[113]: - oracle[21][3] + 50.26666666666667 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[114]: - oracle[21][4] + 47.03333333333332 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[115]: - oracle[21][5] + 50.96666666666667 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[116]: - oracle[21][6] + 36.63333333333333 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[117]: - oracle[21][7] + 73.1 gamma[21] + nu[21][11] >= 0
 dro_dual_ctr[118]: - oracle[21][8] + 31.86666666666665 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[119]: - oracle[21][9] + 39.83333333333333 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[1110]: - oracle[21][10] + 48.86666666666666 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[1111]: - oracle[21][11] + nu[21][11] >= 0
 dro_dual_ctr[1112]: - oracle[21][12] + 59.06666666666666 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[1113]: - oracle[21][13] + 54.4 gamma[21] + nu[21][11] >= 0
 dro_dual_ctr[1114]: - oracle[21][14] + 49.16666666666667 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[1115]: - oracle[21][15] + 32.63333333333333 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[1116]: - oracle[21][16] + 22.13333333333333 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[1117]: - oracle[21][17] + 49.13333333333334 gamma[21]
   + nu[21][11] >= 0
 dro_dual_ctr[1118]: - oracle[21][18] + 68.3 gamma[21] + nu[21][11] >= 0
 dro_dual_ctr[1119]: - oracle[21][19] + 30.9 gamma[21] + nu[21][11] >= 0
 dro_dual_ctr[120]: - oracle[21][0] + 29.03333333333334 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[121]: - oracle[21][1] + 19.63333333333333 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[122]: - oracle[21][2] + 40.23333333333333 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[123]: - oracle[21][3] + 14.4 gamma[21] + nu[21][12] >= 0
 dro_dual_ctr[124]: - oracle[21][4] + 27.36666666666667 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[125]: - oracle[21][5] + 32.7 gamma[21] + nu[21][12] >= 0
 dro_dual_ctr[126]: - oracle[21][6] + 22.43333333333334 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[127]: - oracle[21][7] + 20.63333333333334 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[128]: - oracle[21][8] + 27.2 gamma[21] + nu[21][12] >= 0
 dro_dual_ctr[129]: - oracle[21][9] + 19.23333333333333 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[1210]: - oracle[21][10] + 14.6 gamma[21] + nu[21][12] >= 0
 dro_dual_ctr[1211]: - oracle[21][11] + 59.06666666666666 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[1212]: - oracle[21][12] + nu[21][12] >= 0
 dro_dual_ctr[1213]: - oracle[21][13] + 31.73333333333333 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[1214]: - oracle[21][14] + 30.9 gamma[21] + nu[21][12] >= 0
 dro_dual_ctr[1215]: - oracle[21][15] + 72.9 gamma[21] + nu[21][12] >= 0
 dro_dual_ctr[1216]: - oracle[21][16] + 36.93333333333333 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[1217]: - oracle[21][17] + 30.86666666666667 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[1218]: - oracle[21][18] + 15.43333333333333 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[1219]: - oracle[21][19] + 28.16666666666666 gamma[21]
   + nu[21][12] >= 0
 dro_dual_ctr[130]: - oracle[21][0] + 24.36666666666667 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[131]: - oracle[21][1] + 19.7 gamma[21] + nu[21][13] >= 0
 dro_dual_ctr[132]: - oracle[21][2] + 8.5 gamma[21] + nu[21][13] >= 0
 dro_dual_ctr[133]: - oracle[21][3] + 17.33333333333333 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[134]: - oracle[21][4] + 7.366666666666675 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[135]: - oracle[21][5] + 11.23333333333333 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[136]: - oracle[21][6] + 27.56666666666667 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[137]: - oracle[21][7] + 18.7 gamma[21] + nu[21][13] >= 0
 dro_dual_ctr[138]: - oracle[21][8] + 42.93333333333334 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[139]: - oracle[21][9] + 26.36666666666666 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[1310]: - oracle[21][10] + 17.13333333333333 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[1311]: - oracle[21][11] + 54.4 gamma[21] + nu[21][13] >= 0
 dro_dual_ctr[1312]: - oracle[21][12] + 31.73333333333333 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[1313]: - oracle[21][13] + nu[21][13] >= 0
 dro_dual_ctr[1314]: - oracle[21][14] + 6.43333333333333 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[1315]: - oracle[21][15] + 82.63333333333334 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[1316]: - oracle[21][16] + 32.26666666666667 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[1317]: - oracle[21][17] + 16.06666666666667 gamma[21]
   + nu[21][13] >= 0
 dro_dual_ctr[1318]: - oracle[21][18] + 16.3 gamma[21] + nu[21][13] >= 0
 dro_dual_ctr[1319]: - oracle[21][19] + 42.5 gamma[21] + nu[21][13] >= 0
 dro_dual_ctr[140]: - oracle[21][0] + 20.13333333333333 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[141]: - oracle[21][1] + 26.13333333333333 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[142]: - oracle[21][2] + 9.333333333333339 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[143]: - oracle[21][3] + 16.5 gamma[21] + nu[21][14] >= 0
 dro_dual_ctr[144]: - oracle[21][4] + 3.533333333333325 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[145]: - oracle[21][5] + 4.8 gamma[21] + nu[21][14] >= 0
 dro_dual_ctr[146]: - oracle[21][6] + 26.73333333333333 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[147]: - oracle[21][7] + 25.13333333333333 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[148]: - oracle[21][8] + 42.1 gamma[21] + nu[21][14] >= 0
 dro_dual_ctr[149]: - oracle[21][9] + 25.53333333333333 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[1410]: - oracle[21][10] + 16.3 gamma[21] + nu[21][14] >= 0
 dro_dual_ctr[1411]: - oracle[21][11] + 49.16666666666667 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[1412]: - oracle[21][12] + 30.9 gamma[21] + nu[21][14] >= 0
 dro_dual_ctr[1413]: - oracle[21][13] + 6.43333333333333 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[1414]: - oracle[21][14] + nu[21][14] >= 0
 dro_dual_ctr[1415]: - oracle[21][15] + 81.8 gamma[21] + nu[21][14] >= 0
 dro_dual_ctr[1416]: - oracle[21][16] + 31.43333333333333 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[1417]: - oracle[21][17] + 9.63333333333334 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[1418]: - oracle[21][18] + 20.33333333333333 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[1419]: - oracle[21][19] + 41.66666666666666 gamma[21]
   + nu[21][14] >= 0
 dro_dual_ctr[150]: - oracle[21][0] + 61.66666666666667 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[151]: - oracle[21][1] + 87.93333333333334 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[152]: - oracle[21][2] + 91.13333333333334 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[153]: - oracle[21][3] + 65.3 gamma[21] + nu[21][15] >= 0
 dro_dual_ctr[154]: - oracle[21][4] + 78.26666666666667 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[155]: - oracle[21][5] + 83.6 gamma[21] + nu[21][15] >= 0
 dro_dual_ctr[156]: - oracle[21][6] + 55.06666666666666 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[157]: - oracle[21][7] + 86.93333333333334 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[158]: - oracle[21][8] + 45.7 gamma[21] + nu[21][15] >= 0
 dro_dual_ctr[159]: - oracle[21][9] + 56.26666666666667 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[1510]: - oracle[21][10] + 65.5 gamma[21] + nu[21][15] >= 0
 dro_dual_ctr[1511]: - oracle[21][11] + 32.63333333333333 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[1512]: - oracle[21][12] + 72.9 gamma[21] + nu[21][15] >= 0
 dro_dual_ctr[1513]: - oracle[21][13] + 82.63333333333334 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[1514]: - oracle[21][14] + 81.8 gamma[21] + nu[21][15] >= 0
 dro_dual_ctr[1515]: - oracle[21][15] + nu[21][15] >= 0
 dro_dual_ctr[1516]: - oracle[21][16] + 50.36666666666667 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[1517]: - oracle[21][17] + 81.76666666666667 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[1518]: - oracle[21][18] + 82.13333333333333 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[1519]: - oracle[21][19] + 44.73333333333333 gamma[21]
   + nu[21][15] >= 0
 dro_dual_ctr[160]: - oracle[21][0] + 11.3 gamma[21] + nu[21][16] >= 0
 dro_dual_ctr[161]: - oracle[21][1] + 51.96666666666666 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[162]: - oracle[21][2] + 40.76666666666667 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[163]: - oracle[21][3] + 28.13333333333333 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[164]: - oracle[21][4] + 27.9 gamma[21] + nu[21][16] >= 0
 dro_dual_ctr[165]: - oracle[21][5] + 33.23333333333333 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[166]: - oracle[21][6] + 14.5 gamma[21] + nu[21][16] >= 0
 dro_dual_ctr[167]: - oracle[21][7] + 50.96666666666665 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[168]: - oracle[21][8] + 10.66666666666667 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[169]: - oracle[21][9] + 17.7 gamma[21] + nu[21][16] >= 0
 dro_dual_ctr[1610]: - oracle[21][10] + 26.73333333333333 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[1611]: - oracle[21][11] + 22.13333333333333 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[1612]: - oracle[21][12] + 36.93333333333333 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[1613]: - oracle[21][13] + 32.26666666666667 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[1614]: - oracle[21][14] + 31.43333333333333 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[1615]: - oracle[21][15] + 50.36666666666667 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[1616]: - oracle[21][16] + nu[21][16] >= 0
 dro_dual_ctr[1617]: - oracle[21][17] + 31.4 gamma[21] + nu[21][16] >= 0
 dro_dual_ctr[1618]: - oracle[21][18] + 46.16666666666666 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[1619]: - oracle[21][19] + 10.23333333333333 gamma[21]
   + nu[21][16] >= 0
 dro_dual_ctr[170]: - oracle[21][0] + 20.1 gamma[21] + nu[21][17] >= 0
 dro_dual_ctr[171]: - oracle[21][1] + 35.76666666666667 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[172]: - oracle[21][2] + 15.16666666666667 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[173]: - oracle[21][3] + 16.46666666666667 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[174]: - oracle[21][4] + 8.7 gamma[21] + nu[21][17] >= 0
 dro_dual_ctr[175]: - oracle[21][5] + 4.833333333333334 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[176]: - oracle[21][6] + 26.7 gamma[21] + nu[21][17] >= 0
 dro_dual_ctr[177]: - oracle[21][7] + 34.76666666666667 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[178]: - oracle[21][8] + 42.06666666666667 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[179]: - oracle[21][9] + 25.5 gamma[21] + nu[21][17] >= 0
 dro_dual_ctr[1710]: - oracle[21][10] + 16.26666666666667 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[1711]: - oracle[21][11] + 49.13333333333334 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[1712]: - oracle[21][12] + 30.86666666666667 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[1713]: - oracle[21][13] + 16.06666666666667 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[1714]: - oracle[21][14] + 9.63333333333334 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[1715]: - oracle[21][15] + 81.76666666666667 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[1716]: - oracle[21][16] + 31.4 gamma[21] + nu[21][17] >= 0
 dro_dual_ctr[1717]: - oracle[21][17] + nu[21][17] >= 0
 dro_dual_ctr[1718]: - oracle[21][18] + 29.96666666666667 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[1719]: - oracle[21][19] + 41.63333333333333 gamma[21]
   + nu[21][17] >= 0
 dro_dual_ctr[180]: - oracle[21][0] + 38.26666666666667 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[181]: - oracle[21][1] + 5.8 gamma[21] + nu[21][18] >= 0
 dro_dual_ctr[182]: - oracle[21][2] + 24.8 gamma[21] + nu[21][18] >= 0
 dro_dual_ctr[183]: - oracle[21][3] + 18.03333333333333 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[184]: - oracle[21][4] + 21.26666666666667 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[185]: - oracle[21][5] + 25.13333333333333 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[186]: - oracle[21][6] + 31.66666666666667 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[187]: - oracle[21][7] + 5.2 gamma[21] + nu[21][18] >= 0
 dro_dual_ctr[188]: - oracle[21][8] + 36.43333333333334 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[189]: - oracle[21][9] + 28.46666666666667 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[1810]: - oracle[21][10] + 19.43333333333333 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[1811]: - oracle[21][11] + 68.3 gamma[21] + nu[21][18] >= 0
 dro_dual_ctr[1812]: - oracle[21][12] + 15.43333333333333 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[1813]: - oracle[21][13] + 16.3 gamma[21] + nu[21][18] >= 0
 dro_dual_ctr[1814]: - oracle[21][14] + 20.33333333333333 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[1815]: - oracle[21][15] + 82.13333333333333 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[1816]: - oracle[21][16] + 46.16666666666666 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[1817]: - oracle[21][17] + 29.96666666666667 gamma[21]
   + nu[21][18] >= 0
 dro_dual_ctr[1818]: - oracle[21][18] + nu[21][18] >= 0
 dro_dual_ctr[1819]: - oracle[21][19] + 37.4 gamma[21] + nu[21][18] >= 0
 dro_dual_ctr[190]: - oracle[21][0] + 21.53333333333333 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[191]: - oracle[21][1] + 43.2 gamma[21] + nu[21][19] >= 0
 dro_dual_ctr[192]: - oracle[21][2] + 51 gamma[21] + nu[21][19] >= 0
 dro_dual_ctr[193]: - oracle[21][3] + 25.16666666666667 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[194]: - oracle[21][4] + 38.13333333333334 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[195]: - oracle[21][5] + 43.46666666666667 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[196]: - oracle[21][6] + 14.93333333333333 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[197]: - oracle[21][7] + 42.2 gamma[21] + nu[21][19] >= 0
 dro_dual_ctr[198]: - oracle[21][8] + 0.9666666666666544 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[199]: - oracle[21][9] + 16.13333333333333 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[1910]: - oracle[21][10] + 25.36666666666667 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[1911]: - oracle[21][11] + 30.9 gamma[21] + nu[21][19] >= 0
 dro_dual_ctr[1912]: - oracle[21][12] + 28.16666666666666 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[1913]: - oracle[21][13] + 42.5 gamma[21] + nu[21][19] >= 0
 dro_dual_ctr[1914]: - oracle[21][14] + 41.66666666666666 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[1915]: - oracle[21][15] + 44.73333333333333 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[1916]: - oracle[21][16] + 10.23333333333333 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[1917]: - oracle[21][17] + 41.63333333333333 gamma[21]
   + nu[21][19] >= 0
 dro_dual_ctr[1918]: - oracle[21][18] + 37.4 gamma[21] + nu[21][19] >= 0
 dro_dual_ctr[1919]: - oracle[21][19] + nu[21][19] >= 0
 cut[21,0,0]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][0] >= 28485.44679282424
 cut[21,0,1]: 2.8530385363795772e+07 P + 0.0001 W + 50.74041793474803 C
   + oracle[21][1] >= 28485.46179288937
 cut[21,0,2]: 2.8530385363795772e+07 P + 0.0001 W + 50.74041793474803 C
   + oracle[21][2] >= 28485.46179288937
 cut[21,0,3]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][3] >= 28485.44679282424
 cut[21,0,4]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][4] >= 28485.44679282424
 cut[21,0,5]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][5] >= 28485.44679282424
 cut[21,0,6]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][6] >= 28485.44679282424
 cut[21,0,7]: 2.8530385363795772e+07 P + 0.0001 W + 50.74041793474803 C
   + oracle[21][7] >= 28485.46179288937
 cut[21,0,8]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][8] >= 28485.44679282424
 cut[21,0,9]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][9] >= 28485.44679282424
 cut[21,0,10]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][10] >= 28485.44679282424
 cut[21,0,11]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][11] >= 28485.44679282424
 cut[21,0,12]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][12] >= 28485.44679282424
 cut[21,0,13]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][13] >= 28485.44679282424
 cut[21,0,14]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][14] >= 28485.44679282424
 cut[21,0,15]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][15] >= 28485.44679282424
 cut[21,0,16]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][16] >= 28485.44679282424
 cut[21,0,17]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][17] >= 28485.44679282424
 cut[21,0,18]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][18] >= 28485.44679282424
 cut[21,0,19]: 2.8530385363802288e+07 P + 50.74041793474803 C
   + oracle[21][19] >= 28485.44679282424
Bounds
 M = -10000
 P0 = 0
 Q0 = 0
 W0 = 150
 C0 = 0
 M0 = -10000
 evapotranspiration = 35.2
 rainfall = 18.53333333333333
 oracle[21][0] >= -1e+08
 oracle[21][1] >= -1e+08
 oracle[21][2] >= -1e+08
 oracle[21][3] >= -1e+08
 oracle[21][4] >= -1e+08
 oracle[21][5] >= -1e+08
 oracle[21][6] >= -1e+08
 oracle[21][7] >= -1e+08
 oracle[21][8] >= -1e+08
 oracle[21][9] >= -1e+08
 oracle[21][10] >= -1e+08
 oracle[21][11] >= -1e+08
 oracle[21][12] >= -1e+08
 oracle[21][13] >= -1e+08
 oracle[21][14] >= -1e+08
 oracle[21][15] >= -1e+08
 oracle[21][16] >= -1e+08
 oracle[21][17] >= -1e+08
 oracle[21][18] >= -1e+08
 oracle[21][19] >= -1e+08
 nu[21][0] free
 nu[21][1] free
 nu[21][2] free
 nu[21][3] free
 nu[21][4] free
 nu[21][5] free
 nu[21][6] free
 nu[21][7] free
 nu[21][8] free
 nu[21][9] free
 nu[21][10] free
 nu[21][11] free
 nu[21][12] free
 nu[21][13] free
 nu[21][14] free
 nu[21][15] free
 nu[21][16] free
 nu[21][17] free
 nu[21][18] free
 nu[21][19] free
End
