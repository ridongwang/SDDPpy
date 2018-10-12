\ Model MilkPowder_22
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.0001 W + 0.5 b + 0.275 h + 0.5 i + C19 + 100 C20 + 10 gamma[22]
   + 0.05 nu[22][0] + 0.05 nu[22][1] + 0.05 nu[22][2] + 0.05 nu[22][3]
   + 0.05 nu[22][4] + 0.05 nu[22][5] + 0.05 nu[22][6] + 0.05 nu[22][7]
   + 0.05 nu[22][8] + 0.05 nu[22][9] + 0.05 nu[22][10] + 0.05 nu[22][11]
   + 0.05 nu[22][12] + 0.05 nu[22][13] + 0.05 nu[22][14] + 0.05 nu[22][15]
   + 0.05 nu[22][16] + 0.05 nu[22][17] + 0.05 nu[22][18] + 0.05 nu[22][19]
Subject To
 Pasture_Cover: P - P0 + h + fp - 7 gr <= 0
 Supplement_Storage: Q - Q0 - 0.9 h + fs <= 0
 Supplement_Storage: Q - Q0 - 0.9 h + fs <= 0
 CowsAmount: C - C0 <= 0
 Soil_Limit: W <= 150
 milkRel: - 0.0123305320298615 mlk + milk <= 0
 milk_BALANCE: M - M0 - milk <= 0
 EnergyBalance: 5.578509914972297 C0 + 11 b + 11 fs + 11 fp - mlk
   >= 1210.478508713518
 maximum_milk: - 1028.696260432342 C0 + mlk <= 0
 minimum_milk: - 6.3 C0 + milk >= 0
 pasture_growth_constraints1: - 1.140581258817171 ev + gr <= 0
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
 dro_dual_ctr[00]: - oracle[22][0] + nu[22][0] >= 0
 dro_dual_ctr[01]: - oracle[22][1] + 46.43333333333334 gamma[22]
   + nu[22][0] >= 0
 dro_dual_ctr[02]: - oracle[22][2] + 25.1 gamma[22] + nu[22][0] >= 0
 dro_dual_ctr[03]: - oracle[22][3] + 43.33333333333334 gamma[22]
   + nu[22][0] >= 0
 dro_dual_ctr[04]: - oracle[22][4] + 35.83333333333333 gamma[22]
   + nu[22][0] >= 0
 dro_dual_ctr[05]: - oracle[22][5] + 32.43333333333334 gamma[22]
   + nu[22][0] >= 0
 dro_dual_ctr[06]: - oracle[22][6] + 129.7666666666667 gamma[22]
   + nu[22][0] >= 0
 dro_dual_ctr[07]: - oracle[22][7] + 44.7 gamma[22] + nu[22][0] >= 0
 dro_dual_ctr[08]: - oracle[22][8] + 34.13333333333333 gamma[22]
   + nu[22][0] >= 0
 dro_dual_ctr[09]: - oracle[22][9] + 45.53333333333333 gamma[22]
   + nu[22][0] >= 0
 dro_dual_ctr[010]: - oracle[22][10] + 61.9 gamma[22] + nu[22][0] >= 0
 dro_dual_ctr[011]: - oracle[22][11] + 41.83333333333334 gamma[22]
   + nu[22][0] >= 0
 dro_dual_ctr[012]: - oracle[22][12] + 32.33333333333333 gamma[22]
   + nu[22][0] >= 0
 dro_dual_ctr[013]: - oracle[22][13] + 26.7 gamma[22] + nu[22][0] >= 0
 dro_dual_ctr[014]: - oracle[22][14] + 47.1 gamma[22] + nu[22][0] >= 0
 dro_dual_ctr[015]: - oracle[22][15] + 7.7 gamma[22] + nu[22][0] >= 0
 dro_dual_ctr[016]: - oracle[22][16] + 50.03333333333333 gamma[22]
   + nu[22][0] >= 0
 dro_dual_ctr[017]: - oracle[22][17] + 37.76666666666667 gamma[22]
   + nu[22][0] >= 0
 dro_dual_ctr[018]: - oracle[22][18] + 37.3 gamma[22] + nu[22][0] >= 0
 dro_dual_ctr[019]: - oracle[22][19] + 37.06666666666668 gamma[22]
   + nu[22][0] >= 0
 dro_dual_ctr[10]: - oracle[22][0] + 46.43333333333334 gamma[22]
   + nu[22][1] >= 0
 dro_dual_ctr[11]: - oracle[22][1] + nu[22][1] >= 0
 dro_dual_ctr[12]: - oracle[22][2] + 71.53333333333333 gamma[22]
   + nu[22][1] >= 0
 dro_dual_ctr[13]: - oracle[22][3] + 4.9 gamma[22] + nu[22][1] >= 0
 dro_dual_ctr[14]: - oracle[22][4] + 10.6 gamma[22] + nu[22][1] >= 0
 dro_dual_ctr[15]: - oracle[22][5] + 78.86666666666667 gamma[22]
   + nu[22][1] >= 0
 dro_dual_ctr[16]: - oracle[22][6] + 176.2 gamma[22] + nu[22][1] >= 0
 dro_dual_ctr[17]: - oracle[22][7] + 1.733333333333335 gamma[22]
   + nu[22][1] >= 0
 dro_dual_ctr[18]: - oracle[22][8] + 13.5 gamma[22] + nu[22][1] >= 0
 dro_dual_ctr[19]: - oracle[22][9] + 2.7 gamma[22] + nu[22][1] >= 0
 dro_dual_ctr[110]: - oracle[22][10] + 108.3333333333333 gamma[22]
   + nu[22][1] >= 0
 dro_dual_ctr[111]: - oracle[22][11] + 4.6 gamma[22] + nu[22][1] >= 0
 dro_dual_ctr[112]: - oracle[22][12] + 17.7 gamma[22] + nu[22][1] >= 0
 dro_dual_ctr[113]: - oracle[22][13] + 19.73333333333334 gamma[22]
   + nu[22][1] >= 0
 dro_dual_ctr[114]: - oracle[22][14] + 3.133333333333332 gamma[22]
   + nu[22][1] >= 0
 dro_dual_ctr[115]: - oracle[22][15] + 39.93333333333334 gamma[22]
   + nu[22][1] >= 0
 dro_dual_ctr[116]: - oracle[22][16] + 3.6 gamma[22] + nu[22][1] >= 0
 dro_dual_ctr[117]: - oracle[22][17] + 8.666666666666671 gamma[22]
   + nu[22][1] >= 0
 dro_dual_ctr[118]: - oracle[22][18] + 13.73333333333334 gamma[22]
   + nu[22][1] >= 0
 dro_dual_ctr[119]: - oracle[22][19] + 9.366666666666665 gamma[22]
   + nu[22][1] >= 0
 dro_dual_ctr[20]: - oracle[22][0] + 25.1 gamma[22] + nu[22][2] >= 0
 dro_dual_ctr[21]: - oracle[22][1] + 71.53333333333333 gamma[22]
   + nu[22][2] >= 0
 dro_dual_ctr[22]: - oracle[22][2] + nu[22][2] >= 0
 dro_dual_ctr[23]: - oracle[22][3] + 68.43333333333334 gamma[22]
   + nu[22][2] >= 0
 dro_dual_ctr[24]: - oracle[22][4] + 60.93333333333333 gamma[22]
   + nu[22][2] >= 0
 dro_dual_ctr[25]: - oracle[22][5] + 7.333333333333336 gamma[22]
   + nu[22][2] >= 0
 dro_dual_ctr[26]: - oracle[22][6] + 104.6666666666667 gamma[22]
   + nu[22][2] >= 0
 dro_dual_ctr[27]: - oracle[22][7] + 69.8 gamma[22] + nu[22][2] >= 0
 dro_dual_ctr[28]: - oracle[22][8] + 59.23333333333333 gamma[22]
   + nu[22][2] >= 0
 dro_dual_ctr[29]: - oracle[22][9] + 70.63333333333333 gamma[22]
   + nu[22][2] >= 0
 dro_dual_ctr[210]: - oracle[22][10] + 36.8 gamma[22] + nu[22][2] >= 0
 dro_dual_ctr[211]: - oracle[22][11] + 66.93333333333334 gamma[22]
   + nu[22][2] >= 0
 dro_dual_ctr[212]: - oracle[22][12] + 57.43333333333333 gamma[22]
   + nu[22][2] >= 0
 dro_dual_ctr[213]: - oracle[22][13] + 51.8 gamma[22] + nu[22][2] >= 0
 dro_dual_ctr[214]: - oracle[22][14] + 72.2 gamma[22] + nu[22][2] >= 0
 dro_dual_ctr[215]: - oracle[22][15] + 31.6 gamma[22] + nu[22][2] >= 0
 dro_dual_ctr[216]: - oracle[22][16] + 75.13333333333333 gamma[22]
   + nu[22][2] >= 0
 dro_dual_ctr[217]: - oracle[22][17] + 62.86666666666667 gamma[22]
   + nu[22][2] >= 0
 dro_dual_ctr[218]: - oracle[22][18] + 62.4 gamma[22] + nu[22][2] >= 0
 dro_dual_ctr[219]: - oracle[22][19] + 62.16666666666667 gamma[22]
   + nu[22][2] >= 0
 dro_dual_ctr[30]: - oracle[22][0] + 43.33333333333334 gamma[22]
   + nu[22][3] >= 0
 dro_dual_ctr[31]: - oracle[22][1] + 4.9 gamma[22] + nu[22][3] >= 0
 dro_dual_ctr[32]: - oracle[22][2] + 68.43333333333334 gamma[22]
   + nu[22][3] >= 0
 dro_dual_ctr[33]: - oracle[22][3] + nu[22][3] >= 0
 dro_dual_ctr[34]: - oracle[22][4] + 7.5 gamma[22] + nu[22][3] >= 0
 dro_dual_ctr[35]: - oracle[22][5] + 75.76666666666667 gamma[22]
   + nu[22][3] >= 0
 dro_dual_ctr[36]: - oracle[22][6] + 173.1 gamma[22] + nu[22][3] >= 0
 dro_dual_ctr[37]: - oracle[22][7] + 6.433333333333332 gamma[22]
   + nu[22][3] >= 0
 dro_dual_ctr[38]: - oracle[22][8] + 18.4 gamma[22] + nu[22][3] >= 0
 dro_dual_ctr[39]: - oracle[22][9] + 2.2 gamma[22] + nu[22][3] >= 0
 dro_dual_ctr[310]: - oracle[22][10] + 105.2333333333333 gamma[22]
   + nu[22][3] >= 0
 dro_dual_ctr[311]: - oracle[22][11] + 2.7 gamma[22] + nu[22][3] >= 0
 dro_dual_ctr[312]: - oracle[22][12] + 22.6 gamma[22] + nu[22][3] >= 0
 dro_dual_ctr[313]: - oracle[22][13] + 18.03333333333333 gamma[22]
   + nu[22][3] >= 0
 dro_dual_ctr[314]: - oracle[22][14] + 8.033333333333331 gamma[22]
   + nu[22][3] >= 0
 dro_dual_ctr[315]: - oracle[22][15] + 39.03333333333335 gamma[22]
   + nu[22][3] >= 0
 dro_dual_ctr[316]: - oracle[22][16] + 6.7 gamma[22] + nu[22][3] >= 0
 dro_dual_ctr[317]: - oracle[22][17] + 12.96666666666666 gamma[22]
   + nu[22][3] >= 0
 dro_dual_ctr[318]: - oracle[22][18] + 18.63333333333334 gamma[22]
   + nu[22][3] >= 0
 dro_dual_ctr[319]: - oracle[22][19] + 7.466666666666668 gamma[22]
   + nu[22][3] >= 0
 dro_dual_ctr[40]: - oracle[22][0] + 35.83333333333333 gamma[22]
   + nu[22][4] >= 0
 dro_dual_ctr[41]: - oracle[22][1] + 10.6 gamma[22] + nu[22][4] >= 0
 dro_dual_ctr[42]: - oracle[22][2] + 60.93333333333333 gamma[22]
   + nu[22][4] >= 0
 dro_dual_ctr[43]: - oracle[22][3] + 7.5 gamma[22] + nu[22][4] >= 0
 dro_dual_ctr[44]: - oracle[22][4] + nu[22][4] >= 0
 dro_dual_ctr[45]: - oracle[22][5] + 68.26666666666667 gamma[22]
   + nu[22][4] >= 0
 dro_dual_ctr[46]: - oracle[22][6] + 165.6 gamma[22] + nu[22][4] >= 0
 dro_dual_ctr[47]: - oracle[22][7] + 8.866666666666671 gamma[22]
   + nu[22][4] >= 0
 dro_dual_ctr[48]: - oracle[22][8] + 11.7 gamma[22] + nu[22][4] >= 0
 dro_dual_ctr[49]: - oracle[22][9] + 9.7 gamma[22] + nu[22][4] >= 0
 dro_dual_ctr[410]: - oracle[22][10] + 97.73333333333332 gamma[22]
   + nu[22][4] >= 0
 dro_dual_ctr[411]: - oracle[22][11] + 6 gamma[22] + nu[22][4] >= 0
 dro_dual_ctr[412]: - oracle[22][12] + 15.9 gamma[22] + nu[22][4] >= 0
 dro_dual_ctr[413]: - oracle[22][13] + 11.33333333333334 gamma[22]
   + nu[22][4] >= 0
 dro_dual_ctr[414]: - oracle[22][14] + 11.26666666666667 gamma[22]
   + nu[22][4] >= 0
 dro_dual_ctr[415]: - oracle[22][15] + 32.33333333333335 gamma[22]
   + nu[22][4] >= 0
 dro_dual_ctr[416]: - oracle[22][16] + 14.2 gamma[22] + nu[22][4] >= 0
 dro_dual_ctr[417]: - oracle[22][17] + 6.266666666666669 gamma[22]
   + nu[22][4] >= 0
 dro_dual_ctr[418]: - oracle[22][18] + 11.93333333333334 gamma[22]
   + nu[22][4] >= 0
 dro_dual_ctr[419]: - oracle[22][19] + 1.23333333333334 gamma[22]
   + nu[22][4] >= 0
 dro_dual_ctr[50]: - oracle[22][0] + 32.43333333333334 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[51]: - oracle[22][1] + 78.86666666666667 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[52]: - oracle[22][2] + 7.333333333333336 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[53]: - oracle[22][3] + 75.76666666666667 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[54]: - oracle[22][4] + 68.26666666666667 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[55]: - oracle[22][5] + nu[22][5] >= 0
 dro_dual_ctr[56]: - oracle[22][6] + 97.33333333333334 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[57]: - oracle[22][7] + 77.13333333333333 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[58]: - oracle[22][8] + 66.56666666666666 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[59]: - oracle[22][9] + 77.96666666666667 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[510]: - oracle[22][10] + 38.26666666666667 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[511]: - oracle[22][11] + 74.26666666666668 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[512]: - oracle[22][12] + 64.76666666666667 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[513]: - oracle[22][13] + 59.13333333333333 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[514]: - oracle[22][14] + 79.53333333333333 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[515]: - oracle[22][15] + 38.93333333333333 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[516]: - oracle[22][16] + 82.46666666666667 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[517]: - oracle[22][17] + 70.2 gamma[22] + nu[22][5] >= 0
 dro_dual_ctr[518]: - oracle[22][18] + 69.73333333333335 gamma[22]
   + nu[22][5] >= 0
 dro_dual_ctr[519]: - oracle[22][19] + 69.5 gamma[22] + nu[22][5] >= 0
 dro_dual_ctr[60]: - oracle[22][0] + 129.7666666666667 gamma[22]
   + nu[22][6] >= 0
 dro_dual_ctr[61]: - oracle[22][1] + 176.2 gamma[22] + nu[22][6] >= 0
 dro_dual_ctr[62]: - oracle[22][2] + 104.6666666666667 gamma[22]
   + nu[22][6] >= 0
 dro_dual_ctr[63]: - oracle[22][3] + 173.1 gamma[22] + nu[22][6] >= 0
 dro_dual_ctr[64]: - oracle[22][4] + 165.6 gamma[22] + nu[22][6] >= 0
 dro_dual_ctr[65]: - oracle[22][5] + 97.33333333333334 gamma[22]
   + nu[22][6] >= 0
 dro_dual_ctr[66]: - oracle[22][6] + nu[22][6] >= 0
 dro_dual_ctr[67]: - oracle[22][7] + 174.4666666666667 gamma[22]
   + nu[22][6] >= 0
 dro_dual_ctr[68]: - oracle[22][8] + 163.9 gamma[22] + nu[22][6] >= 0
 dro_dual_ctr[69]: - oracle[22][9] + 175.3 gamma[22] + nu[22][6] >= 0
 dro_dual_ctr[610]: - oracle[22][10] + 67.86666666666667 gamma[22]
   + nu[22][6] >= 0
 dro_dual_ctr[611]: - oracle[22][11] + 171.6 gamma[22] + nu[22][6] >= 0
 dro_dual_ctr[612]: - oracle[22][12] + 162.1 gamma[22] + nu[22][6] >= 0
 dro_dual_ctr[613]: - oracle[22][13] + 156.4666666666667 gamma[22]
   + nu[22][6] >= 0
 dro_dual_ctr[614]: - oracle[22][14] + 176.8666666666667 gamma[22]
   + nu[22][6] >= 0
 dro_dual_ctr[615]: - oracle[22][15] + 136.2666666666667 gamma[22]
   + nu[22][6] >= 0
 dro_dual_ctr[616]: - oracle[22][16] + 179.8 gamma[22] + nu[22][6] >= 0
 dro_dual_ctr[617]: - oracle[22][17] + 167.5333333333333 gamma[22]
   + nu[22][6] >= 0
 dro_dual_ctr[618]: - oracle[22][18] + 167.0666666666667 gamma[22]
   + nu[22][6] >= 0
 dro_dual_ctr[619]: - oracle[22][19] + 166.8333333333333 gamma[22]
   + nu[22][6] >= 0
 dro_dual_ctr[70]: - oracle[22][0] + 44.7 gamma[22] + nu[22][7] >= 0
 dro_dual_ctr[71]: - oracle[22][1] + 1.733333333333335 gamma[22]
   + nu[22][7] >= 0
 dro_dual_ctr[72]: - oracle[22][2] + 69.8 gamma[22] + nu[22][7] >= 0
 dro_dual_ctr[73]: - oracle[22][3] + 6.433333333333332 gamma[22]
   + nu[22][7] >= 0
 dro_dual_ctr[74]: - oracle[22][4] + 8.866666666666671 gamma[22]
   + nu[22][7] >= 0
 dro_dual_ctr[75]: - oracle[22][5] + 77.13333333333333 gamma[22]
   + nu[22][7] >= 0
 dro_dual_ctr[76]: - oracle[22][6] + 174.4666666666667 gamma[22]
   + nu[22][7] >= 0
 dro_dual_ctr[77]: - oracle[22][7] + nu[22][7] >= 0
 dro_dual_ctr[78]: - oracle[22][8] + 11.96666666666666 gamma[22]
   + nu[22][7] >= 0
 dro_dual_ctr[79]: - oracle[22][9] + 4.233333333333336 gamma[22]
   + nu[22][7] >= 0
 dro_dual_ctr[710]: - oracle[22][10] + 106.6 gamma[22] + nu[22][7] >= 0
 dro_dual_ctr[711]: - oracle[22][11] + 3.73333333333333 gamma[22]
   + nu[22][7] >= 0
 dro_dual_ctr[712]: - oracle[22][12] + 16.16666666666667 gamma[22]
   + nu[22][7] >= 0
 dro_dual_ctr[713]: - oracle[22][13] + 18 gamma[22] + nu[22][7] >= 0
 dro_dual_ctr[714]: - oracle[22][14] + 2.4 gamma[22] + nu[22][7] >= 0
 dro_dual_ctr[715]: - oracle[22][15] + 38.2 gamma[22] + nu[22][7] >= 0
 dro_dual_ctr[716]: - oracle[22][16] + 5.333333333333336 gamma[22]
   + nu[22][7] >= 0
 dro_dual_ctr[717]: - oracle[22][17] + 6.933333333333336 gamma[22]
   + nu[22][7] >= 0
 dro_dual_ctr[718]: - oracle[22][18] + 12.2 gamma[22] + nu[22][7] >= 0
 dro_dual_ctr[719]: - oracle[22][19] + 7.63333333333333 gamma[22]
   + nu[22][7] >= 0
 dro_dual_ctr[80]: - oracle[22][0] + 34.13333333333333 gamma[22]
   + nu[22][8] >= 0
 dro_dual_ctr[81]: - oracle[22][1] + 13.5 gamma[22] + nu[22][8] >= 0
 dro_dual_ctr[82]: - oracle[22][2] + 59.23333333333333 gamma[22]
   + nu[22][8] >= 0
 dro_dual_ctr[83]: - oracle[22][3] + 18.4 gamma[22] + nu[22][8] >= 0
 dro_dual_ctr[84]: - oracle[22][4] + 11.7 gamma[22] + nu[22][8] >= 0
 dro_dual_ctr[85]: - oracle[22][5] + 66.56666666666666 gamma[22]
   + nu[22][8] >= 0
 dro_dual_ctr[86]: - oracle[22][6] + 163.9 gamma[22] + nu[22][8] >= 0
 dro_dual_ctr[87]: - oracle[22][7] + 11.96666666666666 gamma[22]
   + nu[22][8] >= 0
 dro_dual_ctr[88]: - oracle[22][8] + nu[22][8] >= 0
 dro_dual_ctr[89]: - oracle[22][9] + 16.2 gamma[22] + nu[22][8] >= 0
 dro_dual_ctr[810]: - oracle[22][10] + 96.03333333333332 gamma[22]
   + nu[22][8] >= 0
 dro_dual_ctr[811]: - oracle[22][11] + 15.7 gamma[22] + nu[22][8] >= 0
 dro_dual_ctr[812]: - oracle[22][12] + 4.2 gamma[22] + nu[22][8] >= 0
 dro_dual_ctr[813]: - oracle[22][13] + 7.433333333333334 gamma[22]
   + nu[22][8] >= 0
 dro_dual_ctr[814]: - oracle[22][14] + 12.96666666666667 gamma[22]
   + nu[22][8] >= 0
 dro_dual_ctr[815]: - oracle[22][15] + 27.63333333333334 gamma[22]
   + nu[22][8] >= 0
 dro_dual_ctr[816]: - oracle[22][16] + 15.9 gamma[22] + nu[22][8] >= 0
 dro_dual_ctr[817]: - oracle[22][17] + 5.43333333333333 gamma[22]
   + nu[22][8] >= 0
 dro_dual_ctr[818]: - oracle[22][18] + 3.166666666666675 gamma[22]
   + nu[22][8] >= 0
 dro_dual_ctr[819]: - oracle[22][19] + 10.93333333333332 gamma[22]
   + nu[22][8] >= 0
 dro_dual_ctr[90]: - oracle[22][0] + 45.53333333333333 gamma[22]
   + nu[22][9] >= 0
 dro_dual_ctr[91]: - oracle[22][1] + 2.7 gamma[22] + nu[22][9] >= 0
 dro_dual_ctr[92]: - oracle[22][2] + 70.63333333333333 gamma[22]
   + nu[22][9] >= 0
 dro_dual_ctr[93]: - oracle[22][3] + 2.2 gamma[22] + nu[22][9] >= 0
 dro_dual_ctr[94]: - oracle[22][4] + 9.7 gamma[22] + nu[22][9] >= 0
 dro_dual_ctr[95]: - oracle[22][5] + 77.96666666666667 gamma[22]
   + nu[22][9] >= 0
 dro_dual_ctr[96]: - oracle[22][6] + 175.3 gamma[22] + nu[22][9] >= 0
 dro_dual_ctr[97]: - oracle[22][7] + 4.233333333333336 gamma[22]
   + nu[22][9] >= 0
 dro_dual_ctr[98]: - oracle[22][8] + 16.2 gamma[22] + nu[22][9] >= 0
 dro_dual_ctr[99]: - oracle[22][9] + nu[22][9] >= 0
 dro_dual_ctr[910]: - oracle[22][10] + 107.4333333333333 gamma[22]
   + nu[22][9] >= 0
 dro_dual_ctr[911]: - oracle[22][11] + 3.7 gamma[22] + nu[22][9] >= 0
 dro_dual_ctr[912]: - oracle[22][12] + 20.4 gamma[22] + nu[22][9] >= 0
 dro_dual_ctr[913]: - oracle[22][13] + 18.83333333333333 gamma[22]
   + nu[22][9] >= 0
 dro_dual_ctr[914]: - oracle[22][14] + 5.833333333333337 gamma[22]
   + nu[22][9] >= 0
 dro_dual_ctr[915]: - oracle[22][15] + 39.03333333333334 gamma[22]
   + nu[22][9] >= 0
 dro_dual_ctr[916]: - oracle[22][16] + 4.5 gamma[22] + nu[22][9] >= 0
 dro_dual_ctr[917]: - oracle[22][17] + 10.76666666666667 gamma[22]
   + nu[22][9] >= 0
 dro_dual_ctr[918]: - oracle[22][18] + 16.43333333333334 gamma[22]
   + nu[22][9] >= 0
 dro_dual_ctr[919]: - oracle[22][19] + 8.466666666666661 gamma[22]
   + nu[22][9] >= 0
 dro_dual_ctr[100]: - oracle[22][0] + 61.9 gamma[22] + nu[22][10] >= 0
 dro_dual_ctr[101]: - oracle[22][1] + 108.3333333333333 gamma[22]
   + nu[22][10] >= 0
 dro_dual_ctr[102]: - oracle[22][2] + 36.8 gamma[22] + nu[22][10] >= 0
 dro_dual_ctr[103]: - oracle[22][3] + 105.2333333333333 gamma[22]
   + nu[22][10] >= 0
 dro_dual_ctr[104]: - oracle[22][4] + 97.73333333333332 gamma[22]
   + nu[22][10] >= 0
 dro_dual_ctr[105]: - oracle[22][5] + 38.26666666666667 gamma[22]
   + nu[22][10] >= 0
 dro_dual_ctr[106]: - oracle[22][6] + 67.86666666666667 gamma[22]
   + nu[22][10] >= 0
 dro_dual_ctr[107]: - oracle[22][7] + 106.6 gamma[22] + nu[22][10] >= 0
 dro_dual_ctr[108]: - oracle[22][8] + 96.03333333333332 gamma[22]
   + nu[22][10] >= 0
 dro_dual_ctr[109]: - oracle[22][9] + 107.4333333333333 gamma[22]
   + nu[22][10] >= 0
 dro_dual_ctr[1010]: - oracle[22][10] + nu[22][10] >= 0
 dro_dual_ctr[1011]: - oracle[22][11] + 103.7333333333333 gamma[22]
   + nu[22][10] >= 0
 dro_dual_ctr[1012]: - oracle[22][12] + 94.23333333333332 gamma[22]
   + nu[22][10] >= 0
 dro_dual_ctr[1013]: - oracle[22][13] + 88.6 gamma[22] + nu[22][10] >= 0
 dro_dual_ctr[1014]: - oracle[22][14] + 109 gamma[22] + nu[22][10] >= 0
 dro_dual_ctr[1015]: - oracle[22][15] + 68.4 gamma[22] + nu[22][10] >= 0
 dro_dual_ctr[1016]: - oracle[22][16] + 111.9333333333333 gamma[22]
   + nu[22][10] >= 0
 dro_dual_ctr[1017]: - oracle[22][17] + 99.66666666666666 gamma[22]
   + nu[22][10] >= 0
 dro_dual_ctr[1018]: - oracle[22][18] + 99.2 gamma[22] + nu[22][10] >= 0
 dro_dual_ctr[1019]: - oracle[22][19] + 98.96666666666667 gamma[22]
   + nu[22][10] >= 0
 dro_dual_ctr[110]: - oracle[22][0] + 41.83333333333334 gamma[22]
   + nu[22][11] >= 0
 dro_dual_ctr[111]: - oracle[22][1] + 4.6 gamma[22] + nu[22][11] >= 0
 dro_dual_ctr[112]: - oracle[22][2] + 66.93333333333334 gamma[22]
   + nu[22][11] >= 0
 dro_dual_ctr[113]: - oracle[22][3] + 2.7 gamma[22] + nu[22][11] >= 0
 dro_dual_ctr[114]: - oracle[22][4] + 6 gamma[22] + nu[22][11] >= 0
 dro_dual_ctr[115]: - oracle[22][5] + 74.26666666666668 gamma[22]
   + nu[22][11] >= 0
 dro_dual_ctr[116]: - oracle[22][6] + 171.6 gamma[22] + nu[22][11] >= 0
 dro_dual_ctr[117]: - oracle[22][7] + 3.73333333333333 gamma[22]
   + nu[22][11] >= 0
 dro_dual_ctr[118]: - oracle[22][8] + 15.7 gamma[22] + nu[22][11] >= 0
 dro_dual_ctr[119]: - oracle[22][9] + 3.7 gamma[22] + nu[22][11] >= 0
 dro_dual_ctr[1110]: - oracle[22][10] + 103.7333333333333 gamma[22]
   + nu[22][11] >= 0
 dro_dual_ctr[1111]: - oracle[22][11] + nu[22][11] >= 0
 dro_dual_ctr[1112]: - oracle[22][12] + 19.9 gamma[22] + nu[22][11] >= 0
 dro_dual_ctr[1113]: - oracle[22][13] + 15.33333333333333 gamma[22]
   + nu[22][11] >= 0
 dro_dual_ctr[1114]: - oracle[22][14] + 5.33333333333333 gamma[22]
   + nu[22][11] >= 0
 dro_dual_ctr[1115]: - oracle[22][15] + 36.33333333333334 gamma[22]
   + nu[22][11] >= 0
 dro_dual_ctr[1116]: - oracle[22][16] + 8.2 gamma[22] + nu[22][11] >= 0
 dro_dual_ctr[1117]: - oracle[22][17] + 10.26666666666666 gamma[22]
   + nu[22][11] >= 0
 dro_dual_ctr[1118]: - oracle[22][18] + 15.93333333333334 gamma[22]
   + nu[22][11] >= 0
 dro_dual_ctr[1119]: - oracle[22][19] + 4.766666666666666 gamma[22]
   + nu[22][11] >= 0
 dro_dual_ctr[120]: - oracle[22][0] + 32.33333333333333 gamma[22]
   + nu[22][12] >= 0
 dro_dual_ctr[121]: - oracle[22][1] + 17.7 gamma[22] + nu[22][12] >= 0
 dro_dual_ctr[122]: - oracle[22][2] + 57.43333333333333 gamma[22]
   + nu[22][12] >= 0
 dro_dual_ctr[123]: - oracle[22][3] + 22.6 gamma[22] + nu[22][12] >= 0
 dro_dual_ctr[124]: - oracle[22][4] + 15.9 gamma[22] + nu[22][12] >= 0
 dro_dual_ctr[125]: - oracle[22][5] + 64.76666666666667 gamma[22]
   + nu[22][12] >= 0
 dro_dual_ctr[126]: - oracle[22][6] + 162.1 gamma[22] + nu[22][12] >= 0
 dro_dual_ctr[127]: - oracle[22][7] + 16.16666666666667 gamma[22]
   + nu[22][12] >= 0
 dro_dual_ctr[128]: - oracle[22][8] + 4.2 gamma[22] + nu[22][12] >= 0
 dro_dual_ctr[129]: - oracle[22][9] + 20.4 gamma[22] + nu[22][12] >= 0
 dro_dual_ctr[1210]: - oracle[22][10] + 94.23333333333332 gamma[22]
   + nu[22][12] >= 0
 dro_dual_ctr[1211]: - oracle[22][11] + 19.9 gamma[22] + nu[22][12] >= 0
 dro_dual_ctr[1212]: - oracle[22][12] + nu[22][12] >= 0
 dro_dual_ctr[1213]: - oracle[22][13] + 5.633333333333333 gamma[22]
   + nu[22][12] >= 0
 dro_dual_ctr[1214]: - oracle[22][14] + 14.76666666666667 gamma[22]
   + nu[22][12] >= 0
 dro_dual_ctr[1215]: - oracle[22][15] + 25.83333333333334 gamma[22]
   + nu[22][12] >= 0
 dro_dual_ctr[1216]: - oracle[22][16] + 17.7 gamma[22] + nu[22][12] >= 0
 dro_dual_ctr[1217]: - oracle[22][17] + 9.633333333333336 gamma[22]
   + nu[22][12] >= 0
 dro_dual_ctr[1218]: - oracle[22][18] + 4.966666666666676 gamma[22]
   + nu[22][12] >= 0
 dro_dual_ctr[1219]: - oracle[22][19] + 15.13333333333333 gamma[22]
   + nu[22][12] >= 0
 dro_dual_ctr[130]: - oracle[22][0] + 26.7 gamma[22] + nu[22][13] >= 0
 dro_dual_ctr[131]: - oracle[22][1] + 19.73333333333334 gamma[22]
   + nu[22][13] >= 0
 dro_dual_ctr[132]: - oracle[22][2] + 51.8 gamma[22] + nu[22][13] >= 0
 dro_dual_ctr[133]: - oracle[22][3] + 18.03333333333333 gamma[22]
   + nu[22][13] >= 0
 dro_dual_ctr[134]: - oracle[22][4] + 11.33333333333334 gamma[22]
   + nu[22][13] >= 0
 dro_dual_ctr[135]: - oracle[22][5] + 59.13333333333333 gamma[22]
   + nu[22][13] >= 0
 dro_dual_ctr[136]: - oracle[22][6] + 156.4666666666667 gamma[22]
   + nu[22][13] >= 0
 dro_dual_ctr[137]: - oracle[22][7] + 18 gamma[22] + nu[22][13] >= 0
 dro_dual_ctr[138]: - oracle[22][8] + 7.433333333333334 gamma[22]
   + nu[22][13] >= 0
 dro_dual_ctr[139]: - oracle[22][9] + 18.83333333333333 gamma[22]
   + nu[22][13] >= 0
 dro_dual_ctr[1310]: - oracle[22][10] + 88.6 gamma[22] + nu[22][13] >= 0
 dro_dual_ctr[1311]: - oracle[22][11] + 15.33333333333333 gamma[22]
   + nu[22][13] >= 0
 dro_dual_ctr[1312]: - oracle[22][12] + 5.633333333333333 gamma[22]
   + nu[22][13] >= 0
 dro_dual_ctr[1313]: - oracle[22][13] + nu[22][13] >= 0
 dro_dual_ctr[1314]: - oracle[22][14] + 20.4 gamma[22] + nu[22][13] >= 0
 dro_dual_ctr[1315]: - oracle[22][15] + 21 gamma[22] + nu[22][13] >= 0
 dro_dual_ctr[1316]: - oracle[22][16] + 23.33333333333334 gamma[22]
   + nu[22][13] >= 0
 dro_dual_ctr[1317]: - oracle[22][17] + 11.06666666666667 gamma[22]
   + nu[22][13] >= 0
 dro_dual_ctr[1318]: - oracle[22][18] + 10.6 gamma[22] + nu[22][13] >= 0
 dro_dual_ctr[1319]: - oracle[22][19] + 10.56666666666666 gamma[22]
   + nu[22][13] >= 0
 dro_dual_ctr[140]: - oracle[22][0] + 47.1 gamma[22] + nu[22][14] >= 0
 dro_dual_ctr[141]: - oracle[22][1] + 3.133333333333332 gamma[22]
   + nu[22][14] >= 0
 dro_dual_ctr[142]: - oracle[22][2] + 72.2 gamma[22] + nu[22][14] >= 0
 dro_dual_ctr[143]: - oracle[22][3] + 8.033333333333331 gamma[22]
   + nu[22][14] >= 0
 dro_dual_ctr[144]: - oracle[22][4] + 11.26666666666667 gamma[22]
   + nu[22][14] >= 0
 dro_dual_ctr[145]: - oracle[22][5] + 79.53333333333333 gamma[22]
   + nu[22][14] >= 0
 dro_dual_ctr[146]: - oracle[22][6] + 176.8666666666667 gamma[22]
   + nu[22][14] >= 0
 dro_dual_ctr[147]: - oracle[22][7] + 2.4 gamma[22] + nu[22][14] >= 0
 dro_dual_ctr[148]: - oracle[22][8] + 12.96666666666667 gamma[22]
   + nu[22][14] >= 0
 dro_dual_ctr[149]: - oracle[22][9] + 5.833333333333337 gamma[22]
   + nu[22][14] >= 0
 dro_dual_ctr[1410]: - oracle[22][10] + 109 gamma[22] + nu[22][14] >= 0
 dro_dual_ctr[1411]: - oracle[22][11] + 5.33333333333333 gamma[22]
   + nu[22][14] >= 0
 dro_dual_ctr[1412]: - oracle[22][12] + 14.76666666666667 gamma[22]
   + nu[22][14] >= 0
 dro_dual_ctr[1413]: - oracle[22][13] + 20.4 gamma[22] + nu[22][14] >= 0
 dro_dual_ctr[1414]: - oracle[22][14] + nu[22][14] >= 0
 dro_dual_ctr[1415]: - oracle[22][15] + 40.6 gamma[22] + nu[22][14] >= 0
 dro_dual_ctr[1416]: - oracle[22][16] + 2.933333333333336 gamma[22]
   + nu[22][14] >= 0
 dro_dual_ctr[1417]: - oracle[22][17] + 9.333333333333336 gamma[22]
   + nu[22][14] >= 0
 dro_dual_ctr[1418]: - oracle[22][18] + 10.6 gamma[22] + nu[22][14] >= 0
 dro_dual_ctr[1419]: - oracle[22][19] + 10.03333333333333 gamma[22]
   + nu[22][14] >= 0
 dro_dual_ctr[150]: - oracle[22][0] + 7.7 gamma[22] + nu[22][15] >= 0
 dro_dual_ctr[151]: - oracle[22][1] + 39.93333333333334 gamma[22]
   + nu[22][15] >= 0
 dro_dual_ctr[152]: - oracle[22][2] + 31.6 gamma[22] + nu[22][15] >= 0
 dro_dual_ctr[153]: - oracle[22][3] + 39.03333333333335 gamma[22]
   + nu[22][15] >= 0
 dro_dual_ctr[154]: - oracle[22][4] + 32.33333333333335 gamma[22]
   + nu[22][15] >= 0
 dro_dual_ctr[155]: - oracle[22][5] + 38.93333333333333 gamma[22]
   + nu[22][15] >= 0
 dro_dual_ctr[156]: - oracle[22][6] + 136.2666666666667 gamma[22]
   + nu[22][15] >= 0
 dro_dual_ctr[157]: - oracle[22][7] + 38.2 gamma[22] + nu[22][15] >= 0
 dro_dual_ctr[158]: - oracle[22][8] + 27.63333333333334 gamma[22]
   + nu[22][15] >= 0
 dro_dual_ctr[159]: - oracle[22][9] + 39.03333333333334 gamma[22]
   + nu[22][15] >= 0
 dro_dual_ctr[1510]: - oracle[22][10] + 68.4 gamma[22] + nu[22][15] >= 0
 dro_dual_ctr[1511]: - oracle[22][11] + 36.33333333333334 gamma[22]
   + nu[22][15] >= 0
 dro_dual_ctr[1512]: - oracle[22][12] + 25.83333333333334 gamma[22]
   + nu[22][15] >= 0
 dro_dual_ctr[1513]: - oracle[22][13] + 21 gamma[22] + nu[22][15] >= 0
 dro_dual_ctr[1514]: - oracle[22][14] + 40.6 gamma[22] + nu[22][15] >= 0
 dro_dual_ctr[1515]: - oracle[22][15] + nu[22][15] >= 0
 dro_dual_ctr[1516]: - oracle[22][16] + 43.53333333333334 gamma[22]
   + nu[22][15] >= 0
 dro_dual_ctr[1517]: - oracle[22][17] + 31.26666666666667 gamma[22]
   + nu[22][15] >= 0
 dro_dual_ctr[1518]: - oracle[22][18] + 30.8 gamma[22] + nu[22][15] >= 0
 dro_dual_ctr[1519]: - oracle[22][19] + 31.56666666666668 gamma[22]
   + nu[22][15] >= 0
 dro_dual_ctr[160]: - oracle[22][0] + 50.03333333333333 gamma[22]
   + nu[22][16] >= 0
 dro_dual_ctr[161]: - oracle[22][1] + 3.6 gamma[22] + nu[22][16] >= 0
 dro_dual_ctr[162]: - oracle[22][2] + 75.13333333333333 gamma[22]
   + nu[22][16] >= 0
 dro_dual_ctr[163]: - oracle[22][3] + 6.7 gamma[22] + nu[22][16] >= 0
 dro_dual_ctr[164]: - oracle[22][4] + 14.2 gamma[22] + nu[22][16] >= 0
 dro_dual_ctr[165]: - oracle[22][5] + 82.46666666666667 gamma[22]
   + nu[22][16] >= 0
 dro_dual_ctr[166]: - oracle[22][6] + 179.8 gamma[22] + nu[22][16] >= 0
 dro_dual_ctr[167]: - oracle[22][7] + 5.333333333333336 gamma[22]
   + nu[22][16] >= 0
 dro_dual_ctr[168]: - oracle[22][8] + 15.9 gamma[22] + nu[22][16] >= 0
 dro_dual_ctr[169]: - oracle[22][9] + 4.5 gamma[22] + nu[22][16] >= 0
 dro_dual_ctr[1610]: - oracle[22][10] + 111.9333333333333 gamma[22]
   + nu[22][16] >= 0
 dro_dual_ctr[1611]: - oracle[22][11] + 8.2 gamma[22] + nu[22][16] >= 0
 dro_dual_ctr[1612]: - oracle[22][12] + 17.7 gamma[22] + nu[22][16] >= 0
 dro_dual_ctr[1613]: - oracle[22][13] + 23.33333333333334 gamma[22]
   + nu[22][16] >= 0
 dro_dual_ctr[1614]: - oracle[22][14] + 2.933333333333336 gamma[22]
   + nu[22][16] >= 0
 dro_dual_ctr[1615]: - oracle[22][15] + 43.53333333333334 gamma[22]
   + nu[22][16] >= 0
 dro_dual_ctr[1616]: - oracle[22][16] + nu[22][16] >= 0
 dro_dual_ctr[1617]: - oracle[22][17] + 12.26666666666667 gamma[22]
   + nu[22][16] >= 0
 dro_dual_ctr[1618]: - oracle[22][18] + 13.13333333333334 gamma[22]
   + nu[22][16] >= 0
 dro_dual_ctr[1619]: - oracle[22][19] + 12.96666666666667 gamma[22]
   + nu[22][16] >= 0
 dro_dual_ctr[170]: - oracle[22][0] + 37.76666666666667 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[171]: - oracle[22][1] + 8.666666666666671 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[172]: - oracle[22][2] + 62.86666666666667 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[173]: - oracle[22][3] + 12.96666666666666 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[174]: - oracle[22][4] + 6.266666666666669 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[175]: - oracle[22][5] + 70.2 gamma[22] + nu[22][17] >= 0
 dro_dual_ctr[176]: - oracle[22][6] + 167.5333333333333 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[177]: - oracle[22][7] + 6.933333333333336 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[178]: - oracle[22][8] + 5.43333333333333 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[179]: - oracle[22][9] + 10.76666666666667 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[1710]: - oracle[22][10] + 99.66666666666666 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[1711]: - oracle[22][11] + 10.26666666666666 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[1712]: - oracle[22][12] + 9.633333333333336 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[1713]: - oracle[22][13] + 11.06666666666667 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[1714]: - oracle[22][14] + 9.333333333333336 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[1715]: - oracle[22][15] + 31.26666666666667 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[1716]: - oracle[22][16] + 12.26666666666667 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[1717]: - oracle[22][17] + nu[22][17] >= 0
 dro_dual_ctr[1718]: - oracle[22][18] + 5.666666666666675 gamma[22]
   + nu[22][17] >= 0
 dro_dual_ctr[1719]: - oracle[22][19] + 5.5 gamma[22] + nu[22][17] >= 0
 dro_dual_ctr[180]: - oracle[22][0] + 37.3 gamma[22] + nu[22][18] >= 0
 dro_dual_ctr[181]: - oracle[22][1] + 13.73333333333334 gamma[22]
   + nu[22][18] >= 0
 dro_dual_ctr[182]: - oracle[22][2] + 62.4 gamma[22] + nu[22][18] >= 0
 dro_dual_ctr[183]: - oracle[22][3] + 18.63333333333334 gamma[22]
   + nu[22][18] >= 0
 dro_dual_ctr[184]: - oracle[22][4] + 11.93333333333334 gamma[22]
   + nu[22][18] >= 0
 dro_dual_ctr[185]: - oracle[22][5] + 69.73333333333335 gamma[22]
   + nu[22][18] >= 0
 dro_dual_ctr[186]: - oracle[22][6] + 167.0666666666667 gamma[22]
   + nu[22][18] >= 0
 dro_dual_ctr[187]: - oracle[22][7] + 12.2 gamma[22] + nu[22][18] >= 0
 dro_dual_ctr[188]: - oracle[22][8] + 3.166666666666675 gamma[22]
   + nu[22][18] >= 0
 dro_dual_ctr[189]: - oracle[22][9] + 16.43333333333334 gamma[22]
   + nu[22][18] >= 0
 dro_dual_ctr[1810]: - oracle[22][10] + 99.2 gamma[22] + nu[22][18] >= 0
 dro_dual_ctr[1811]: - oracle[22][11] + 15.93333333333334 gamma[22]
   + nu[22][18] >= 0
 dro_dual_ctr[1812]: - oracle[22][12] + 4.966666666666676 gamma[22]
   + nu[22][18] >= 0
 dro_dual_ctr[1813]: - oracle[22][13] + 10.6 gamma[22] + nu[22][18] >= 0
 dro_dual_ctr[1814]: - oracle[22][14] + 10.6 gamma[22] + nu[22][18] >= 0
 dro_dual_ctr[1815]: - oracle[22][15] + 30.8 gamma[22] + nu[22][18] >= 0
 dro_dual_ctr[1816]: - oracle[22][16] + 13.13333333333334 gamma[22]
   + nu[22][18] >= 0
 dro_dual_ctr[1817]: - oracle[22][17] + 5.666666666666675 gamma[22]
   + nu[22][18] >= 0
 dro_dual_ctr[1818]: - oracle[22][18] + nu[22][18] >= 0
 dro_dual_ctr[1819]: - oracle[22][19] + 11.16666666666667 gamma[22]
   + nu[22][18] >= 0
 dro_dual_ctr[190]: - oracle[22][0] + 37.06666666666668 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[191]: - oracle[22][1] + 9.366666666666665 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[192]: - oracle[22][2] + 62.16666666666667 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[193]: - oracle[22][3] + 7.466666666666668 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[194]: - oracle[22][4] + 1.23333333333334 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[195]: - oracle[22][5] + 69.5 gamma[22] + nu[22][19] >= 0
 dro_dual_ctr[196]: - oracle[22][6] + 166.8333333333333 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[197]: - oracle[22][7] + 7.63333333333333 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[198]: - oracle[22][8] + 10.93333333333332 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[199]: - oracle[22][9] + 8.466666666666661 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[1910]: - oracle[22][10] + 98.96666666666667 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[1911]: - oracle[22][11] + 4.766666666666666 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[1912]: - oracle[22][12] + 15.13333333333333 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[1913]: - oracle[22][13] + 10.56666666666666 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[1914]: - oracle[22][14] + 10.03333333333333 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[1915]: - oracle[22][15] + 31.56666666666668 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[1916]: - oracle[22][16] + 12.96666666666667 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[1917]: - oracle[22][17] + 5.5 gamma[22] + nu[22][19] >= 0
 dro_dual_ctr[1918]: - oracle[22][18] + 11.16666666666667 gamma[22]
   + nu[22][19] >= 0
 dro_dual_ctr[1919]: - oracle[22][19] + nu[22][19] >= 0
 cut[22,0,0]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][0] >= 126009.6071794694
 cut[22,0,1]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][1] >= 126009.6071794694
 cut[22,0,2]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][2] >= 126009.6071794694
 cut[22,0,3]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][3] >= 126009.6071794694
 cut[22,0,4]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][4] >= 126009.6071794694
 cut[22,0,5]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][5] >= 126009.6071794694
 cut[22,0,6]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][6] >= 126009.6071794694
 cut[22,0,7]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][7] >= 126009.6071794694
 cut[22,0,8]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][8] >= 126009.6071794694
 cut[22,0,9]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][9] >= 126009.6071794694
 cut[22,0,10]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][10] >= 126009.6071794694
 cut[22,0,11]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][11] >= 126009.6071794694
 cut[22,0,12]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][12] >= 126009.6071794694
 cut[22,0,13]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][13] >= 126009.6071794694
 cut[22,0,14]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][14] >= 126009.6071794694
 cut[22,0,15]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][15] >= 126009.6071794694
 cut[22,0,16]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][16] >= 126009.6071794694
 cut[22,0,17]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][17] >= 126009.6071794694
 cut[22,0,18]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][18] >= 126009.6071794694
 cut[22,0,19]: 1.8769990387538821e+07 P + 1.5 Q + 50.1064963535012 C
   + oracle[22][19] >= 126009.6071794694
Bounds
 M = -10000
 P0 = 0
 Q0 = 0
 W0 = 150
 C0 = 0
 M0 = -10000
 evapotranspiration = 40.1
 rainfall = 26.46666666666667
 oracle[22][0] >= -1e+08
 oracle[22][1] >= -1e+08
 oracle[22][2] >= -1e+08
 oracle[22][3] >= -1e+08
 oracle[22][4] >= -1e+08
 oracle[22][5] >= -1e+08
 oracle[22][6] >= -1e+08
 oracle[22][7] >= -1e+08
 oracle[22][8] >= -1e+08
 oracle[22][9] >= -1e+08
 oracle[22][10] >= -1e+08
 oracle[22][11] >= -1e+08
 oracle[22][12] >= -1e+08
 oracle[22][13] >= -1e+08
 oracle[22][14] >= -1e+08
 oracle[22][15] >= -1e+08
 oracle[22][16] >= -1e+08
 oracle[22][17] >= -1e+08
 oracle[22][18] >= -1e+08
 oracle[22][19] >= -1e+08
 nu[22][0] free
 nu[22][1] free
 nu[22][2] free
 nu[22][3] free
 nu[22][4] free
 nu[22][5] free
 nu[22][6] free
 nu[22][7] free
 nu[22][8] free
 nu[22][9] free
 nu[22][10] free
 nu[22][11] free
 nu[22][12] free
 nu[22][13] free
 nu[22][14] free
 nu[22][15] free
 nu[22][16] free
 nu[22][17] free
 nu[22][18] free
 nu[22][19] free
End
