\ Model Hydrovalley
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  1000 pour[0] + 1000 pour[1] + 1000 pour[2] + 1000 pour[3] + 1000 pour[4]
   + 1000 pour[5] + 1000 pour[6] + 1000 pour[7] + 1000 pour[8]
   + 1000 pour[9] - generation + lambda[0] + 0.1333333333333333 gamma[0]
   - 0.0333333333333333 nu[0][0] - 0.0333333333333333 nu[0][1]
   - 0.0333333333333333 nu[0][2] - 0.0333333333333333 nu[0][3]
   - 0.0333333333333333 nu[0][4] - 0.0333333333333333 nu[0][5]
   - 0.0333333333333333 nu[0][6] - 0.0333333333333333 nu[0][7]
   - 0.0333333333333333 nu[0][8] - 0.0333333333333333 nu[0][9]
   - 0.0333333333333333 nu[0][10] - 0.0333333333333333 nu[0][11]
   - 0.0333333333333333 nu[0][12] - 0.0333333333333333 nu[0][13]
   - 0.0333333333333333 nu[0][14] - 0.0333333333333333 nu[0][15]
   - 0.0333333333333333 nu[0][16] - 0.0333333333333333 nu[0][17]
   - 0.0333333333333333 nu[0][18] - 0.0333333333333333 nu[0][19]
   - 0.0333333333333333 nu[0][20] - 0.0333333333333333 nu[0][21]
   - 0.0333333333333333 nu[0][22] - 0.0333333333333333 nu[0][23]
   - 0.0333333333333333 nu[0][24] - 0.0333333333333333 nu[0][25]
   - 0.0333333333333333 nu[0][26] - 0.0333333333333333 nu[0][27]
   - 0.0333333333333333 nu[0][28] - 0.0333333333333333 nu[0][29]
Subject To
 AR_model_1[0]: inflow[0,1] - 0.982026172983832 inflow0[0,1]
   - innovations[0] = 0
 AR_model_1[1]: inflow[1,1] - 0.5893689920528696 inflow0[0,1]
   - 0.3000786041836117 inflow0[1,1] - innovations[1] = 0
 AR_model_1[2]: inflow[2,1] - 1.033778995074984 inflow0[1,1]
   - 1.220446599600729 inflow0[2,1] - innovations[2] = 0
 AR_model_1[3]: inflow[3,1] - 0.5750442087627947 inflow0[2,1]
   + 0.3886389399382055 inflow0[3,1] - innovations[3] = 0
 AR_model_1[4]: inflow[4,1] - 0.0483905741032211 inflow0[3,1]
   - 0.0243213958511511 inflow0[4,1] - innovations[4] = 0
 AR_model_1[5]: inflow[5,1] - 0.172021785580439 inflow0[4,1]
   - 0.3052992509691862 inflow0[5,1] - innovations[5] = 0
 AR_model_1[6]: inflow[6,1] - 0.4805188625734967 inflow0[5,1]
   - 0.8271367534814875 inflow0[6,1] - innovations[6] = 0
 AR_model_1[7]: inflow[7,1] - 0.3219316163727128 inflow0[6,1]
   - 0.1608375082464142 inflow0[7,1] - innovations[7] = 0
 AR_model_1[8]: inflow[8,1] - 0.847039536578803 inflow0[7,1]
   - 0.2668371636871334 inflow0[8,1] - innovations[8] = 0
 AR_model_1[9]: inflow[9,1] - 0.2565338508254507 inflow0[8,1]
   + 0.00257913188290043 inflow0[9,1] - innovations[9] = 0
 balance[0]: reservoir_level[0] - reservoir_level0[0]
   - 0.982026172983832 inflow0[0,1] - innovations[0] + outflow[0]
   + spill[0] - pour[0] = 0
 balance[1]: reservoir_level[1] - reservoir_level0[1]
   - 0.5893689920528696 inflow0[0,1] - 0.3000786041836117 inflow0[1,1]
   - innovations[1] - outflow[0] + outflow[1] - spill[0] + spill[1]
   - pour[1] = 0
 balance[2]: reservoir_level[2] - reservoir_level0[2]
   - 1.033778995074984 inflow0[1,1] - 1.220446599600729 inflow0[2,1]
   - innovations[2] - outflow[1] + outflow[2] - spill[1] + spill[2]
   - pour[2] = 0
 balance[3]: reservoir_level[3] - reservoir_level0[3]
   - 0.5750442087627947 inflow0[2,1] + 0.3886389399382055 inflow0[3,1]
   - innovations[3] - outflow[2] + outflow[3] - spill[2] + spill[3]
   - pour[3] = 0
 balance[4]: reservoir_level[4] - reservoir_level0[4]
   - 0.0483905741032211 inflow0[3,1] - 0.0243213958511511 inflow0[4,1]
   - innovations[4] - outflow[3] + outflow[4] - spill[3] + spill[4]
   - pour[4] = 0
 balance[5]: reservoir_level[5] - reservoir_level0[5]
   - 0.172021785580439 inflow0[4,1] - 0.3052992509691862 inflow0[5,1]
   - innovations[5] - outflow[4] + outflow[5] - spill[4] + spill[5]
   - pour[5] = 0
 balance[6]: reservoir_level[6] - reservoir_level0[6]
   - 0.4805188625734967 inflow0[5,1] - 0.8271367534814875 inflow0[6,1]
   - innovations[6] - outflow[5] + outflow[6] - spill[5] + spill[6]
   - pour[6] = 0
 balance[7]: reservoir_level[7] - reservoir_level0[7]
   - 0.3219316163727128 inflow0[6,1] - 0.1608375082464142 inflow0[7,1]
   - innovations[7] - outflow[6] + outflow[7] - spill[6] + spill[7]
   - pour[7] = 0
 balance[8]: reservoir_level[8] - reservoir_level0[8]
   - 0.847039536578803 inflow0[7,1] - 0.2668371636871334 inflow0[8,1]
   - innovations[8] - outflow[7] + outflow[8] - spill[7] + spill[8]
   - pour[8] = 0
 balance[9]: reservoir_level[9] - reservoir_level0[9]
   - 0.2565338508254507 inflow0[8,1] + 0.00257913188290043 inflow0[9,1]
   - innovations[9] - outflow[8] + outflow[9] - spill[8] + spill[9]
   - pour[9] = 0
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
 dro_dual_ctr[0]: - oracle[0][0] + lambda[0] - nu[0][0] >= 0
 dro_dual_ctr[1]: - oracle[0][1] + lambda[0] - nu[0][1] >= 0
 dro_dual_ctr[2]: - oracle[0][2] + lambda[0] - nu[0][2] >= 0
 dro_dual_ctr[3]: - oracle[0][3] + lambda[0] - nu[0][3] >= 0
 dro_dual_ctr[4]: - oracle[0][4] + lambda[0] - nu[0][4] >= 0
 dro_dual_ctr[5]: - oracle[0][5] + lambda[0] - nu[0][5] >= 0
 dro_dual_ctr[6]: - oracle[0][6] + lambda[0] - nu[0][6] >= 0
 dro_dual_ctr[7]: - oracle[0][7] + lambda[0] - nu[0][7] >= 0
 dro_dual_ctr[8]: - oracle[0][8] + lambda[0] - nu[0][8] >= 0
 dro_dual_ctr[9]: - oracle[0][9] + lambda[0] - nu[0][9] >= 0
 dro_dual_ctr[10]: - oracle[0][10] + lambda[0] - nu[0][10] >= 0
 dro_dual_ctr[11]: - oracle[0][11] + lambda[0] - nu[0][11] >= 0
 dro_dual_ctr[12]: - oracle[0][12] + lambda[0] - nu[0][12] >= 0
 dro_dual_ctr[13]: - oracle[0][13] + lambda[0] - nu[0][13] >= 0
 dro_dual_ctr[14]: - oracle[0][14] + lambda[0] - nu[0][14] >= 0
 dro_dual_ctr[15]: - oracle[0][15] + lambda[0] - nu[0][15] >= 0
 dro_dual_ctr[16]: - oracle[0][16] + lambda[0] - nu[0][16] >= 0
 dro_dual_ctr[17]: - oracle[0][17] + lambda[0] - nu[0][17] >= 0
 dro_dual_ctr[18]: - oracle[0][18] + lambda[0] - nu[0][18] >= 0
 dro_dual_ctr[19]: - oracle[0][19] + lambda[0] - nu[0][19] >= 0
 dro_dual_ctr[20]: - oracle[0][20] + lambda[0] - nu[0][20] >= 0
 dro_dual_ctr[21]: - oracle[0][21] + lambda[0] - nu[0][21] >= 0
 dro_dual_ctr[22]: - oracle[0][22] + lambda[0] - nu[0][22] >= 0
 dro_dual_ctr[23]: - oracle[0][23] + lambda[0] - nu[0][23] >= 0
 dro_dual_ctr[24]: - oracle[0][24] + lambda[0] - nu[0][24] >= 0
 dro_dual_ctr[25]: - oracle[0][25] + lambda[0] - nu[0][25] >= 0
 dro_dual_ctr[26]: - oracle[0][26] + lambda[0] - nu[0][26] >= 0
 dro_dual_ctr[27]: - oracle[0][27] + lambda[0] - nu[0][27] >= 0
 dro_dual_ctr[28]: - oracle[0][28] + lambda[0] - nu[0][28] >= 0
 dro_dual_ctr[29]: - oracle[0][29] + lambda[0] - nu[0][29] >= 0
 cut[0,0,0]: 1000 reservoir_level[0] + 998.1080000191745 reservoir_level[1]
   + 909.1345314236662 reservoir_level[2]
   + 372.854612206289 reservoir_level[3]
   + 370.6565523250109 reservoir_level[4]
   + 350.8053965594901 reservoir_level[5]
   + 6.496913426963605 reservoir_level[6]
   + 4.298853539917175 reservoir_level[7]
   + 2.406853560186514 reservoir_level[8]
   + 1.546853567549811 reservoir_level[9] - 421.2864262475495 inflow[0,1]
   + 177.7316761590393 inflow[1,1] + 852.1289976313595 inflow[2,1]
   - 285.7794470842946 inflow[3,1] - 102.7873443340154 inflow[4,1]
   + 8.433523926385535 inflow[5,1] - 6.98413816893737 inflow[6,1]
   + 0.2769869180379113 inflow[7,1] + 1.960936902767379 inflow[8,1]
   + 0.6799193338195268 inflow[9,1] + oracle[0][0] >= 16596.69286642363
 cut[0,0,1]: 999.9999999998513 reservoir_level[0]
   + 998.1080000060742 reservoir_level[1]
   + 898.9855489299578 reservoir_level[2]
   + 299.746801738152 reservoir_level[3]
   + 297.5487566019856 reservoir_level[4]
   + 277.8412247850889 reservoir_level[5]
   + 6.277273529651282 reservoir_level[6]
   + 4.079228387703813 reservoir_level[7]
   + 2.187228395288049 reservoir_level[8]
   + 1.327228398090486 reservoir_level[9] - 419.2391894487191 inflow[0,1]
   + 158.7250803923665 inflow[1,1] + 795.7701477744199 inflow[2,1]
   - 224.5423675849742 inflow[3,1] - 89.2247579322252 inflow[4,1]
   + 7.621021149446354 inflow[5,1] - 6.71529583337955 inflow[6,1]
   + 0.2369618300438851 inflow[7,1] + 1.768906983533664 inflow[8,1]
   + 0.5904697460873212 inflow[9,1] + oracle[0][1] >= 963.9014965896131
 cut[0,0,2]: 1000.000000002093 reservoir_level[0]
   + 998.1079999189909 reservoir_level[1]
   + 893.9221553635372 reservoir_level[2]
   + 266.3363677734797 reservoir_level[3]
   + 264.1371991155768 reservoir_level[4]
   + 244.8273001004796 reservoir_level[5]
   + 6.168834661468886 reservoir_level[6]
   + 3.969665969870608 reservoir_level[7]
   + 2.077665857070285 reservoir_level[8]
   + 1.217665810628648 reservoir_level[9] - 418.0281370301312 inflow[0,1]
   + 149.5384892456317 inflow[1,1] + 769.3499748162027 inflow[2,1]
   - 197.1111523274355 inflow[3,1] - 82.26047844551007 inflow[4,1]
   + 7.241732898250456 inflow[5,1] - 6.580028380207473 inflow[6,1]
   + 0.2174026976120812 inflow[7,1] + 1.673584407818326 inflow[8,1]
   + 0.5456603999488894 inflow[9,1] + oracle[0][2] >= -2363.307261802427
 cut[0,0,3]: 1000 reservoir_level[0] + 998.1079999896408 reservoir_level[1]
   + 888.6684663485056 reservoir_level[2]
   + 244.0985459556908 reservoir_level[3]
   + 241.9016646110968 reservoir_level[4]
   + 223.9556846201467 reservoir_level[5]
   + 6.082509898840147 reservoir_level[6]
   + 3.885628544748227 reservoir_level[7]
   + 1.993628531809218 reservoir_level[8]
   + 0.860000013837158 reservoir_level[9] - 414.9145570699746 inflow[0,1]
   + 141.9996496589193 inflow[1,1] + 749.5731217142828 inflow[2,1]
   - 180.4193551347392 inflow[3,1] - 75.56875936887195 inflow[4,1]
   + 6.975225135574983 inflow[5,1] - 6.464258026126013 inflow[6,1]
   + 0.2021926180651664 inflow[7,1] + 1.553258143605711 inflow[8,1]
   + 0.4352527866691449 inflow[9,1] + oracle[0][3] >= -6040.729096260773
 cut[0,0,4]: 1000.000000006668 reservoir_level[0]
   + 998.1080000609995 reservoir_level[1]
   + 888.5115415125492 reservoir_level[2]
   + 240.2256353091028 reservoir_level[3]
   + 238.0271007145262 reservoir_level[4]
   + 220.3473041554921 reservoir_level[5]
   + 6.07215890736099 reservoir_level[6]
   + 3.873624296065807 reservoir_level[7]
   + 1.981624326731242 reservoir_level[8]
   + 0.8599996762002277 reservoir_level[9] - 414.7436340006553 inflow[0,1]
   + 140.8982310927561 inflow[1,1] + 746.5863317235454 inflow[2,1]
   - 177.5543469655448 inflow[3,1] - 74.34717453254412 inflow[4,1]
   + 6.932433672682616 inflow[5,1] - 6.448518950036807 inflow[6,1]
   + 0.2003486791427136 inflow[7,1] + 1.545296609065743 inflow[8,1]
   + 0.4339230957495531 inflow[9,1] + oracle[0][4] >= -9218.502848908702
 cut[0,0,5]: 911.4981203307572 reservoir_level[0]
   + 909.6061203554355 reservoir_level[1]
   + 886.0970060987695 reservoir_level[2]
   + 229.2813490891015 reservoir_level[3]
   + 227.0832660273353 reservoir_level[4]
   + 210.1499836603183 reservoir_level[5]
   + 6.030569193499087 reservoir_level[6]
   + 3.83248612795478 reservoir_level[7]
   + 1.940486155488152 reservoir_level[8]
   + 0.8599999600500097 reservoir_level[9] - 419.2276123124012 inflow[0,1]
   + 154.8175572693843 inflow[1,1] + 736.8861873437301 inflow[2,1]
   - 169.4570732488325 inflow[3,1] - 70.89956760447069 inflow[4,1]
   + 6.801826086877041 inflow[5,1] - 6.39172632685248 inflow[6,1]
   + 0.1930176573424356 inflow[7,1] + 1.516550368974762 inflow[8,1]
   + 0.4299663011305619 inflow[9,1] + oracle[0][5] >= -12541.42374617035
 cut[0,0,6]: 905.5631696814271 reservoir_level[0]
   + 903.6711695977111 reservoir_level[1]
   + 883.1354303784397 reservoir_level[2]
   + 214.2050638625981 reservoir_level[3]
   + 212.0065608197099 reservoir_level[4]
   + 196.0827595441055 reservoir_level[5]
   + 5.975322690562367 reservoir_level[6]
   + 3.776819501725081 reservoir_level[7]
   + 1.884819304518555 reservoir_level[8]
   + 0.8599995719917509 reservoir_level[9] - 417.6519194258393 inflow[0,1]
   + 150.9866711048306 inflow[1,1] + 723.7555444485838 inflow[2,1]
   - 158.2999384661385 inflow[3,1] - 66.14569064541109 inflow[4,1]
   + 6.623503430073717 inflow[5,1] - 6.315704827657184 inflow[6,1]
   + 0.1832332144406888 inflow[7,1] + 1.477917234205348 inflow[8,1]
   + 0.4246121899471654 inflow[9,1] + oracle[0][6] >= -15021.22997432682
 cut[0,0,7]: 904.4358359447234 reservoir_level[0]
   + 902.5438359405532 reservoir_level[1]
   + 882.4808414338156 reservoir_level[2]
   + 212.1499665251102 reservoir_level[3]
   + 209.9520463217442 reservoir_level[4]
   + 194.1744497663636 reservoir_level[5]
   + 5.966448923213449 reservoir_level[6]
   + 3.768528711111378 reservoir_level[7]
   + 1.876528704646636 reservoir_level[8]
   + 0.8600000033252915 reservoir_level[9] - 417.3196403880738 inflow[0,1]
   + 150.4504007231154 inflow[1,1] + 721.8052512889234 inflow[2,1]
   - 156.7797812594141 inflow[3,1] - 65.50521040148733 inflow[4,1]
   + 6.598228560646563 inflow[5,1] - 6.304017334208341 inflow[6,1]
   + 0.1816833394042239 inflow[7,1] + 1.472019830910646 inflow[8,1]
   + 0.4238590338221525 inflow[9,1] + oracle[0][7] >= -15868.75091847764
 cut[0,0,8]: 900.4170845888375 reservoir_level[0]
   + 898.5250845988106 reservoir_level[1]
   + 880.4640213159283 reservoir_level[2]
   + 202.9967893423775 reservoir_level[3]
   + 200.7990629777904 reservoir_level[4]
   + 185.598534462378 reservoir_level[5]
   + 5.931004297968162 reservoir_level[6]
   + 3.733277931547933 reservoir_level[7]
   + 1.841277945945966 reservoir_level[8]
   + 0.8599999929768255 reservoir_level[9] - 416.2619097627223 inflow[0,1]
   + 148.1397936739725 inflow[1,1] + 713.6928439388905 inflow[2,1]
   - 150.0019554103778 inflow[3,1] - 62.5998773323725 inflow[4,1]
   + 6.488224907662833 inflow[5,1] - 6.255997438959968 inflow[6,1]
   + 0.175359872805342 inflow[7,1] + 1.447433218361795 inflow[8,1]
   + 0.4206161012438318 inflow[9,1] + oracle[0][8] >= -17766.91693933532
 cut[0,0,9]: 898.2014994349262 reservoir_level[0]
   + 896.3095004735927 reservoir_level[1]
   + 879.1954674875099 reservoir_level[2]
   + 199.1860550409257 reservoir_level[3]
   + 196.9895250145583 reservoir_level[4]
   + 182.0259520424215 reservoir_level[5]
   + 5.913535239480789 reservoir_level[6]
   + 3.717005545795102 reservoir_level[7]
   + 1.825006866714313 reservoir_level[8]
   + 0.8599991818244742 reservoir_level[9] - 415.6206703715328 inflow[0,1]
   + 147.1534269418187 inflow[1,1] + 710.0411728997302 inflow[2,1]
   - 147.1796181341435 inflow[3,1] - 61.39058025440937 inflow[4,1]
   + 6.440297953722459 inflow[5,1] - 6.233317299026833 inflow[6,1]
   + 0.1722596138779656 inflow[7,1] + 1.435823083958929 inflow[8,1]
   + 0.4192279282759816 inflow[9,1] + oracle[0][9] >= -18750.98257129714
 cut[0,0,10]: 897.4458556564258 reservoir_level[0]
   + 895.5538556668889 reservoir_level[1]
   + 879.0701750406635 reservoir_level[2]
   + 196.733075341263 reservoir_level[3]
   + 194.5354955681066 reservoir_level[4]
   + 179.7230536612112 reservoir_level[5]
   + 5.906571250977203 reservoir_level[6]
   + 3.708991476831963 reservoir_level[7]
   + 1.816991493278194 reservoir_level[8]
   + 0.8599999921967825 reservoir_level[9] - 415.5335388411921 inflow[0,1]
   + 146.5984398964565 inflow[1,1] + 708.132102936544 inflow[2,1]
   - 145.3626725684145 inflow[3,1] - 60.61100590862888 inflow[4,1]
   + 6.412819362382996 inflow[5,1] - 6.223035588404041 inflow[6,1]
   + 0.1709894422368608 inflow[7,1] + 1.43049983793062 inflow[8,1]
   + 0.4184199364969129 inflow[9,1] + oracle[0][10] >= -19404.2005051981
 cut[0,0,11]: 896.4733080652378 reservoir_level[0]
   + 894.5813080996046 reservoir_level[1]
   + 878.5714086128406 reservoir_level[2]
   + 194.8324301347333 reservoir_level[3]
   + 192.6350772125918 reservoir_level[4]
   + 177.9370549239943 reservoir_level[5]
   + 5.898629952289822 reservoir_level[6]
   + 3.701277035292499 reservoir_level[7]
   + 1.809277083948579 reservoir_level[8]
   + 0.8599999744456434 reservoir_level[9] - 415.276802651999 inflow[0,1]
   + 146.1214518573245 inflow[1,1] + 706.3962889233921 inflow[2,1]
   - 143.9545821723376 inflow[3,1] - 60.00557655697509 inflow[4,1]
   + 6.389492893737586 inflow[5,1] - 6.212498382443621 inflow[6,1]
   + 0.1695667174691056 inflow[7,1] + 1.425073771457136 inflow[8,1]
   + 0.4177457036191757 inflow[9,1] + oracle[0][11] >= -19796.46302940493
 cut[0,0,12]: 895.6582332607665 reservoir_level[0]
   + 880.564537063174 reservoir_level[1]
   + 878.3753774148346 reservoir_level[2]
   + 192.5318200385617 reservoir_level[3]
   + 190.333661147942 reservoir_level[4]
   + 175.7738017314566 reservoir_level[5]
   + 5.891604365168286 reservoir_level[6]
   + 3.693445454127281 reservoir_level[7]
   + 1.801445391435318 reservoir_level[8]
   + 0.8600000261892934 reservoir_level[9] - 409.1126929249433 inflow[0,1]
   + 148.2293319824535 inflow[1,1] + 704.5555583010785 inflow[2,1]
   - 142.250080809216 inflow[3,1] - 59.27267845204517 inflow[4,1]
   + 6.363308509971437 inflow[5,1] - 6.202363904050804 inflow[6,1]
   + 0.1682887326933754 inflow[7,1] + 1.419822801237003 inflow[8,1]
   + 0.4169829306602449 inflow[9,1] + oracle[0][12] >= -20959.61405551908
 cut[0,0,13]: 894.4523689122312 reservoir_level[0]
   + 879.8866564251684 reservoir_level[1]
   + 877.6984557276453 reservoir_level[2]
   + 190.6702757950231 reservoir_level[3]
   + 188.4728334248642 reservoir_level[4]
   + 174.0271589552362 reservoir_level[5]
   + 5.88270040938079 reservoir_level[6]
   + 3.685258069028705 reservoir_level[7]
   + 1.793258204029254 reservoir_level[8]
   + 0.8599999213517434 reservoir_level[9] - 409.0149140482108 inflow[0,1]
   + 147.649192703425 inflow[1,1] + 702.7349593525964 inflow[2,1]
   - 140.871053133905 inflow[3,1] - 58.68250405997964 inflow[4,1]
   + 6.339613013560993 inflow[5,1] - 6.190936190939714 inflow[6,1]
   + 0.16670561390475 inflow[7,1] + 1.413954846665361 inflow[8,1]
   + 0.416306704578732 inflow[9,1] + oracle[0][13] >= -21522.85753352258
 cut[0,0,14]: 893.9387155070548 reservoir_level[0]
   + 879.6803839310555 reservoir_level[1]
   + 877.4922486632354 reservoir_level[2]
   + 189.5306005228594 reservoir_level[3]
   + 187.3330384418826 reservoir_level[4]
   + 172.9560205101755 reservoir_level[5]
   + 5.87853080746427 reservoir_level[6]
   + 3.680968727541055 reservoir_level[7]
   + 1.788968753347358 reservoir_level[8]
   + 0.8599999864542024 reservoir_level[9] - 409.0456760999189 inflow[0,1]
   + 147.3167476876592 inflow[1,1] + 701.7533630638735 inflow[2,1]
   - 140.0266613834408 inflow[3,1] - 58.32007327931619 inflow[4,1]
   + 6.326100657601775 inflow[5,1] - 6.185222083298084 inflow[6,1]
   + 0.165952030850723 inflow[7,1] + 1.410998056263259 inflow[8,1]
   + 0.4159167020912093 inflow[9,1] + oracle[0][14] >= -21873.13574093341
 cut[0,0,15]: 893.1129755133484 reservoir_level[0]
   + 879.2560274424241 reservoir_level[1]
   + 877.0683923046723 reservoir_level[2]
   + 188.0894743910296 reservoir_level[3]
   + 185.8921887644285 reservoir_level[4]
   + 171.5986301228612 reservoir_level[5]
   + 5.872188351956353 reservoir_level[6]
   + 3.674903005879348 reservoir_level[7]
   + 1.782903887867207 reservoir_level[8]
   + 0.8599993770175854 reservoir_level[9] - 409.0137218734392 inflow[0,1]
   + 146.8785663352575 inflow[1,1] + 700.4079163050615 inflow[2,1]
   - 138.9585850904886 inflow[3,1] - 57.85974112957499 inflow[4,1]
   + 6.308146829001021 inflow[5,1] - 6.176923513381052 inflow[6,1]
   + 0.1648139096486956 inflow[7,1] + 1.406710221204586 inflow[8,1]
   + 0.4154058362951389 inflow[9,1] + oracle[0][15] >= -22353.71351483881
 cut[0,0,16]: 892.6904344661644 reservoir_level[0]
   + 879.3942061985076 reservoir_level[1]
   + 877.2054911642806 reservoir_level[2]
   + 186.4240076722735 reservoir_level[3]
   + 184.2254094336016 reservoir_level[4]
   + 170.0316913262695 reservoir_level[5]
   + 5.868825436017187 reservoir_level[6]
   + 3.670227195425877 reservoir_level[7]
   + 1.778227211176301 reservoir_level[8]
   + 0.8599999924918237 reservoir_level[9] - 409.31714405946 inflow[0,1]
   + 146.4405942123304 inflow[1,1] + 699.253340958586 inflow[2,1]
   - 137.7244222855933 inflow[3,1] - 57.33088030890633 inflow[4,1]
   + 6.290587778227755 inflow[5,1] - 6.171367383941494 inflow[6,1]
   + 0.1641789010540943 inflow[7,1] + 1.403781450290283 inflow[8,1]
   + 0.4148992712022339 inflow[9,1] + oracle[0][16] >= -23136.38163594768
 cut[0,0,17]: 889.8403876896405 reservoir_level[0]
   + 877.9957154431909 reservoir_level[1]
   + 875.8086444934294 reservoir_level[2]
   + 181.5805001327097 reservoir_level[3]
   + 179.3827092576424 reservoir_level[4]
   + 165.4354328519381 reservoir_level[5]
   + 5.847057767427966 reservoir_level[6]
   + 3.649268423770271 reservoir_level[7]
   + 1.757273535190744 reservoir_level[8]
   + 0.8599967835235829 reservoir_level[9] - 409.2707308730825 inflow[0,1]
   + 144.9709301320271 inflow[1,1] + 694.7487154142437 inflow[2,1]
   - 134.1311018386663 inflow[3,1] - 55.7627968913764 inflow[4,1]
   + 6.229640645801136 inflow[5,1] - 6.143046012488137 inflow[6,1]
   + 0.1602252074897154 inflow[7,1] + 1.388986306966795 inflow[8,1]
   + 0.4132049146793749 inflow[9,1] + oracle[0][17] >= -24774.90939912668
 cut[0,0,18]: 888.3761395581851 reservoir_level[0]
   + 877.2775511163336 reservoir_level[1]
   + 875.0913932616685 reservoir_level[2]
   + 179.2575479374681 reservoir_level[3]
   + 177.0602473040334 reservoir_level[4]
   + 163.2275422671042 reservoir_level[5]
   + 5.83630431629895 reservoir_level[6]
   + 3.639003491149486 reservoir_level[7]
   + 1.747002693791474 reservoir_level[8]
   + 0.6922096573393089 reservoir_level[9] - 409.2495225493718 inflow[0,1]
   + 144.2609128333544 inflow[1,1] + 692.5579383822339 inflow[2,1]
   - 132.4072878395791 inflow[3,1] - 55.00924974473105 inflow[4,1]
   + 6.20014716197513 inflow[5,1] - 6.12916227103412 inflow[6,1]
   + 0.1582670514690616 inflow[7,1] + 1.35301571246517 inflow[8,1]
   + 0.3645051268988924 inflow[9,1] + oracle[0][18] >= -26235.59170674399
 cut[0,0,19]: 887.8300557780186 reservoir_level[0]
   + 877.124221904661 reservoir_level[1]
   + 874.9380426987021 reservoir_level[2]
   + 178.272014049864 reservoir_level[3]
   + 176.0745296697719 reservoir_level[4]
   + 162.2913920322642 reservoir_level[5]
   + 5.832685701104213 reservoir_level[6]
   + 3.635201322293395 reservoir_level[7]
   + 1.743201367584933 reservoir_level[8]
   + 0.6908910008252975 reservoir_level[9] - 409.3471799413821 inflow[0,1]
   + 143.9761538973308 inflow[1,1] + 691.7247294273882 inflow[2,1]
   - 131.6758877251835 inflow[3,1] - 54.69153350918015 inflow[4,1]
   + 6.188413818091634 inflow[5,1] - 6.124254534168625 inflow[6,1]
   + 0.1575977421111638 inflow[7,1] + 1.350190843954406 inflow[8,1]
   + 0.3638108281110565 inflow[9,1] + oracle[0][19] >= -26729.13500457643
 cut[0,0,20]: 887.1292271974377 reservoir_level[0]
   + 876.7580059867769 reservoir_level[1]
   + 874.5724171479494 reservoir_level[2]
   + 177.4181232184957 reservoir_level[3]
   + 175.2210716041945 reservoir_level[4]
   + 161.4845110849688 reservoir_level[5]
   + 5.828229814742395 reservoir_level[6]
   + 3.631178204092566 reservoir_level[7]
   + 1.739178283077837 reservoir_level[8]
   + 0.6896326663689631 reservoir_level[9] - 409.3195795159932 inflow[0,1]
   + 143.7040286176759 inflow[1,1] + 690.8541392240496 inflow[2,1]
   - 131.0425420249878 inflow[3,1] - 54.41899291879059 inflow[4,1]
   + 6.17723348998813 inflow[5,1] - 6.118670612780649 inflow[6,1]
   + 0.1568004930762105 inflow[7,1] + 1.347078861334286 inflow[8,1]
   + 0.3631482194940166 inflow[9,1] + oracle[0][20] >= -28262.00714168244
 cut[0,0,21]: 886.6383590966868 reservoir_level[0]
   + 876.5498252128792 reservoir_level[1]
   + 874.3644602511325 reservoir_level[2]
   + 176.6243458109015 reservoir_level[3]
   + 174.4273388456981 reservoir_level[4]
   + 160.7304504032655 reservoir_level[5]
   + 5.82485557840784 reservoir_level[6]
   + 3.627848595910884 reservoir_level[7]
   + 1.735848561519167 reservoir_level[8]
   + 0.6885171925983851 reservoir_level[9] - 409.3427895213655 inflow[0,1]
   + 143.4654622898375 inflow[1,1] + 690.1291709615626 inflow[2,1]
   - 130.4534784854687 inflow[3,1] - 54.1624732853891 inflow[4,1]
   + 6.167403150524708 inflow[5,1] - 6.114236226520875 inflow[6,1]
   + 0.1561856504211463 inflow[7,1] + 1.344563095179384 inflow[8,1]
   + 0.3625608005834502 inflow[9,1] + oracle[0][21] >= -29853.01281239374
 cut[0,0,22]: 886.2928736347553 reservoir_level[0]
   + 876.356675373865 reservoir_level[1]
   + 874.1716416436881 reservoir_level[2]
   + 176.1931516031198 reservoir_level[3]
   + 173.9963804730681 reservoir_level[4]
   + 160.3177530094017 reservoir_level[5]
   + 5.822445515120531 reservoir_level[6]
   + 3.625674391544469 reservoir_level[7]
   + 1.733674448680659 reservoir_level[8]
   + 0.6878745972384207 reservoir_level[9] - 409.3172353647022 inflow[0,1]
   + 143.3271537487261 inflow[1,1] + 689.6843300753019 inflow[2,1]
   - 130.1331601688302 inflow[3,1] - 54.02101113363244 inflow[4,1]
   + 6.161581310691405 inflow[5,1] - 6.11124580602668 inflow[6,1]
   + 0.1557474576196277 inflow[7,1] + 1.342881302943126 inflow[8,1]
   + 0.3622224168341493 inflow[9,1] + oracle[0][22] >= -30425.93941104715
 cut[0,0,23]: 885.3673268087452 reservoir_level[0]
   + 875.9785977464779 reservoir_level[1]
   + 873.7939670270845 reservoir_level[2]
   + 174.7382259011686 reservoir_level[3]
   + 172.541532611213 reservoir_level[4]
   + 158.9325523068855 reservoir_level[5]
   + 5.816183981235106 reservoir_level[6]
   + 3.619490462006262 reservoir_level[7]
   + 1.727489602685498 reservoir_level[8]
   + 0.6858484200502121 reservoir_level[9] - 409.3751206551716 inflow[0,1]
   + 142.8902635876752 inflow[1,1] + 688.3576570289729 inflow[2,1]
   - 129.0530795602203 inflow[3,1] - 53.54947936215558 inflow[4,1]
   + 6.143488580678965 inflow[5,1] - 6.103052658184721 inflow[6,1]
   + 0.1545986511302559 inflow[7,1] + 1.338212963272511 inflow[8,1]
   + 0.3611553676140462 inflow[9,1] + oracle[0][23] >= -32462.98993235871
 cut[0,0,24]: 885.4694139869337 reservoir_level[0]
   + 876.3482897658653 reservoir_level[1]
   + 874.162444345328 reservoir_level[2]
   + 174.0865452707329 reservoir_level[3]
   + 171.8884769941054 reservoir_level[4]
   + 158.307993837661 reservoir_level[5]
   + 5.816675304998657 reservoir_level[6]
   + 3.618606458136242 reservoir_level[7]
   + 1.72660455970061 reservoir_level[8]
   + 0.6852139712707181 reservoir_level[9] - 409.6645993613877 inflow[0,1]
   + 142.7532543089343 inflow[1,1] + 688.1068398874797 inflow[2,1]
   - 128.5690079418532 inflow[3,1] - 53.33604053419528 inflow[4,1]
   + 6.137966507100599 inflow[5,1] - 6.102705708346143 inflow[6,1]
   + 0.1546419339699507 inflow[7,1] + 1.337820643701307 inflow[8,1]
   + 0.3608211858964829 inflow[9,1] + oracle[0][24] >= -33304.84642109794
 cut[0,0,25]: 884.9757728467187 reservoir_level[0]
   + 876.2915837355063 reservoir_level[1]
   + 874.105509082411 reservoir_level[2]
   + 173.1634283612333 reservoir_level[3]
   + 170.9649085377189 reservoir_level[4]
   + 157.4313753743556 reservoir_level[5]
   + 5.813959547276908 reservoir_level[6]
   + 3.615439723405365 reservoir_level[7]
   + 1.544059375573175 reservoir_level[8]
   + 0.6840593775545762 reservoir_level[9] - 409.8291775277373 inflow[0,1]
   + 142.4959951052211 inflow[1,1] + 687.3815966519246 inflow[2,1]
   - 127.8838161226767 inflow[3,1] - 53.04028433946974 inflow[4,1]
   + 6.127530321365803 inflow[5,1] - 6.098807667789835 inflow[6,1]
   + 0.0748030784516465 inflow[7,1] + 1.16219692702065 inflow[8,1]
   + 0.3602133537761327 inflow[9,1] + oracle[0][25] >= -36139.34732252005
 cut[0,0,26]: 884.7368270755643 reservoir_level[0]
   + 876.1879781955784 reservoir_level[1]
   + 874.002051844237 reservoir_level[2]
   + 172.8478392647426 reservoir_level[3]
   + 170.6493891824861 reservoir_level[4]
   + 157.1270140399415 reservoir_level[5]
   + 5.812379234389366 reservoir_level[6]
   + 3.613929159797488 reservoir_level[7]
   + 1.543621686964517 reservoir_level[8]
   + 0.6836217075117844 reservoir_level[9] - 409.8393494320845 inflow[0,1]
   + 142.3988791017364 inflow[1,1] + 687.0800001704783 inflow[2,1]
   - 127.6490964372927 inflow[3,1] - 52.93581010053323 inflow[4,1]
   + 6.123404183027969 inflow[5,1] - 6.096812451198028 inflow[6,1]
   + 0.0749828516676477 inflow[7,1] + 1.162087756005362 inflow[8,1]
   + 0.3599828833146708 inflow[9,1] + oracle[0][26] >= -36451.59386713216
 cut[0,0,27]: 884.4043237755973 reservoir_level[0]
   + 876.1359294443812 reservoir_level[1]
   + 873.949904017789 reservoir_level[2]
   + 172.2882397445965 reservoir_level[3]
   + 170.0895497879395 reservoir_level[4]
   + 156.5956118001958 reservoir_level[5]
   + 5.810684411689837 reservoir_level[6]
   + 3.611994451283982 reservoir_level[7]
   + 1.542919762336497 reservoir_level[8]
   + 0.6829197632476082 reservoir_level[9] - 409.9380977467185 inflow[0,1]
   + 142.241001562145 inflow[1,1] + 686.6290867847804 inflow[2,1]
   - 127.2337240561035 inflow[3,1] - 52.75653862755051 inflow[4,1]
   + 6.117042054231953 inflow[5,1] - 6.094404187826647 inflow[6,1]
   + 0.0752160756360067 inflow[7,1] + 1.161884864067835 inflow[8,1]
   + 0.3596132675887727 inflow[9,1] + oracle[0][27] >= -39606.69206449418
 cut[0,0,28]: 882.9025464254453 reservoir_level[0]
   + 875.4346842649264 reservoir_level[1]
   + 873.2498668261994 reservoir_level[2]
   + 170.7791518376694 reservoir_level[3]
   + 168.5813152106286 reservoir_level[4]
   + 155.135874607523 reservoir_level[5]
   + 5.801831901284802 reservoir_level[6]
   + 3.603995582937506 reservoir_level[7]
   + 1.540799185338206 reservoir_level[8]
   + 0.6807997708818689 reservoir_level[9] - 409.9617136974976 inflow[0,1]
   + 141.7541964824678 inflow[1,1] + 685.0554184867002 inflow[2,1]
   - 126.1105992367183 inflow[3,1] - 52.25644813346385 inflow[4,1]
   + 6.096299989110294 inflow[5,1] - 6.083616978474161 inflow[6,1]
   + 0.0761818462599935 inflow[7,1] + 1.161428368638135 inflow[8,1]
   + 0.3584970547198658 inflow[9,1] + oracle[0][28] >= -43418.56498127655
 cut[0,0,29]: 882.1839871651838 reservoir_level[0]
   + 875.4412472669766 reservoir_level[1]
   + 873.2559330937494 reservoir_level[2]
   + 169.6036573703993 reservoir_level[3]
   + 167.4050412024505 reservoir_level[4]
   + 153.9887640550691 reservoir_level[5]
   + 5.798173576423554 reservoir_level[6]
   + 3.599557406138442 reservoir_level[7]
   + 1.539469001925525 reservoir_level[8]
   + 0.6794690050499499 reservoir_level[9] - 410.2870940143481 inflow[0,1]
   + 141.4351779954649 inflow[1,1] + 684.1815872511224 inflow[2,1]
   - 125.2347160187828 inflow[3,1] - 51.86219965882951 inflow[4,1]
   + 6.082657122185346 inflow[5,1] - 6.078507962019065 inflow[6,1]
   + 0.0768189652153755 inflow[7,1] + 1.161299539442233 inflow[8,1]
   + 0.3577963211271326 inflow[9,1] + oracle[0][29] >= -45455.79441053069
 norm_dro_ctr: [ - gamma[0] ^2 + nu[0][0] ^2 + nu[0][1] ^2 + nu[0][2] ^2
   + nu[0][3] ^2 + nu[0][4] ^2 + nu[0][5] ^2 + nu[0][6] ^2 + nu[0][7] ^2
   + nu[0][8] ^2 + nu[0][9] ^2 + nu[0][10] ^2 + nu[0][11] ^2 + nu[0][12] ^2
   + nu[0][13] ^2 + nu[0][14] ^2 + nu[0][15] ^2 + nu[0][16] ^2
   + nu[0][17] ^2 + nu[0][18] ^2 + nu[0][19] ^2 + nu[0][20] ^2
   + nu[0][21] ^2 + nu[0][22] ^2 + nu[0][23] ^2 + nu[0][24] ^2
   + nu[0][25] ^2 + nu[0][26] ^2 + nu[0][27] ^2 + nu[0][28] ^2
   + nu[0][29] ^2 ] <= 0
Bounds
 reservoir_level[0] <= 200
 reservoir_level[1] <= 200
 reservoir_level[2] <= 200
 reservoir_level[3] <= 200
 reservoir_level[4] <= 200
 reservoir_level[5] <= 200
 reservoir_level[6] <= 200
 reservoir_level[7] <= 200
 reservoir_level[8] <= 200
 reservoir_level[9] <= 200
 reservoir_level0[0] = 20
 reservoir_level0[1] = 20
 reservoir_level0[2] = 20
 reservoir_level0[3] = 20
 reservoir_level0[4] = 20
 reservoir_level0[5] = 20
 reservoir_level0[6] = 20
 reservoir_level0[7] = 20
 reservoir_level0[8] = 20
 reservoir_level0[9] = 20
 inflow[0,1] = 9.790388199339375
 inflow[1,1] = 8.81232246338927
 inflow[2,1] = 6.276819214639752
 inflow[3,1] = 9.304865334026656
 inflow[4,1] = 7.247899531024389
 inflow[5,1] = 7.166778860553562
 inflow[6,1] = 0.4107856658418385
 inflow[7,1] = 4.394817751674749
 inflow[8,1] = 2.820697831160878
 inflow[9,1] = 3.349959688919273
 inflow0[0,1] = 1.915194503788923
 inflow0[1,1] = 6.221087710398319
 inflow0[2,1] = 4.377277390071145
 inflow0[3,1] = 7.853585837137692
 inflow0[4,1] = 7.799758081188035
 inflow0[5,1] = 2.725926052826416
 inflow0[6,1] = 2.764642551430967
 inflow0[7,1] = 8.018721775350192
 inflow0[8,1] = 9.581393536837052
 inflow0[9,1] = 8.759326347420947
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
 oracle[0][20] >= -1e+08
 oracle[0][21] >= -1e+08
 oracle[0][22] >= -1e+08
 oracle[0][23] >= -1e+08
 oracle[0][24] >= -1e+08
 oracle[0][25] >= -1e+08
 oracle[0][26] >= -1e+08
 oracle[0][27] >= -1e+08
 oracle[0][28] >= -1e+08
 oracle[0][29] >= -1e+08
 lambda[0] free
 nu[0][0] free
 nu[0][1] free
 nu[0][2] free
 nu[0][3] free
 nu[0][4] free
 nu[0][5] free
 nu[0][6] free
 nu[0][7] free
 nu[0][8] free
 nu[0][9] free
 nu[0][10] free
 nu[0][11] free
 nu[0][12] free
 nu[0][13] free
 nu[0][14] free
 nu[0][15] free
 nu[0][16] free
 nu[0][17] free
 nu[0][18] free
 nu[0][19] free
 nu[0][20] free
 nu[0][21] free
 nu[0][22] free
 nu[0][23] free
 nu[0][24] free
 nu[0][25] free
 nu[0][26] free
 nu[0][27] free
 nu[0][28] free
 nu[0][29] free
End
