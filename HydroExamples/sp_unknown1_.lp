\ Model Hydrovalley_1
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  1000 pour[0] + 1000 pour[1] - 5.46 generation + 10 thermal_gen
   + 0.01 gamma[1] + 0.2 nu[1][0] + 0.2 nu[1][1] + 0.2 nu[1][2]
   + 0.2 nu[1][3] + 0.2 nu[1][4]
Subject To
 AR_model_1[0]: inflow[0,1] - 0.0637258834012862 inflow0[0,1]
   - innovations[0] = 0
 AR_model_1[1]: inflow[1,1] - 0.0640446838459459 inflow0[0,1]
   - 0.0327539552224049 inflow0[1,1] - innovations[1] = 0
 balance[0]: reservoir_level[0] - reservoir_level0[0] - inflow[0,1]
   + outflow[0] + spill[0] - pour[0] = 0
 balance[1]: reservoir_level[1] - reservoir_level0[1] - inflow[1,1]
   - outflow[0] + outflow[1] - spill[0] + spill[1] - pour[1] = 0
 generationCtr: generation - 55 dispatch[0,0] - 65 dispatch[0,1]
   - 70 dispatch[0,2] - 55 dispatch[1,0] - 65 dispatch[1,1]
   - 70 dispatch[1,2] = 0
 R5: generation + thermal_gen >= 500
 outflowCtr[0]: outflow[0] - 50 dispatch[0,0] - 60 dispatch[0,1]
   - 70 dispatch[0,2] = 0
 outflowCtr[1]: outflow[1] - 50 dispatch[1,0] - 60 dispatch[1,1]
   - 70 dispatch[1,2] = 0
 dispatchCtr[0]: dispatch[0,0] + dispatch[0,1] + dispatch[0,2] <= 1
 dispatchCtr[1]: dispatch[1,0] + dispatch[1,1] + dispatch[1,2] <= 1
 dro_dual_ctr[00]: - oracle[1][0] + nu[1][0] >= 0
 dro_dual_ctr[01]: - oracle[1][1] + 1.449569101930528 gamma[1] + nu[1][0]
   >= 0
 dro_dual_ctr[02]: - oracle[1][2] + 1.389303930123527 gamma[1] + nu[1][0]
   >= 0
 dro_dual_ctr[03]: - oracle[1][3] + 5.001800074637798 gamma[1] + nu[1][0]
   >= 0
 dro_dual_ctr[04]: - oracle[1][4] + 2.455705837660444 gamma[1] + nu[1][0]
   >= 0
 dro_dual_ctr[10]: - oracle[1][0] + 1.449569101930528 gamma[1] + nu[1][1]
   >= 0
 dro_dual_ctr[11]: - oracle[1][1] + nu[1][1] >= 0
 dro_dual_ctr[12]: - oracle[1][2] + 2.766882820714683 gamma[1] + nu[1][1]
   >= 0
 dro_dual_ctr[13]: - oracle[1][3] + 6.379378965228954 gamma[1] + nu[1][1]
   >= 0
 dro_dual_ctr[14]: - oracle[1][4] + 3.032791554542609 gamma[1] + nu[1][1]
   >= 0
 dro_dual_ctr[20]: - oracle[1][0] + 1.389303930123527 gamma[1] + nu[1][2]
   >= 0
 dro_dual_ctr[21]: - oracle[1][1] + 2.766882820714683 gamma[1] + nu[1][2]
   >= 0
 dro_dual_ctr[22]: - oracle[1][2] + nu[1][2] >= 0
 dro_dual_ctr[23]: - oracle[1][3] + 3.612496144514271 gamma[1] + nu[1][2]
   >= 0
 dro_dual_ctr[24]: - oracle[1][4] + 1.521367653891211 gamma[1] + nu[1][2]
   >= 0
 dro_dual_ctr[30]: - oracle[1][0] + 5.001800074637798 gamma[1] + nu[1][3]
   >= 0
 dro_dual_ctr[31]: - oracle[1][1] + 6.379378965228954 gamma[1] + nu[1][3]
   >= 0
 dro_dual_ctr[32]: - oracle[1][2] + 3.612496144514271 gamma[1] + nu[1][3]
   >= 0
 dro_dual_ctr[33]: - oracle[1][3] + nu[1][3] >= 0
 dro_dual_ctr[34]: - oracle[1][4] + 3.346587410686345 gamma[1] + nu[1][3]
   >= 0
 dro_dual_ctr[40]: - oracle[1][0] + 2.455705837660444 gamma[1] + nu[1][4]
   >= 0
 dro_dual_ctr[41]: - oracle[1][1] + 3.032791554542609 gamma[1] + nu[1][4]
   >= 0
 dro_dual_ctr[42]: - oracle[1][2] + 1.521367653891211 gamma[1] + nu[1][4]
   >= 0
 dro_dual_ctr[43]: - oracle[1][3] + 3.346587410686345 gamma[1] + nu[1][4]
   >= 0
 dro_dual_ctr[44]: - oracle[1][4] + nu[1][4] >= 0
 cut[1,0,0]: 54.912 reservoir_level[0] + 27.456 reservoir_level[1]
   + 13.60651547586231 inflow[0,1] + 5.530631076459557 inflow[1,1]
   + oracle[1][0] >= 49675.21186678657
 cut[1,0,1]: 54.912 reservoir_level[0] + 27.456 reservoir_level[1]
   + 13.60651547586231 inflow[0,1] + 5.530631076459557 inflow[1,1]
   + oracle[1][1] >= 49636.69126440095
 cut[1,0,2]: 54.912 reservoir_level[0] + 27.456 reservoir_level[1]
   + 13.60651547586231 inflow[0,1] + 5.530631076459557 inflow[1,1]
   + oracle[1][2] >= 49751.19098334396
 cut[1,0,3]: 54.912 reservoir_level[0] + 27.456 reservoir_level[1]
   + 13.60651547586231 inflow[0,1] + 5.530631076459557 inflow[1,1]
   + oracle[1][3] >= 49911.53618640096
 cut[1,0,4]: 54.912 reservoir_level[0] + 27.456 reservoir_level[1]
   + 13.60651547586231 inflow[0,1] + 5.530631076459557 inflow[1,1]
   + oracle[1][4] >= 49786.53626024566
 cut[1,1,0]: 54.33226620904721 reservoir_level[0]
   + 26.994 reservoir_level[1] + 13.43537781243553 inflow[0,1]
   + 5.437567572540843 inflow[1,1] + oracle[1][0] >= 49620.39678120045
 cut[1,1,1]: 54.33226620904721 reservoir_level[0]
   + 26.994 reservoir_level[1] + 13.43537781243553 inflow[0,1]
   + 5.437567572540843 inflow[1,1] + oracle[1][1] >= 49582.53725587865
 cut[1,1,2]: 54.33226620904721 reservoir_level[0]
   + 26.994 reservoir_level[1] + 13.43537781243553 inflow[0,1]
   + 5.437567572540843 inflow[1,1] + oracle[1][2] >= 49695.5135802881
 cut[1,1,3]: 54.33226620904721 reservoir_level[0]
   + 26.994 reservoir_level[1] + 13.43537781243553 inflow[0,1]
   + 5.437567572540843 inflow[1,1] + oracle[1][3] >= 49853.81263791311
 cut[1,1,4]: 54.33226620904721 reservoir_level[0]
   + 26.994 reservoir_level[1] + 13.43537781243553 inflow[0,1]
   + 5.437567572540843 inflow[1,1] + oracle[1][4] >= 49730.58421611808
 cut[1,2,0]: 53.988 reservoir_level[0] + 26.994 reservoir_level[1]
   + 13.37755968724097 inflow[0,1] + 5.43756757311384 inflow[1,1]
   + oracle[1][0] >= 49610.8361385367
 cut[1,2,1]: 53.988 reservoir_level[0] + 26.994 reservoir_level[1]
   + 13.37755968724097 inflow[0,1] + 5.43756757311384 inflow[1,1]
   + oracle[1][1] >= 49572.96371936184
 cut[1,2,2]: 53.988 reservoir_level[0] + 26.994 reservoir_level[1]
   + 13.37755968724097 inflow[0,1] + 5.43756757311384 inflow[1,1]
   + oracle[1][2] >= 49685.53676034753
 cut[1,2,3]: 53.988 reservoir_level[0] + 26.994 reservoir_level[1]
   + 13.37755968724097 inflow[0,1] + 5.43756757311384 inflow[1,1]
   + oracle[1][3] >= 49843.18384701453
 cut[1,2,4]: 53.988 reservoir_level[0] + 26.994 reservoir_level[1]
   + 13.37755968724097 inflow[0,1] + 5.43756757311384 inflow[1,1]
   + oracle[1][4] >= 49720.28728499468
 cut[1,3,0]: 53.53337572132897 reservoir_level[0]
   + 26.53937572132897 reservoir_level[1] + 13.22930931340951 inflow[0,1]
   + 5.346464840784465 inflow[1,1] + oracle[1][0] >= 49559.24280986878
 cut[1,3,1]: 53.53337572132897 reservoir_level[0]
   + 26.53937572132897 reservoir_level[1] + 13.22930931340951 inflow[0,1]
   + 5.346464840784465 inflow[1,1] + oracle[1][1] >= 49522.02183933045
 cut[1,3,2]: 53.53337572132897 reservoir_level[0]
   + 26.53937572132897 reservoir_level[1] + 13.22930931340951 inflow[0,1]
   + 5.346464840784465 inflow[1,1] + oracle[1][2] >= 49633.23807200013
 cut[1,3,3]: 53.53337572132897 reservoir_level[0]
   + 26.53937572132897 reservoir_level[1] + 13.22930931340951 inflow[0,1]
   + 5.346464840784465 inflow[1,1] + oracle[1][3] >= 49789.09954729985
 cut[1,3,4]: 53.53337572132897 reservoir_level[0]
   + 26.53937572132897 reservoir_level[1] + 13.22930931340951 inflow[0,1]
   + 5.346464840784465 inflow[1,1] + oracle[1][4] >= 49667.82649476841
 cut[1,4,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10727823976897 inflow[0,1] + 5.271434831063891 inflow[1,1]
   + oracle[1][0] >= 49510.85923002064
 cut[1,4,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10727823976897 inflow[0,1] + 5.271434831063891 inflow[1,1]
   + oracle[1][1] >= 49474.17479119456
 cut[1,4,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10727823976897 inflow[0,1] + 5.271434831063891 inflow[1,1]
   + oracle[1][2] >= 49584.27403817524
 cut[1,4,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10727823976897 inflow[0,1] + 5.271434831063891 inflow[1,1]
   + oracle[1][3] >= 49738.66565274705
 cut[1,4,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10727823976897 inflow[0,1] + 5.271434831063891 inflow[1,1]
   + oracle[1][4] >= 49618.72931404442
 cut[1,5,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10727823976897 inflow[0,1] + 5.271434831063891 inflow[1,1]
   + oracle[1][0] >= 49510.85923002064
 cut[1,5,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10727823976897 inflow[0,1] + 5.271434831063891 inflow[1,1]
   + oracle[1][1] >= 49474.17479119455
 cut[1,5,2]: 53.57441538993874 reservoir_level[0]
   + 26.58041538993874 reservoir_level[1] + 13.24260445590609 inflow[0,1]
   + 5.35463763306698 inflow[1,1] + oracle[1][2] >= 49639.4664377945
 cut[1,5,3]: 53.82046852307567 reservoir_level[0]
   + 26.82646852307566 reservoir_level[1] + 13.32353832539762 inflow[0,1]
   + 5.404398314410672 inflow[1,1] + oracle[1][3] >= 49829.16355134212
 cut[1,5,4]: 53.57441538993874 reservoir_level[0]
   + 26.58041538993874 reservoir_level[1] + 13.24260445590609 inflow[0,1]
   + 5.35463763306698 inflow[1,1] + oracle[1][4] >= 49674.06937861436
 cut[1,6,0]: 53.24494397581032 reservoir_level[0]
   + 26.25094397581032 reservoir_level[1] + 13.13412749245813 inflow[0,1]
   + 5.287949502938988 inflow[1,1] + oracle[1][0] >= 49523.60530012063
 cut[1,6,1]: 53.24494397581032 reservoir_level[0]
   + 26.25094397581032 reservoir_level[1] + 13.13412749245813 inflow[0,1]
   + 5.287949502938988 inflow[1,1] + oracle[1][1] >= 49486.80276422397
 cut[1,6,2]: 53.24494397581032 reservoir_level[0]
   + 26.25094397581032 reservoir_level[1] + 13.13412749245813 inflow[0,1]
   + 5.287949502938988 inflow[1,1] + oracle[1][2] >= 49597.1477935006
 cut[1,6,3]: 53.24494397581032 reservoir_level[0]
   + 26.25094397581032 reservoir_level[1] + 13.13412749245813 inflow[0,1]
   + 5.287949502938988 inflow[1,1] + oracle[1][3] >= 49751.86281477096
 cut[1,6,4]: 53.24494397581032 reservoir_level[0]
   + 26.25094397581032 reservoir_level[1] + 13.13412749245813 inflow[0,1]
   + 5.287949502938988 inflow[1,1] + oracle[1][4] >= 49631.63231672088
 cut[1,7,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49513.43514643209
 cut[1,7,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49476.75064649891
 cut[1,7,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49586.84993647868
 cut[1,7,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][3] >= 49741.24158388683
 cut[1,7,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49621.30516420903
 cut[1,8,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][0] >= 49513.69244655768
 cut[1,8,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][1] >= 49477.00794662449
 cut[1,8,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][2] >= 49587.10723660428
 cut[1,8,3]: 53.29449370461514 reservoir_level[0]
   + 26.30049370461513 reservoir_level[1] + 13.15053185980782 inflow[0,1]
   + 5.298034081599159 inflow[1,1] + oracle[1][3] >= 49759.14087608092
 cut[1,8,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][4] >= 49621.56246433462
 cut[1,9,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49513.96446456941
 cut[1,9,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.27996463624
 cut[1,9,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.379254616
 cut[1,9,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][3] >= 49741.77090202415
 cut[1,9,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49621.83448234635
 cut[1,10,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.00931751145
 cut[1,10,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.32481757827
 cut[1,10,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.42410755803
 cut[1,10,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][3] >= 49741.81575496616
 cut[1,10,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49621.87933528838
 cut[1,11,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.03361596599
 cut[1,11,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.34911603281
 cut[1,11,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.44840601258
 cut[1,11,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][3] >= 49741.84005342072
 cut[1,11,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49621.90363374291
 cut[1,12,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][0] >= 49514.07215404622
 cut[1,12,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][1] >= 49477.38765411304
 cut[1,12,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][2] >= 49587.48694409281
 cut[1,12,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][3] >= 49741.87859150095
 cut[1,12,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][4] >= 49621.94217182317
 cut[1,13,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][0] >= 49514.09810732028
 cut[1,13,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][1] >= 49477.4136073871
 cut[1,13,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][2] >= 49587.51289736685
 cut[1,13,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][3] >= 49741.90454477499
 cut[1,13,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.27144302339628 inflow[1,1]
   + oracle[1][4] >= 49621.96812509721
 cut[1,14,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.11339141578
 cut[1,14,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.4288914826
 cut[1,14,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.52818146239
 cut[1,14,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][3] >= 49741.91982887052
 cut[1,14,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49621.98340919273
 cut[1,15,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.12762967558
 cut[1,15,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.44312974239
 cut[1,15,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.54241972216
 cut[1,15,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][3] >= 49741.93406713029
 cut[1,15,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49621.99764745252
 cut[1,16,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.1340902251
 cut[1,16,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.44959029194
 cut[1,16,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.54888027169
 cut[1,16,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][3] >= 49741.94052767984
 cut[1,16,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49622.00410800205
 cut[1,17,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.13409022511
 cut[1,17,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.44959029193
 cut[1,17,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.5488802717
 cut[1,17,3]: 53.4193550721597 reservoir_level[0]
   + 26.42535507215971 reservoir_level[1] + 13.19174250495646 inflow[0,1]
   + 5.323369276394431 inflow[1,1] + oracle[1][3] >= 49776.59512392159
 cut[1,17,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49622.00410800205
 cut[1,18,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.13485130554
 cut[1,18,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.45035137236
 cut[1,18,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.54964135213
 cut[1,18,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][3] >= 49741.94128876027
 cut[1,18,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49622.00486908248
 cut[1,19,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.1480680368
 cut[1,19,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.46356810362
 cut[1,19,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.56285808339
 cut[1,19,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][3] >= 49741.95450549152
 cut[1,19,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49622.01808581375
 cut[1,20,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.15568449118
 cut[1,20,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.471184558
 cut[1,20,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.57047453777
 cut[1,20,3]: 53.4193550721597 reservoir_level[0]
   + 26.42535507215971 reservoir_level[1] + 13.19174250495646 inflow[0,1]
   + 5.323369276394431 inflow[1,1] + oracle[1][3] >= 49776.59512392159
 cut[1,20,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49622.02570226813
 cut[1,21,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.15821763237
 cut[1,21,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.47371769919
 cut[1,21,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.57300767896
 cut[1,21,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][3] >= 49741.9646550871
 cut[1,21,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49622.02823540932
 cut[1,22,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.1588359513
 cut[1,22,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.47433601812
 cut[1,22,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.57362599788
 cut[1,22,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][3] >= 49741.96527340602
 cut[1,22,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49622.02885372824
 cut[1,23,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.1588359513
 cut[1,23,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.47433601812
 cut[1,23,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.57362599788
 cut[1,23,3]: 53.19278027964091 reservoir_level[0]
   + 26.19878027964091 reservoir_level[1] + 13.11709206319125 inflow[0,1]
   + 5.27747519713859 inflow[1,1] + oracle[1][3] >= 49746.02338250267
 cut[1,23,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49622.02885372824
 cut[1,24,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][0] >= 49514.16333440738
 cut[1,24,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][1] >= 49477.4788344742
 cut[1,24,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][2] >= 49587.57812445396
 cut[1,24,3]: 53.19278027964091 reservoir_level[0]
   + 26.19878027964091 reservoir_level[1] + 13.11709206319125 inflow[0,1]
   + 5.27747519713859 inflow[1,1] + oracle[1][3] >= 49746.0260815763
 cut[1,24,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728024341037 inflow[0,1] + 5.271443023396283 inflow[1,1]
   + oracle[1][4] >= 49622.03335218432
 cut[1,25,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.16683238281
 cut[1,25,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.48233259011
 cut[1,25,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.58162247103
 cut[1,25,3]: 53.18206859392648 reservoir_level[0]
   + 26.18806859392647 reservoir_level[1] + 13.11356284047122 inflow[0,1]
   + 5.275305466069222 inflow[1,1] + oracle[1][3] >= 49744.5718020388
 cut[1,25,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.03685031204
 cut[1,26,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.1668323828
 cut[1,26,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.48233259009
 cut[1,26,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.58162247101
 cut[1,26,3]: 53.46607636735672 reservoir_level[0]
   + 26.47207636735672 reservoir_level[1] + 13.20713594430392 inflow[0,1]
   + 5.332832954153419 inflow[1,1] + oracle[1][3] >= 49782.94080273609
 cut[1,26,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.03685031203
 cut[1,27,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.16935761252
 cut[1,27,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.48485781981
 cut[1,27,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.58414770073
 cut[1,27,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.97579503337
 cut[1,27,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.03937554174
 cut[1,28,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][0] >= 49514.17117807125
 cut[1,28,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][1] >= 49477.48667827857
 cut[1,28,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][2] >= 49587.58596815947
 cut[1,28,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][3] >= 49741.97761549213
 cut[1,28,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][4] >= 49622.04119600049
 cut[1,29,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17178644968
 cut[1,29,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.487286657
 cut[1,29,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.58657653791
 cut[1,29,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.97822387056
 cut[1,29,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04180437892
 cut[1,30,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][0] >= 49514.17179523739
 cut[1,30,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][1] >= 49477.48729544469
 cut[1,30,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][2] >= 49587.5865853256
 cut[1,30,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][3] >= 49741.97823265825
 cut[1,30,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][4] >= 49622.04181316661
 cut[1,31,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17373911592
 cut[1,31,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.48923932323
 cut[1,31,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.58852920415
 cut[1,31,3]: 53.18206859392648 reservoir_level[0]
   + 26.18806859392647 reservoir_level[1] + 13.11356284047122 inflow[0,1]
   + 5.275305466069222 inflow[1,1] + oracle[1][3] >= 49744.5773274253
 cut[1,31,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04375704515
 cut[1,32,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17379845416
 cut[1,32,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.48929866147
 cut[1,32,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.58858854238
 cut[1,32,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98023587504
 cut[1,32,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.0438163834
 cut[1,33,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17508827168
 cut[1,33,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49058847898
 cut[1,33,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.5898783599
 cut[1,33,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98152569256
 cut[1,33,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04510620092
 cut[1,34,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][0] >= 49514.17542550128
 cut[1,34,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][1] >= 49477.49092570857
 cut[1,34,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][2] >= 49587.59021558949
 cut[1,34,3]: 53.46607636735672 reservoir_level[0]
   + 26.47207636735672 reservoir_level[1] + 13.20713594430392 inflow[0,1]
   + 5.332832954153419 inflow[1,1] + oracle[1][3] >= 49782.94080273609
 cut[1,34,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][4] >= 49622.04544343051
 cut[1,35,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17551727735
 cut[1,35,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49101748467
 cut[1,35,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59030736558
 cut[1,35,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98195469823
 cut[1,35,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04553520658
 cut[1,36,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17556607965
 cut[1,36,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49106628694
 cut[1,36,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59035616786
 cut[1,36,3]: 53.18206859392648 reservoir_level[0]
   + 26.18806859392647 reservoir_level[1] + 13.11356284047122 inflow[0,1]
   + 5.275305466069222 inflow[1,1] + oracle[1][3] >= 49744.57878899628
 cut[1,36,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04558400888
 cut[1,37,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17633550536
 cut[1,37,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49183571265
 cut[1,37,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59112559356
 cut[1,37,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98277292623
 cut[1,37,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04635343459
 cut[1,38,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17633550535
 cut[1,38,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49183571264
 cut[1,38,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59112559356
 cut[1,38,3]: 53.55902487858334 reservoir_level[0]
   + 26.56502487858334 reservoir_level[1] + 13.23772463917251 inflow[0,1]
   + 5.351638739430357 inflow[1,1] + oracle[1][3] >= 49795.44989620717
 cut[1,38,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04635343458
 cut[1,39,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][0] >= 49514.17702201
 cut[1,39,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][1] >= 49477.4925222173
 cut[1,39,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][2] >= 49587.59181209822
 cut[1,39,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][3] >= 49741.98345943088
 cut[1,39,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][4] >= 49622.04703993924
 cut[1,40,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17815819487
 cut[1,40,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49365840217
 cut[1,40,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59294828309
 cut[1,40,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98459561574
 cut[1,40,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.0481761241
 cut[1,41,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][0] >= 49514.17840518338
 cut[1,41,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][1] >= 49477.49390539069
 cut[1,41,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][2] >= 49587.59319527159
 cut[1,41,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][3] >= 49741.98484260424
 cut[1,41,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][4] >= 49622.0484231126
 cut[1,42,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][0] >= 49514.17856156772
 cut[1,42,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][1] >= 49477.49406177501
 cut[1,42,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][2] >= 49587.59335165593
 cut[1,42,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][3] >= 49741.98499898858
 cut[1,42,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][4] >= 49622.04857949694
 cut[1,43,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][0] >= 49514.1785626682
 cut[1,43,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][1] >= 49477.4940628755
 cut[1,43,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][2] >= 49587.59335275643
 cut[1,43,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][3] >= 49741.98500008907
 cut[1,43,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][4] >= 49622.04858059743
 cut[1,44,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17867643965
 cut[1,44,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49417664695
 cut[1,44,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59346652785
 cut[1,44,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.9851138605
 cut[1,44,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04869436887
 cut[1,45,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.1787236034
 cut[1,45,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.4942238107
 cut[1,45,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59351369163
 cut[1,45,3]: 53.20337943767198 reservoir_level[0]
   + 26.20937943767198 reservoir_level[1] + 13.12058420362599 inflow[0,1]
   + 5.279622104919433 inflow[1,1] + oracle[1][3] >= 49747.4815574789
 cut[1,45,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04874153264
 cut[1,46,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][0] >= 49514.17883041345
 cut[1,46,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][1] >= 49477.49433062075
 cut[1,46,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][2] >= 49587.59362050165
 cut[1,46,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][3] >= 49741.98526783431
 cut[1,46,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562975 inflow[1,1]
   + oracle[1][4] >= 49622.04884834267
 cut[1,47,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17896891163
 cut[1,47,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49446911894
 cut[1,47,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59375899986
 cut[1,47,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.9854063325
 cut[1,47,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04898684086
 cut[1,48,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17909695833
 cut[1,48,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49459716563
 cut[1,48,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59388704655
 cut[1,48,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98553437921
 cut[1,48,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04911488757
 cut[1,49,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17911426609
 cut[1,49,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.4946144734
 cut[1,49,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59390435431
 cut[1,49,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98555168696
 cut[1,49,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04913219532
 cut[1,50,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17911426608
 cut[1,50,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49461447339
 cut[1,50,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.5939043543
 cut[1,50,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98555168696
 cut[1,50,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04913219532
 cut[1,51,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17915157615
 cut[1,51,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49465178345
 cut[1,51,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59394166437
 cut[1,51,3]: 53.55902487858334 reservoir_level[0]
   + 26.56502487858334 reservoir_level[1] + 13.23772463917251 inflow[0,1]
   + 5.351638739430357 inflow[1,1] + oracle[1][3] >= 49795.44989620717
 cut[1,51,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04916950539
 cut[1,52,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17924673036
 cut[1,52,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49474693767
 cut[1,52,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59403681859
 cut[1,52,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98568415124
 cut[1,52,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.0492646596
 cut[1,53,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17942141154
 cut[1,53,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49492161885
 cut[1,53,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59421149977
 cut[1,53,3]: 53.55902487858334 reservoir_level[0]
   + 26.56502487858334 reservoir_level[1] + 13.23772463917251 inflow[0,1]
   + 5.351638739430357 inflow[1,1] + oracle[1][3] >= 49795.44989620717
 cut[1,53,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04943934078
 cut[1,54,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17946912349
 cut[1,54,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49496933079
 cut[1,54,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59425921171
 cut[1,54,3]: 53.55902487858334 reservoir_level[0]
   + 26.56502487858334 reservoir_level[1] + 13.23772463917251 inflow[0,1]
   + 5.351638739430357 inflow[1,1] + oracle[1][3] >= 49795.44989620717
 cut[1,54,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04948705273
 cut[1,55,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17947111485
 cut[1,55,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49497132216
 cut[1,55,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59426120307
 cut[1,55,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98590853572
 cut[1,55,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04948904408
 cut[1,56,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17958780792
 cut[1,56,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49508801521
 cut[1,56,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59437789613
 cut[1,56,3]: 53.55902487858334 reservoir_level[0]
   + 26.56502487858334 reservoir_level[1] + 13.23772463917251 inflow[0,1]
   + 5.351638739430357 inflow[1,1] + oracle[1][3] >= 49795.44989620717
 cut[1,56,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04960573716
 cut[1,57,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562976 inflow[1,1]
   + oracle[1][0] >= 49514.17959477725
 cut[1,57,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562976 inflow[1,1]
   + oracle[1][1] >= 49477.49509498454
 cut[1,57,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562976 inflow[1,1]
   + oracle[1][2] >= 49587.59438486546
 cut[1,57,3]: 53.20191639859922 reservoir_level[0]
   + 26.20791639859922 reservoir_level[1] + 13.12010217068716 inflow[0,1]
   + 5.279325757605536 inflow[1,1] + oracle[1][3] >= 49747.28322587473
 cut[1,57,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562976 inflow[1,1]
   + oracle[1][4] >= 49622.04961270648
 cut[1,58,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17971140124
 cut[1,58,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49521160854
 cut[1,58,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59450148946
 cut[1,58,3]: 53.55902487858334 reservoir_level[0]
   + 26.56502487858334 reservoir_level[1] + 13.23772463917251 inflow[0,1]
   + 5.351638739430357 inflow[1,1] + oracle[1][3] >= 49795.44989620717
 cut[1,58,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04972933047
 cut[1,59,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.17971839848
 cut[1,59,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49521860579
 cut[1,59,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.5945084867
 cut[1,59,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98615581936
 cut[1,59,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.04973632772
 cut[1,60,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.1802226041
 cut[1,60,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49572281141
 cut[1,60,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59501269232
 cut[1,60,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98666002497
 cut[1,60,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.05024053333
 cut[1,61,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.18031103139
 cut[1,61,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.4958112387
 cut[1,61,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59510111962
 cut[1,61,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98674845227
 cut[1,61,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.05032896063
 cut[1,62,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.18034946052
 cut[1,62,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49584966782
 cut[1,62,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59513954874
 cut[1,62,3]: 53.55902487858334 reservoir_level[0]
   + 26.56502487858334 reservoir_level[1] + 13.23772463917251 inflow[0,1]
   + 5.351638739430357 inflow[1,1] + oracle[1][3] >= 49795.44989620717
 cut[1,62,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.05036738975
 cut[1,63,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.18035073167
 cut[1,63,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49585093898
 cut[1,63,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59514081989
 cut[1,63,3]: 53.55902487858334 reservoir_level[0]
   + 26.56502487858334 reservoir_level[1] + 13.23772463917251 inflow[0,1]
   + 5.351638739430357 inflow[1,1] + oracle[1][3] >= 49795.44989620717
 cut[1,63,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.05036866091
 cut[1,64,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.18035073168
 cut[1,64,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49585093898
 cut[1,64,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59514081991
 cut[1,64,3]: 53.19902846263717 reservoir_level[0]
   + 26.20502846263717 reservoir_level[1] + 13.11915067165693 inflow[0,1]
   + 5.278740788905113 inflow[1,1] + oracle[1][3] >= 49746.8917398332
 cut[1,64,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.05036866092
 cut[1,65,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.18039945066
 cut[1,65,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49589965796
 cut[1,65,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59518953886
 cut[1,65,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.9868368715
 cut[1,65,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.05041737988
 cut[1,66,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.18039945064
 cut[1,66,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49589965794
 cut[1,66,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59518953885
 cut[1,66,3]: 53.55902487858334 reservoir_level[0]
   + 26.56502487858334 reservoir_level[1] + 13.23772463917251 inflow[0,1]
   + 5.351638739430357 inflow[1,1] + oracle[1][3] >= 49795.44989620717
 cut[1,66,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.05041737988
 cut[1,67,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.1804231896
 cut[1,67,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.4959233969
 cut[1,67,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.59521327782
 cut[1,67,3]: 53.19902846263717 reservoir_level[0]
   + 26.20502846263717 reservoir_level[1] + 13.11915067165693 inflow[0,1]
   + 5.278740788905113 inflow[1,1] + oracle[1][3] >= 49746.89178330795
 cut[1,67,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.05044111884
 cut[1,68,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.18067213899
 cut[1,68,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.4961723463
 cut[1,68,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.5954622272
 cut[1,68,3]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][3] >= 49741.98710955984
 cut[1,68,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.05069006822
 cut[1,69,0]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][0] >= 49514.18067213898
 cut[1,69,1]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][1] >= 49477.49617234628
 cut[1,69,2]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][2] >= 49587.5954622272
 cut[1,69,3]: 53.56577951312128 reservoir_level[0]
   + 26.57177951312128 reservoir_level[1] + 13.23994714246457 inflow[0,1]
   + 5.353005127019888 inflow[1,1] + oracle[1][3] >= 49796.38332897174
 cut[1,69,4]: 53.163 reservoir_level[0] + 26.169 reservoir_level[1]
   + 13.10728023880421 inflow[0,1] + 5.271443004562977 inflow[1,1]
   + oracle[1][4] >= 49622.05069006822
Bounds
 30 <= reservoir_level[0] <= 200
 30 <= reservoir_level[1] <= 200
 reservoir_level0[0] = 63.03691084341838
 reservoir_level0[1] = 56.54260135843519
 inflow[0,1] free
 inflow[1,1] free
 inflow0[0,1] = 13.03691084341839
 inflow0[1,1] = 6.54260135843519
 innovations[0] = 10.9901797184652
 innovations[1] = 5.058349588095452
 dispatch[0,0] <= 1
 dispatch[0,1] <= 1
 dispatch[0,2] <= 1
 dispatch[1,0] <= 1
 dispatch[1,1] <= 1
 dispatch[1,2] <= 1
 oracle[1][0] >= -1e+08
 oracle[1][1] >= -1e+08
 oracle[1][2] >= -1e+08
 oracle[1][3] >= -1e+08
 oracle[1][4] >= -1e+08
 nu[1][0] free
 nu[1][1] free
 nu[1][2] free
 nu[1][3] free
 nu[1][4] free
End
