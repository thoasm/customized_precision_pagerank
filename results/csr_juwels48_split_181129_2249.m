% Optimizations: vectorized own_csr own_ell
% all kernel runtimes averaged over 50 executions (and measured in ns); PageRank executed with damping factor of 8.500000e-01; Sparse format: CSR; Mode: Using completely separate approach;

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% adaptive
272444  6998372  173511  883892  162322  378854  4704807  8180188  117126  172063  3913408  353698  507780  686717  931662  3441818  5672129  8462853  10932370  134969  181396  231235  290147
% delaunay_n22
123259  6609963  102346  309511  110179  259705  4924957  8243252  98348  121122  2257483  223046  316576  433774  577377  3716394  6502342  9066337  11814886  105916  142048  167606  179370
% europe_osm
7333448  31020773  4099555  10459808  3578163  7875321  26186939  41668984  1948627  3925581  34692278  2501700  3864839  7744361  8851289  20222571  34690491  48789781  55238761  1788668  2679498  5500399  4837310
% hugebubbles-00020
2508366  24266834  1416044  4048348  879490  2735872  19772103  35842637  576424  1415139  13455607  953880  1481163  2212166  3541207  12969344  28665543  40981378  55256170  532040  948155  1286489  2129204
% rgg_n_2_24_s0
1886594  32699305  931005  3242997  644620  1954486  21521487  32663144  412042  1218470  10782221  761045  1149478  1760979  2562842  16119301  22356549  29136187  36315358  292933  639524  942746  1232034
% road_usa
3131158  16248935  1620541  4917429  1197691  3317309  10676236  17381993  699404  1752299  14812183  1065503  1690223  2468987  3501122  7077605  11721506  15773045  20074621  661897  1035577  1426361  1786934
% stanford
22694  371812  72164  15393  18505  30995  327003  422559  71768  71181  244932  25825  33390  44825  55190  333015  426200  520390  614088  74451  76480  77671  77186
% wb-edu
638897  8718293  365027  1625185  278747  1020536  5644885  8731683  180588  435938  5964667  486981  675604  969776  1371897  4263812  6328168  7872048  9488884  224582  259596  361460  535998
% web-BerkStan
31098  1022243  78863  37570  26741  55902  614816  1029921  78597  78281  407828  50028  63225  84472  108848  609589  944410  1251606  1519254  81336  86878  83911  82860
% web-Google
37752  1108013  84366  53633  34390  73081  805527  1328700  82431  80001  485402  57873  90176  117725  146268  725894  1145197  1555554  2001436  87374  90236  85901  88122
];
data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% adaptive
22107  583655  23623  24297  192495  680186  21764  24573  24028  24087  205462  348123  529926  768015
% delaunay_n22
17085  189193  17507  17749  117942  250644  17217  17887  17632  17165  136444  189009  287967  393047
% europe_osm
21907  7994412  22289  23786  4093117  8020144  20576  22883  23000  23267  2709572  3991392  7070458  8402357
% hugebubbles-00020
21980  2999384  21929  23711  1501544  3003409  19647  22309  22232  21926  1884698  2771476  3465414  3975484
% rgg_n_2_24_s0
17801  2234649  17268  17633  1163147  2573316  18256  18093  18527  17717  730867  1321902  2227771  3335985
% road_usa
21066  3458549  21694  23475  1883436  3638648  20501  23507  23269  22375  1158406  1895626  2746870  3832934
% stanford
23818  72930  23830  24465  66999  72767  25174  24131  24398  24362  69546  77554  84412  81759
% wb-edu
124491  1103444  161550  189290  477857  1226895  181766  151568  182696  306093  322843  764066  1297494  1529859
% web-BerkStan
25135  84162  24732  24501  76053  82386  24698  24468  24649  24615  79736  83664  92336  102774
% web-Google
28980  89023  28977  29579  78795  86391  29270  28593  30346  30463  77091  94204  110797  118067
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% adaptive
1874066  883892  1269975  1664848  2091861  3913408
% delaunay_n22
1140780  309511  773105  1020684  1270645  2257483
% europe_osm
12379288  10459808  8412298  13792204  15197781  34692278
% hugebubbles-00020
5310358  4048348  3675749  4686383  6314179  13455607
% rgg_n_2_24_s0
4088616  3242997  3099203  3931905  5019171  10782221
% road_usa
6068462  4917429  3953867  5151741  6486167  14812183
% stanford
138800  15393  102235  109862  119674  244932
% wb-edu
2641022  1625185  1847944  2402560  3071723  5964667
% web-BerkStan
241412  37570  154829  181364  220553  407828
% web-Google
289398  53633  186538  259069  286064  485402
];
% all PageRank runtimes averaged over 5 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR
data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% adaptive
8164665  5145014  8000688  4072264  4051995  8966152  14324599
% delaunay_n22
7326225  5327233  7195014  4092048  6739187  9597853  8121479
% europe_osm
40123171  33188202  39822937  23546355  42225459  57032303  38103117
% hugebubbles-00020
28091989  22742138  27764792  13978976  16771517  46821687  31856527
% rgg_n_2_24_s0
36058733  23170225  35649372  16901308  18255191  31284788  30839056
% road_usa
20456636  14176898  20028829  8139459  14949877  20198372  16630375
% stanford
496095  415250  496794  448768  513609  599086  513765
% wb-edu
11186287  7509270  11136917  5289763  8244457  10420052  10013459
% web-BerkStan
1221735  740830  1205927  720569  1064097  1081436  1224336
% web-Google
1284172  900848  1264946  844403  1238947  1677818  1422377
];
data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% adaptive
65  530703250  65  7  507958400  84  1  2  83  761719879
% delaunay_n22
64  468878455  64  15  439242079  88  1  22  86  785534794
% europe_osm
118  4734534181  118  50  4423397506  119  1  51  116  6051869894
% hugebubbles-00020
77  2163083157  77  12  2109718487  82  1  2  79  3773701974
% rgg_n_2_24_s0
89  3209227281  89  28  2853880584  93  1  2  91  2920915367
% road_usa
117  2393426527  117  51  2070087673  117  1  51  115  2120133906
% stanford
118  58539327  118  51  55032831  118  1  51  115  67037722
% wb-edu
114  1275236758  114  47  1110885072  115  1  47  112  1110455542
% web-BerkStan
119  145386544  119  50  121270944  119  1  50  51  138883457
% web-Google
116  148964063  116  47  130865572  116  1  47  114  175156310
];
data_pagerank_speedups = [
% double  &  2 segment  &  4 segment
1  1.04478  0.696717
1  1.06747  0.596891
1  1.07034  0.782326
1  1.02529  0.573199
1  1.12451  1.09871
1  1.1562  1.1289
1  1.06372  0.87323
1  1.14795  1.14839
1  1.19886  1.04682
1  1.1383  0.850464
];
