% Optimizations: vectorized own_csr own_ell single_block
% all kernel runtimes averaged over 50 executions (and measured in ns); PageRank executed with damping factor of 8.500000e-01; Sparse format: CSR; Mode: interleaved with 8192 Bytes; Using thread-buffer;

data_basicblocks = [
% read/write double  &  SpMV ELL double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV ELL 2segment 32 bit  &  SpMV ELL 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV ELL 4segment 16 bit  &  SpMV ELL 4segment 32 bit  &  SpMV ELL 4segment 48 bit  &  SpMV ELL 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% adaptive
227523  6381356  141834  4020872  156737  856151  4610087  7659636  118648  146412  4889803  197806  253385  383099  423670  4304870  6464581  8935732  11118194  190537  210471  215689  218896
% delaunay_n22
127796  17583939  120763  13697371  121528  80493  14613467  24784642  104761  119087  13776860  104270  159205  228158  271026  24082342  44789382  67689224  81669455  127407  135319  140859  143710
% europe_osm
6764493  85917315  3028090  97826277  3564037  9029296  85984992  133723182  1603731  3009877  111716773  1728731  3267319  4935960  6236233  126377797  180900174  210996364  225955276  1335028  1676123  2291233  3023859
% hugebubbles-00020
2522497  20930636  1403228  12193600  1350102  3506868  16833022  29676752  742725  1493944  16728323  795181  1230582  1495801  1991694  21338426  34575020  50654738  58199265  569447  752600  1128766  1501961
% rgg_n_2_24_s0
1651847  125777963  833237  95184619  787971  1336824  95762836  141402120  441039  867259  94359974  435677  759441  1062428  1396488  182716030  248747046  330264563  407380813  462353  525291  593205  722616
% road_usa
3200098  30724919  1271215  33741733  1501840  3934920  21341515  41384571  678405  1263294  39418370  849081  1256659  1886514  2499702  15623368  45817797  60767600  66322876  636215  755231  948257  1229719
% stanford
18420  9476883  74167  9495148  16761  45062  8297503  11649868  73050  73246  9365460  19283  24279  30223  34029  16792386  19272583  21458804  22809226  27622  27163  27463  27627
];
data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% adaptive
20052  497083  19996  20438  293636  534841  19454  20536  20450  20651  194578  244883  349105  435849
% delaunay_n22
23403  174985  23261  23099  121059  186012  22242  22261  22242  22120  129111  155714  215946  259027
% europe_osm
18530  6417977  18308  20570  3211329  6450947  19142  20643  20762  20796  1632764  3261562  4897750  6676422
% hugebubbles-00020
16844  2585739  16331  16138  1471998  2939716  16569  16531  16180  16493  627128  1301372  1973053  2642651
% rgg_n_2_24_s0
26245  2008129  26400  26232  964825  1956957  26833  26967  25867  26228  491156  954871  1407533  1816473
% road_usa
23119  2750605  23229  22676  1339973  2776527  23062  22581  22475  22367  705793  1394978  2111734  2861402
% stanford
21027  76402  21204  20793  73432  80813  20855  20813  21106  21402  26680  28170  32930  36177
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% adaptive
383959  4020872  471018  581973  660449  4889803
% delaunay_n22
262907  13697371  316158  377328  427764  13776860
% europe_osm
9005809  97826277  4815628  7066226  9030644  111716773
% hugebubbles-00020
3190140  12193600  1717112  2326607  3056454  16728323
% rgg_n_2_24_s0
2113853  95184619  1307979  1753746  2252590  94359974
% road_usa
3563042  33741733  1933168  2737698  3624158  39418370
% stanford
95779  9495148  51405  59893  65491  9365460
];
% all PageRank runtimes averaged over 5 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: ELL
data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% adaptive
7380252  5439310  7320549  5175002  5121474  9924290  9572460
% delaunay_n22
19195461  15721859  18469591  24986472  46089275  66648227  43040604
% europe_osm
93130367  93216927  92676548  125354443  184826048  218448866  129323685
% hugebubbles-00020
23881973  19438130  23651045  23045750  24190210  53555954  32310588
% rgg_n_2_24_s0
131876268  101387149  130996718  180630388  188438677  331008521  231368600
% road_usa
34774993  24419907  34338831  17254284  48359877  63361276  45001170
% stanford
9197677  8284675  9248146  16488517  19305522  21562620  13660625
];
data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% adaptive
65  479716431  65  7  472511379  84  1  2  83  835515275
% delaunay_n22
64  1228509552  64  15  1170520140  88  1  22  86  5384314225
% europe_osm
118  10989383339  118  50  11162901066  119  1  51  116  24211788472
% hugebubbles-00020
77  1838911956  77  12  1805397482  82  1  2  79  4314850853
% rgg_n_2_24_s0
89  11736987863  89  28  11028326564  93  1  2  91  30571869777
% road_usa
117  4068674193  117  51  3573503172  117  1  51  115  6645340240
% stanford
118  1085325977  118  51  1060020061  118  1  51  115  2428785383
];
data_pagerank_speedups = [
% double  &  2 segment  &  4 segment
1  1.01525  0.574156
1  1.04954  0.228165
1  0.984456  0.453886
1  1.01856  0.426182
1  1.06426  0.383915
1  1.13857  0.61226
1  1.02387  0.44686
];