% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double    SpMV ELL double    norm double    inplace conversion 2segment    read/write 2segment 32 bit    read/write 2segment 64 bit    SpMV ELL 2segment 32 bit    SpMV ELL 2segment 64 bit    norm 2segment 32 bit    norm 2segment 64 bit    inplace conversion 4segment    read/write 4segment 16 bit    read/write 4segment 32 bit    read/write 4segment 48 bit    read/write 4segment 64 bit    SpMV ELL 4segment 16 bit    SpMV ELL 4segment 32 bit    SpMV ELL 4segment 48 bit    SpMV ELL 4segment 64 bit    norm 4segment 16 bit    norm 4segment 32 bit    norm 4segment 48 bit    norm 4segment 64 bit
% adaptive
202762  950818  120374  1012092  101928  203865  648089  1092040  79537  127587  1036219  67104  105190  156611  211929  497937  797266  1193496  1742624  70899  88870  125369  152635
% delaunay_n22
125853  1878033  85906  2236565  63590  126203  1283905  2101986  62923  87564  2288230  42217  65666  97285  131480  1065030  1667562  2387418  3007058  54826  62967  87470  101428
% europe_osm
1505556  16211730  736970  19574112  748070  1507975  11764610  18699619  464275  767265  20014017  483721  768849  1153950  1565106  9806406  14915230  21235959  25144883  389530  529762  861675  1156489
% hugebubbles-00020
625896  4095283  322404  2511609  312691  629170  3272966  6274098  202195  333350  2569350  202546  321646  481863  653278  2674614  5241578  8962656  11748204  169916  235677  357782  477875
% rgg_n_2_24_s0
497436  11262552  259486  13893814  248060  498286  7684573  12773533  163619  268564  14207045  160892  255075  381812  517581  6588388  10242985  14709292  20602441  139546  191507  283907  371380
% road_usa
708972  5466077  358248  6377011  353102  710340  3766493  6211309  226817  374046  6521265  228629  362947  544108  737439  3026738  4682680  6526868  7973439  191155  265390  399515  540133
% Stanford
4764  1806928  29829  2137933  3763  5026  1622969  2468831  29634  29581  2187610  3769  4217  5083  6075  1294083  1927888  2814819  3161149  29677  29476  29428  29535
% %
% NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
% %
% NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
% %
% NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
];

data_norms = [
% double selective-norm    double diff-norm    2segment 32 bit selective-norm    2segment 64 bit selective-norm    2segment 32 bit diff-norm    2segment 64 bit diff-norm    4segment 16 bit selective-norm    4segment 32 bit selective-norm    4segment 48 bit selective-norm    4segment 64 bit selective-norm    4segment 16 bit diff-norm    4segment 32 bit diff-norm    4segment 48 bit diff-norm    4segment 64 bit diff-norm
% adaptive
55242  211088  47258  55555  123666  226601  42700  47442  51554  56472  95252  167987  234421  305504
% delaunay_n22
50797  139306  46609  54842  87641  152870  44542  48652  55121  61768  68214  109394  145944  193842
% europe_osm
162862  1415688  146698  256805  767611  1551257  113916  187318  270042  344387  623554  1145605  1870265  2438911
% hugebubbles-00020
189305  604945  143004  233430  333730  658568  112338  161658  216341  271375  257094  480404  744380  1004161
% rgg_n_2_24_s0
28637  482793  28631  28914  267488  521918  28573  29427  30875  32595  203582  382021  592894  766424
% road_usa
255174  679029  190381  292635  374236  739562  152410  208071  271419  337282  293554  538389  856935  1086841
% Stanford
28640  30017  28955  29122  29780  29961  28882  28823  28884  29410  29532  29604  29786  30597
%%
%NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN 
%%                                                                                 
%NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN 
%%                                                                                 
%NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN 
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: ELL

data_oneiteration = [
% runtime/iter double    runtime/iter 2-segment 32 bit    runtime/iter 2-segment 64 bit    runtime/iter 4-segment 16 bit    runtime/iter 4-segment 32 bit    runtime/iter 4-segment 48 bit    runtime/iter 4-segment 64 bit
% adaptive
1243878  925515  1262721  1253790  1037391  1507875  1951055
% delaunay_n22
2143975  1589403  2215294  2345843  1807897  2572111  2898286
% europe_osm
17937723  13568251  18100575  21077316  16434175  23558535  28031171
% hugebubbles-00020
4961223  4138549  4973879  6080634  5894127  9949744  11407747
% rgg_n_2_24_s0
11881290  8440543  11958179  13729766  10703325  15386340  21498475
% road_usa
6585971  4903513  6797986  6913221  5460531  7696778  11463013
% Stanford
1883049  1716009  1897415  2732005  1979349  2876010  2953384
%%
%NaN  NaN  NaN  NaN  NaN  NaN  NaN 
%%                                       
%NaN  NaN  NaN  NaN  NaN  NaN  NaN 
%%                                       
%NaN  NaN  NaN  NaN  NaN  NaN  NaN 
];

data_pagerank = [
% total iterations double    total runtime PageRank double    total iterations 2segment    switch point 2segment    total runtime 2segment    total iterations 4segment    switch point 4segment 16->32    switch point 4segment 32->48    switch point 4segment 48->64    total runtime 4segment
% adaptive
52  64681694  52  18  59591921  52  1  17  51  71070944
% delaunay_n22
26  55743368  26  12  50087051  26  1  12  25  58568529
% europe_osm
112  2009024990  112  48  1809713008  113  1  48  110  2338206437
% hugebubbles-00020
61  302634645  61  19  287535460  61  1  19  60  531522295
% rgg_n_2_24_s0
106  1259416773  106  41  1123344012  106  1  41  104  1454199242
% road_usa
31  204165104  31  14  184215088  32  1  14  30  223974705
% Stanford
118  222199877  118  51  214643438  118  1  51  115  294624394
%%
%NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN 
%%                                                         
%NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN 
%%                                                         
%NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN 
];

