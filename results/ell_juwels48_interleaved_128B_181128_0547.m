% Optimizations: vectorized own_csr own_ell single_block
% all kernel runtimes averaged over 50 executions (and measured in ns); PageRank executed with damping factor of 8.500000e-01; Sparse format: CSR; Mode: interleaved with 128 Bytes; Using thread-buffer;

data_basicblocks = [
% read/write double  &  SpMV ELL double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV ELL 2segment 32 bit  &  SpMV ELL 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV ELL 4segment 16 bit  &  SpMV ELL 4segment 32 bit  &  SpMV ELL 4segment 48 bit  &  SpMV ELL 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% adaptive
263366  6310153  170961  5048239  462493  1305913  6380334  8022145  210752  253956  7773667  462450  483049  561111  1173671  6334367  8862337  9562260  9914234  292916  366493  410839  449607
% delaunay_n22
138960  18532311  114776  15943006  280266  320035  15157585  24372876  95757  117181  17768561  265909  209070  265796  681122  25937924  43441039  61348553  73774700  157451  181549  206829  195010
% europe_osm
6882022  85827136  3155925  112841072  5053001  10978078  92600058  127011033  3230835  3958444  143067996  4056532  5017783  6421802  10169216  134801057  178762022  207166306  211904411  3156393  3844269  4211920  4728832
% hugebubbles-00020
2686334  20895261  1241463  13580337  2051335  4564197  18756911  29634699  1404150  1517696  27399763  1759613  2141309  2478306  4059792  25424758  35618280  48767015  56891377  1311684  1676674  1841285  1893329
% rgg_n_2_24_s0
1672053  136031922  891455  108390468  1469876  2108524  100100674  145326265  1036550  1151503  120621293  1329618  1590360  1895720  3166008  167734676  246332241  321209662  398644000  1049106  1267742  1404235  1471147
% road_usa
3035561  30877116  1321114  38747566  2280089  5044055  30062964  34215992  1471225  1744727  50878370  1830505  2395237  2883008  4672514  34560363  39538012  39119343  39209565  1479547  1823765  1999518  2122383
% stanford
18047  9138020  70359  11387834  17597  55067  7761499  11461245  70357  71203  11739007  27338  27285  36523  61144  16767901  19050258  20968579  22440651  28458  29186  30311  29978
];
data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% adaptive
18061  477243  18684  20394  528999  755144  18443  20436  20749  20416  509607  688038  923545  1087302
% delaunay_n22
24550  217638  24182  24329  256711  332147  20874  23262  23093  23569  278984  345990  472620  515915
% europe_osm
21806  6432471  21779  20703  5926150  7932981  21827  21163  21044  21152  4155079  6115109  7790423  9403031
% hugebubbles-00020
16888  2541383  16317  16433  2419494  3204942  16522  16342  16182  16110  1798351  2663544  3380991  3860149
% rgg_n_2_24_s0
23180  1961443  22906  23927  1836506  2557794  24007  23630  24285  23597  1354325  1964791  2542640  3063321
% road_usa
18537  2834277  18695  19746  2688035  3663515  17372  19696  19894  19943  2032283  2935139  3784779  4411838
% stanford
21651  64931  21864  21705  68167  71366  22024  22061  21785  21924  28383  22761  26173  37839
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% adaptive
1339875  5048239  803307  964663  1843906  7773667
% delaunay_n22
725075  15943006  350644  437213  1008322  17768561
% europe_osm
13186299  112841072  8666050  10147732  16291921  143067996
% hugebubbles-00020
5548872  13580337  3489502  4066635  6437205  27399763
% rgg_n_2_24_s0
4068483  108390468  2805660  3309314  4933263  120621293
% road_usa
5940438  38747566  4054671  4735537  7222810  50878370
% stanford
106919  11387834  68898  82751  119236  11739007
];
% all PageRank runtimes averaged over 5 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: ELL
data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% adaptive
7349919  7490206  7225116  7138077  8152824  10760327  8274811
% delaunay_n22
18989420  16714498  18641027  25629727  44446529  62926638  40419213
% europe_osm
93472880  101872736  91714704  132884712  187790025  215146914  122765547
% hugebubbles-00020
23931714  22600693  23727669  27033046  29736272  53469795  29193668
% rgg_n_2_24_s0
138021638  102887558  133896947  160237496  175305133  326499253  226715431
% road_usa
34752605  34031904  34399321  35899176  43866572  44718115  32332422
% stanford
9045313  8103474  9066131  16465568  19066343  21014912  12917690
];
data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% adaptive
65  477744769  65  7  485366634  84  1  2  83  913676420
% delaunay_n22
64  1215322934  64  15  1197511207  88  1  22  86  5112345256
% europe_osm
118  11029799878  118  50  11558137380  119  1  51  116  24186290977
% hugebubbles-00020
77  1842741993  77  12  1855237167  82  1  2  79  4329951161
% rgg_n_2_24_s0
89  12283925814  89  28  11263912657  93  1  2  91  30139314345
% road_usa
117  4066054827  117  51  4084702704  117  1  51  115  5258643018
% stanford
118  1067346967  118  51  1040306456  118  1  51  115  2381966062
];
data_pagerank_speedups = [
% double  &  2 segment  &  4 segment
1  0.984297  0.522882
1  1.01487  0.237723
1  0.954289  0.456035
1  0.993265  0.42558
1  1.09056  0.407572
1  0.995435  0.773214
1  1.02599  0.448095
];