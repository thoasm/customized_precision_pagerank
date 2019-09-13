% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% adaptive
136873  788477  84117  676551  68753  136482  530283  896076  54919  82063  980643  40748  69637  102537  136868  383700  584787  855232  1188450  45223  56165  68114  83707
% delaunay_n22
85035  1260203  60218  583121  43052  84740  521232  1317830  41094  58641  790495  25760  43562  63746  84946  340489  625333  1073214  1844286  35955  44278  49204  58265
% europe_osm
1007170  3831869  481502  3143004  499239  1006335  2724157  4424151  275626  483401  5073803  287981  504992  754194  1009351  2313488  3189436  4531380  6037864  230044  273609  378304  483823
% hugebubbles-00020
420925  4783420  212154  1676333  209144  420315  4002795  6863504  125032  212634  2511886  121330  211573  314971  422010  3377794  6108967  11195576  14017764  101810  125275  169034  213754
% rgg_n_2_24_s0
333176  17744748  172350  5564961  165949  332999  12019945  25979487  103205  172373  6238055  96543  167863  249954  334087  7186557  13367516  25618452  45797803  83967  103702  138038  173564
% road_usa
475445  2206273  236719  1614946  235921  474767  1269007  2428072  139181  238188  2552949  136706  238754  355958  476230  993310  1353272  1926328  2794168  112491  138089  188111  238306
% stanford
3470  124842  24850  55203  2969  3473  92546  171931  24963  25229  108178  2841  3071  3391  3942  78682  119602  173424  271651  24923  24864  25385  25010
% wb-edu
196790  4037288  109796  1325744  98296  196383  3159402  4033561  68165  109627  1742207  57550  99575  147519  197153  3286729  3384744  4240454  5366532  55739  68507  89625  110151
% web-BerkStan
7342  7203197  25645  167391  4138  7128  7277649  7918807  25172  25245  225549  3976  4498  5472  6814  6502572  6881911  8033588  9826614  25089  25137  25091  25068
% web-Google
20375  1017983  29418  123462  5368  20326  798243  1236013  25341  27939  199305  4645  11145  12522  20149  703195  1024101  1384929  2079619  24937  26466  26503  28240
];
data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% adaptive
26566  145582  25415  25059  82384  143285  24807  24941  24781  24821  56144  84852  115889  151665
% delaunay_n22
25107  96739  24769  24756  58193  97080  24647  24561  24499  24782  45043  59279  77975  98860
% europe_osm
25322  934336  25030  24978  483914  934821  25051  24719  24887  24753  326054  482630  741008  986381
% hugebubbles-00020
25269  401631  24953  25247  213397  402670  24875  25191  24729  24746  135163  212375  315567  421246
% rgg_n_2_24_s0
13570  322100  13367  13238  172936  321772  13166  13211  13121  13177  109915  172341  256840  336303
% road_usa
25220  451198  25203  25113  238287  451241  25044  24911  24891  24891  152350  237411  357294  475540
% stanford
23797  25367  23007  23517  25230  25338  22565  23268  23237  23382  24978  25119  25323  25077
% wb-edu
87518  198084  67428  95331  110004  197867  56607  70971  86656  102784  73493  109838  158689  209307
% web-BerkStan
23164  35259  22672  23204  25257  34263  22747  23430  23097  23166  25138  25416  27341  34815
% web-Google
32677  39175  24893  33139  25526  38236  24630  28992  28706  33207  25026  27887  35026  37640
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% adaptive
220648  676551  126677  173037  220357  980643
% delaunay_n22
144958  583121  86652  115245  144794  790495
% europe_osm
1491313  3143004  781329  1133150  1493537  5073803
% hugebubbles-00020
633914  1676333  338782  485428  635951  2511886
% rgg_n_2_24_s0
507173  5564961  273126  389822  507552  6238055
% road_usa
714120  1614946  379730  545494  714889  2552949
% stanford
28074  55203  27952  27959  28095  108178
% wb-edu
307544  1325744  170602  238617  308014  1742207
% web-BerkStan
28377  167391  28023  28074  27998  225549
% web-Google
43555  123462  35612  36958  44168  199305
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR
data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% adaptive
948734  636649  932123  433906  488464  990326  728586
% delaunay_n22
1380400  609811  1346643  380472  699879  1178328  1182115
% europe_osm
4800551  3425258  4755230  2624506  3818664  5409514  4429764
% hugebubbles-00020
5249152  4406472  5178250  3605254  3747831  11602438  7102364
% rgg_n_2_24_s0
18282909  12628775  17998127  7499523  7617034  26005393  21819661
% road_usa
2677931  1615601  2646560  1118893  1651871  2351832  2102540
% stanford
178601  141349  172166  118234  167711  224508  175886
% wb-edu
4244155  3396258  4197351  3159304  3596503  4506197  4196662
% web-BerkStan
7208461  7325151  7160474  6509869  6910189  6944284  7172640
% web-Google
1029236  805336  1015723  749296  1098119  1455666  1210524
];
data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% adaptive
65  61667715  65  7  60053634  84  1  2  83  83802103
% delaunay_n22
64  88345649  64  15  76470681  88  1  22  86  94372921
% europe_osm
118  566465061  118  50  502678262  119  1  51  116  569571855
% hugebubbles-00020
77  404184734  77  12  396180789  82  1  2  79  929625336
% rgg_n_2_24_s0
89  1627178913  89  28  1470192478  93  1  2  91  2388144078
% road_usa
117  313317969  117  51  261013422  117  1  51  115  243746859
% stanford
118  21075011  118  51  18968730  118  1  51  115  23710492
% wb-edu
114  483833758  114  47  445876290  115  1  47  112  479710188
% web-BerkStan
119  857806931  119  50  867851364  119  1  50  51  846612655
% web-Google
116  119391383  116  47  108908223  116  1  47  114  152278944
];