% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% Stanford
6023  130003  33142  56182  5992  6024  99114  176427  32784  33098  135673  6103  6034  5992  6054  87035  132375  193515  268782  32919  32820  32839  32791
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% Stanford
43299  32739  43327  43188  33092  32825  43379  43360  43989  44861  32864  33041  32848  32874
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% Stanford
39229  56182  39190  39102  39203  135673
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% Stanford
207895  171978  204363  143795  196297  256515  198621
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% Stanford
118  24531660  118  51  22731792  118  1  51  115  27369613
];

% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% adaptive
137000  659185  93556  407477  68967  136691  364183  771982  63550  92359  738234  40872  69604  102890  137316  266749  425507  668976  951526  56073  64039  76292  92539
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% adaptive
63988  154419  53637  61952  91643  154756  41452  47541  52952  56862  65758  91827  123638  157600
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% adaptive
230578  407477  135485  181666  228969  738234
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% adaptive
887998  533184  863880  366122  578211  863369  580879
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% adaptive
52  46175906  52  18  40141569  52  1  17  51  41204483
];

% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% delaunay_n22
84676  567864  68601  335003  43183  84750  325449  649675  48262  67938  560117  26052  43765  63859  85327  234091  390361  620667  914902  45757  48548  57698  66954
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% delaunay_n22
53589  104653  49187  56737  67521  105369  45877  51419  56848  62010  48925  66803  85369  107827
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% delaunay_n22
151932  335003  93732  123473  152510  560117
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% delaunay_n22
750467  453983  715575  318131  503065  762300  585044
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% delaunay_n22
26  19512148  26  12  16407811  26  1  12  25  17595793
];

% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% europe_osm
1006241  2626275  494800  2075047  499716  1006626  1931531  2859507  290406  494364  4038110  290950  505790  754635  1010754  1698937  2118018  2792129  3528582  250081  288585  388647  494411
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% europe_osm
158225  946923  130289  213449  495059  943894  109700  172230  235117  298803  347876  499970  753528  1001764
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% europe_osm
1502672  2075047  796804  1144929  1505527  4038110
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% europe_osm
3739865  2743502  3707781  2098637  2894554  3892097  3209620
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% europe_osm
112  418864971  112  48  375308392  113  1  48  110  398666670
];

Sorry, this application does not support Market Market type: [matrix array real general]
% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% hugebubbles-00020
420802  2373000  220467  1049575  209356  420156  1836455  3090958  135615  220821  1905050  122274  211799  315194  422142  1432521  2441802  4408054  5864650  115449  136896  177303  221487
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% hugebubbles-00020
165245  409274  122489  190020  221653  410460  94976  132795  171704  210966  147144  221060  325411  427324
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% hugebubbles-00020
643124  1049575  350079  494718  644301  1905050
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% hugebubbles-00020
3011070  2305506  2959463  1701412  2844471  4954302  4268061
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% hugebubbles-00020
61  183675307  61  19  172101305  61  1  19  60  265392927
];

Sorry, this application does not support Market Market type: [matrix array real general]
% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% rgg_n_2_24_s0
333550  6453977  179920  2949230  166232  332965  2591614  7286800  112317  181235  3650913  97278  168087  249679  334386  1434826  3135027  6333917  10867419  94089  111527  148341  182579
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% rgg_n_2_24_s0
37289  329593  37011  38220  182133  330184  36981  38142  38025  38149  121310  183009  266843  347945
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% rgg_n_2_24_s0
515613  2949230  284093  398509  516264  3650913
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% rgg_n_2_24_s0
6894506  3013136  6777585  1598749  3448256  6744856  6516122
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% rgg_n_2_24_s0
106  730817718  106  41  570560634  106  1  41  104  583936791
];

Sorry, this application does not support Market Market type: [matrix array real general]
% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% road_usa
475356  1373836  248835  1046447  236326  474744  869142  1444058  149179  247126  2004923  137981  239039  355944  476796  716411  955124  1260709  1623200  126710  147890  197438  248547
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% road_usa
211028  461689  148897  225796  247846  459924  117138  156269  198736  241677  164598  249613  366942  486156
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% road_usa
722651  1046447  390909  555592  723978  2004923
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% road_usa
2088732  1350279  2007212  975092  1396995  1863785  1545062
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% road_usa
31  64750702  31  14  56146948  32  1  14  30  56698265
];

% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% wb-edu
196396  4008036  120830  1325979  98512  196445  3200157  4110066  77201  117616  1770557  58084  99713  147647  197562  3060211  3412659  4283703  5404163  66367  79602  100299  119946
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% wb-edu
108328  207628  82758  110865  118711  205968  76287  91931  108094  122316  84590  121091  167184  220181
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% wb-edu
316721  1325979  179938  247599  316739  1770557
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% wb-edu
4363872  3472760  4316017  3223666  3636849  4584718  4277965
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% wb-edu
114  497481512  114  47  457509403  115  1  47  112  487098907
];

% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% web-BerkStan
7453  7169026  28115  167380  4900  7121  7273131  7942746  27965  27826  246433  4938  4913  5543  6914  6523756  6837169  8015775  9596949  27921  27877  27927  27955
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% web-BerkStan
36941  39619  36899  37121  27820  40524  36622  36630  37461  37755  27703  28111  37766  39975
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% web-BerkStan
36429  167380  32996  32933  39337  246433
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% web-BerkStan
7223558  7362978  7172847  6555856  6887524  6912369  7186332
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% web-BerkStan
119  859603467  119  50  870643184  119  1  50  51  846536178
];

% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% web-Google
20182  956921  39109  123056  5345  20182  745891  1237130  29296  38844  209393  5350  5353  7616  20149  705518  816928  1321588  2069865  29364  29251  29381  38463
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% web-Google
47157  44355  40161  48522  39118  44018  40371  40165  42452  48978  29088  38789  40751  43439
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% web-Google
54173  123056  34842  40294  54248  209393
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% web-Google
1059324  834691  1047819  769262  899599  1408243  1204259
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% web-Google
116  122881638  116  47  112543018  116  1  47  114  140020757
];

