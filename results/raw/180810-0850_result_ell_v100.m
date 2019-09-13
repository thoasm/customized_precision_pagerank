% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV ELL double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV ELL 2segment 32 bit  &  SpMV ELL 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV ELL 4segment 16 bit  &  SpMV ELL 4segment 32 bit  &  SpMV ELL 4segment 48 bit  &  SpMV ELL 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% Stanford
5391  1160222  31162  1427883  5354  5363  732233  1250715  30866  30981  1501883  5397  5320  5341  5310  548914  958997  1269653  1371246  30752  30809  30717  30887
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% Stanford
41178  30947  40867  40992  30821  30866  40957  40837  40945  40924  30931  30925  30793  30826
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% Stanford
36098  1427883  36512  36012  36199  1501883
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: ELL

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% Stanford
1241871  814890  1228941  614122  1031686  1341803  1249371
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% Stanford
118  146540883  118  51  126177451  118  1  51  115  144046784
];

% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV ELL double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV ELL 2segment 32 bit  &  SpMV ELL 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV ELL 4segment 16 bit  &  SpMV ELL 4segment 32 bit  &  SpMV ELL 4segment 48 bit  &  SpMV ELL 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% adaptive
136399  673118  90339  676090  68831  136453  436882  735098  61628  90570  998023  40849  69581  102588  137336  310328  474223  683870  893225  57535  66253  78249  91008
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% adaptive
59666  152418  51597  59740  90193  151815  44009  51439  55853  59583  67767  91201  124944  159306
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% adaptive
229728  676090  135335  182044  229722  998023
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: ELL

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% adaptive
901138  606682  877261  412851  627888  881567  601802
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% adaptive
52  46859224  52  18  42259753  52  1  17  51  42992246
];

% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV ELL double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV ELL 2segment 32 bit  &  SpMV ELL 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV ELL 4segment 16 bit  &  SpMV ELL 4segment 32 bit  &  SpMV ELL 4segment 48 bit  &  SpMV ELL 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% delaunay_n22
84740  1236043  70414  1492750  43222  84795  789569  1279978  50815  68742  1719321  26102  43811  63945  85300  572320  811824  1095783  1388027  46457  50673  60294  69264
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% delaunay_n22
56514  104627  51879  59871  69394  105103  48857  54893  59697  65302  50983  70173  88926  107897
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% delaunay_n22
153625  1492750  95294  125060  154285  1719321
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: ELL

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% delaunay_n22
1449880  937894  1368823  671551  936285  1239539  1068656
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% delaunay_n22
26  37696899  26  12  33002622  26  1  12  25  30918963
];

% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV ELL double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV ELL 2segment 32 bit  &  SpMV ELL 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV ELL 4segment 16 bit  &  SpMV ELL 4segment 32 bit  &  SpMV ELL 4segment 48 bit  &  SpMV ELL 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% europe_osm
1006144  10357864  494535  13055619  499610  1006587  6843114  10968124  290838  494037  15046720  290564  505638  754465  1010692  5087414  7168395  9380829  11631173  251196  288372  388910  494625
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% europe_osm
155104  950179  132024  214643  496934  945325  110432  172808  235689  298913  348778  495363  754004  1006372
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% europe_osm
1502830  13055619  796781  1145139  1506013  15046720
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: ELL

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% europe_osm
11546128  7797004  11401962  5598760  8016737  10520317  10586907
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% europe_osm
112  1293166369  112  48  1126337373  113  1  48  110  1090499359
];

Sorry, this application does not support Market Market type: [matrix array real general]
% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV ELL double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV ELL 2segment 32 bit  &  SpMV ELL 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV ELL 4segment 16 bit  &  SpMV ELL 4segment 32 bit  &  SpMV ELL 4segment 48 bit  &  SpMV ELL 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% hugebubbles-00020
420066  2670157  222097  1676564  209311  420221  2013092  3372596  136904  223294  2537253  122619  211761  315240  422203  1505694  2473766  4605869  6390009  113944  136257  179187  224254
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% hugebubbles-00020
169003  411216  125922  193742  223921  416894  98278  136256  175390  214500  149574  224014  326147  427317
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% hugebubbles-00020
644599  1676564  350450  497062  645525  2537253
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: ELL

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% hugebubbles-00020
3295036  2482991  3234433  1775981  2871147  5176885  4492397
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% hugebubbles-00020
61  200997213  61  19  187827308  61  1  19  60  276007730
];

Sorry, this application does not support Market Market type: [matrix array real general]
% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV ELL double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV ELL 2segment 32 bit  &  SpMV ELL 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV ELL 4segment 16 bit  &  SpMV ELL 4segment 32 bit  &  SpMV ELL 4segment 48 bit  &  SpMV ELL 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% rgg_n_2_24_s0
333446  6894868  180804  9267227  166149  332868  4525469  7018636  113088  180601  9980156  97317  168064  249842  334365  3407394  4642355  5994850  7379876  93138  112583  149294  181588
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% rgg_n_2_24_s0
38833  330353  38308  38790  181880  330353  38591  39272  39694  39478  120597  182624  267124  347967
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% rgg_n_2_24_s0
515790  9267227  283489  398000  516236  9980156
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: ELL

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% rgg_n_2_24_s0
7340754  4919914  7243036  3634563  4965946  6385387  6540086
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% rgg_n_2_24_s0
106  778119966  106  41  687216895  106  1  41  104  632444566
];

Sorry, this application does not support Market Market type: [matrix array real general]
% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV ELL double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV ELL 2segment 32 bit  &  SpMV ELL 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV ELL 4segment 16 bit  &  SpMV ELL 4segment 32 bit  &  SpMV ELL 4segment 48 bit  &  SpMV ELL 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% road_usa
474449  3549377  255755  4254477  236206  474640  2260983  3635738  150766  247898  5224585  138008  238964  355909  476599  1623574  2288078  3009268  3738410  128911  149652  198427  249140
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% road_usa
207895  460709  151087  228464  249714  470450  119606  158642  201423  244792  167531  249769  370104  491280
];
data_normalization_conversion = [
% 2segment 32->64 bit  &  2segment in-place conversion  &  4segment 16->32 bit  &  4segment 32->48 bit  &  4segment 48->64 bit  &  4segment in-place conversion
% road_usa
723887  4254477  391734  556943  725345  5224585
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: ELL

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% road_usa
4340439  2779458  4147521  1911057  2726186  3601350  3518064
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% road_usa
31  134553609  31  14  117178192  32  1  14  30  110818966
];

