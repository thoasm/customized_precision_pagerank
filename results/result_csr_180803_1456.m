% all kernel runtimes averaged over 100 executions (and measured in ns)

data_basicblocks = [
% read/write double  &  SpMV CSR double  &  norm double  &  inplace conversion 2segment  &  read/write 2segment 32 bit  &  read/write 2segment 64 bit  &  SpMV CSR 2segment 32 bit  &  SpMV CSR 2segment 64 bit  &  norm 2segment 32 bit  &  norm 2segment 64 bit  &  inplace conversion 4segment  &  read/write 4segment 16 bit  &  read/write 4segment 32 bit  &  read/write 4segment 48 bit  &  read/write 4segment 64 bit  &  SpMV CSR 4segment 16 bit  &  SpMV CSR 4segment 32 bit  &  SpMV CSR 4segment 48 bit  &  SpMV CSR 4segment 64 bit  &  norm 4segment 16 bit  &  norm 4segment 32 bit  &  norm 4segment 48 bit  &  norm 4segment 64 bit
% adaptive
202777  1175679  120675  609019  102024  203828  729814  1486572  79626  127296  1200483  66804  105211  156728  211979  526944  861656  1550410  1952132  70979  87649  126522  149640
% delaunay_n22
125772  926436  86339  500525  63590  126199  663976  1258231  62714  87422  892059  41778  65636  97311  131337  512655  875448  1462620  1995290  54664  63098  87294  99663
% europe_osm
1505585  4080601  736154  3108223  748023  1507981  3145048  4531776  464195  767811  7177174  483910  768945  1154065  1565113  2739295  3540831  5401125  6488101  390837  529903  862927  1158383
% hugebubbles-00020
625823  4194266  321537  1571401  312669  629234  3351953  6078153  202845  333797  3302050  202953  321709  481979  653316  2725158  4945084  9242299  12364207  169767  235767  356740  477049
% rgg_n_2_24_s0
497396  10281719  259810  4419960  247893  498411  6270681  15245171  163894  268778  5857959  161048  255111  382024  517779  4231522  8447030  15473885  29920570  139484  191742  284169  373526
% road_usa
706732  2407458  358287  1565613  352808  710524  1719588  2847193  226328  374080  3463051  228727  363054  544133  737511  1443557  2067794  3039897  3442808  190351  264799  399346  540374
% Stanford
4839  322692  30572  80636  3761  5009  225051  488260  29979  29783  162270  3770  4181  5070  6079  189750  324909  616572  1141530  29735  29775  29740  29884
% wb-edu
292958  9372804  161580  1985668  146449  293492  7910431  11265096  105325  169576  2850360  95458  150896  225348  305133  7297470  10335861  14467688  26484177  94844  120524  175328  205720
% web-BerkStan
22272  13464553  37759  249266  6111  22311  11270141  12570350  29819  37942  363453  6209  6990  16224  23861  10049752  11577320  13876293  26532796  29661  29655  37099  37899
% web-Google
29102  2143594  38473  182603  7453  29222  1651445  2529529  29938  38247  317721  7614  8638  23361  31181  1269118  2094412  3796440  6591342  31114  35861  38023  38503
];

data_norms = [
% double selective-norm  &  double diff-norm  &  2segment 32 bit selective-norm  &  2segment 64 bit selective-norm  &  2segment 32 bit diff-norm  &  2segment 64 bit diff-norm  &  4segment 16 bit selective-norm  &  4segment 32 bit selective-norm  &  4segment 48 bit selective-norm  &  4segment 64 bit selective-norm  &  4segment 16 bit diff-norm  &  4segment 32 bit diff-norm  &  4segment 48 bit diff-norm  &  4segment 64 bit diff-norm
% adaptive
38238  213216  37058  37496  125074  226890  36849  36792  36811  37454  94689  167987  228646  302543
% delaunay_n22
37273  139189  37539  37263  87522  151882  37196  36885  36951  37252  66110  108567  144873  185788
% europe_osm
37252  1416952  37013  36917  767485  1550380  36892  36855  37051  36937  621691  1145682  1875301  2433854
% hugebubbles-00020
37460  605458  37776  37705  333938  658039  37836  37735  37478  38000  257547  479199  741696  992781
% rgg_n_2_24_s0
37191  483164  36990  36903  267800  522262  37091  36823  36842  37163  202553  381518  591728  763842
% road_usa
36874  679012  37220  36894  374158  738234  36882  36862  36845  37090  292580  538036  850298  1088776
% Stanford
37539  30596  37120  37027  29829  29851  36932  36863  37016  37046  30129  29837  29975  30496
% wb-edu
36935  293674  37016  36920  168426  316965  37222  36906  36897  36935  126033  228793  337603  459899
% web-BerkStan
36781  45825  37052  37058  38053  46707  37110  36854  36981  36914  29611  38241  38703  46704
% web-Google
37315  46923  36925  37209  38304  50270  36679  37105  36908  37140  37347  41532  46420  54359
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01 and an epsilon of 1.000000e-10; Sparse format: CSR

data_oneiteration = [
% runtime/iter double  &  runtime/iter 2-segment 32 bit  &  runtime/iter 2-segment 64 bit  &  runtime/iter 4-segment 16 bit  &  runtime/iter 4-segment 32 bit  &  runtime/iter 4-segment 48 bit  &  runtime/iter 4-segment 64 bit
% adaptive
1469296  1010400  1441502  1293912  1089577  1849593  1592663
% delaunay_n22
1157039  895503  1118702  1249356  1022920  1646604  1393911
% europe_osm
5698582  4436202  5688511  6686822  5058862  7724388  7244573
% hugebubbles-00020
5063793  4232892  4999837  6210356  5556722  10220537  9813045
% rgg_n_2_24_s0
10894364  6932554  10745199  8971985  8827636  16115461  13936219
% road_usa
3431064  2576947  3376507  3662197  2852980  4195711  4422042
% Stanford
387503  291358  383538  473358  385352  674786  492920
% wb-edu
9874586  8481459  9779642  15128270  10668353  14926914  11689612
% web-BerkStan
13560330  11514022  13444187  20043098  11683721  11379129  13492527
% web-Google
2198993  1774760  2176958  2640869  2176168  3874651  3075173
];

data_pagerank = [
% total iterations double  &  total runtime PageRank double  &  total iterations 2segment  &  switch point 2segment  &  total runtime 2segment  &  total iterations 4segment  &  switch point 4segment 16->32  &  switch point 4segment 32->48  &  switch point 4segment 48->64  &  total runtime 4segment
% adaptive
52  76403433  52  18  67807442  52  1  17  51  84406629
% delaunay_n22
26  30083027  26  12  26908478  26  1  12  25  36193382
% europe_osm
112  638241221  112  48  580110767  113  1  48  110  752276426
% hugebubbles-00020
61  308891400  61  19  291989600  61  1  19  60  538388598
% rgg_n_2_24_s0
106  1154802589  106  41  987092738  106  1  41  104  1411082110
% road_usa
31  106363009  31  14  95043628  32  1  14  30  120189582
% Stanford
118  45725427  118  51  40637107  118  1  51  115  64568443
% wb-edu
114  1125702871  114  47  1055850421  115  1  47  112  1514041544
% web-BerkStan
119  1613679291  119  50  1503599438  119  1  50  51  1521780095
% web-Google
116  255083219  116  47  233806625  116  1  47  114  368814438
];
