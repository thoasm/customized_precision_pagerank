% Optimizations: vectorized own_csr own_ell
% all PageRank runtimes averaged over 10 executions (and measured in ns); PageRank executed with damping factor of 8.500000e-01; Sparse format: CSR; Mode: Using completely separate approach;

runtime_pagerank_1_csr = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  14937216  22947445  134669296  325854839  525047792  739962211
% delaunay_n22
  28910148  72258546  174375003  311952073  472625512  647616095
% europe_osm
  466639208  1461485674  2536155044  3606110824  4715487399  5785104165
% hugebubbles-00020
  51680886  106864860  687557705  1405898550  2154324516  2933717286
% rgg_n_2_24_s0
  142025740  686924893  1443707012  2350338054  3256612729  4192801308
% road_usa
  297014474  783299514  1304674576  1828065983  2371935725  2909692585
% stanford
  6100825  15998597  26582084  37722143  48755168  59944802
% wb-edu
  120228631  385748758  672896964  969219363  1265767308  1573077803
% web-BerkStan
  14517324  45046379  77964317  112194706  146509890  180895652
% web-Google
  13207930  41883027  74234814  106629207  140265411  173726466
];
runtime_pagerank_2_csr = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  10380657  15489511  113744865  303907017  502541454  717521963
% delaunay_n22
  21705727  54051268  146198578  284269314  443389088  616188302
% europe_osm
  392290215  1231319255  2227204731  3310398216  4392158630  5469649082
% hugebubbles-00020
  43011823  88119669  621892428  1343310580  2095758624  2857768115
% rgg_n_2_24_s0
  95615339  456928613  1099716814  1998655004  2899150223  3838028423
% road_usa
  207515981  552238889  1010282418  1536546514  2087139808  2635654610
% stanford
  4929263  12920415  23029253  34073100  45304978  56332860
% wb-edu
  81395703  262158661  508492471  812590330  1115060323  1430587117
% web-BerkStan
  8939480  27514671  54101536  88086437  122299149  156252831
% web-Google
  8992297  28666489  58491697  90861637  124084378  157184427
];
runtime_pagerank_4_csr = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
 % adaptive
  7364728  20965652  250585182  499937963  746715924  983188220
 % delaunay_n22
  18958083  68451001  253236809  487816068  735323251  930939701
 % europe_osm
  438974057  1469076364  2704424678  4112296942  5576931329  6651300808
 % hugebubbles-00020
  84901456  130996745  1175691702  2399719511  3652303187  4447608395
 % rgg_n_2_24_s0
  86120623  586504200  1276563425  2082112520  2909994105  3881685156
 % road_usa
  205838063  565773613  1020723473  1541458688  2093201708  2637328765
 % stanford
  5295293  14057225  24746385  39321167  53208460  63419563
 % wb-edu
  79566079  272557125  509084918  781848312  1072120551  1371024551
 % web-BerkStan
  10313615  32444609  61067781  94962078  128455146  162073776
 % web-Google
  10745893  37028912  75178725  118920940  163931423  197022383
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01; Sparse format: ELL; Mode: Using completely separate approach;
runtime_pagerank_1_ell = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  13748719  21151574  123850528  300516690  484025272  683296803
% delaunay_n22
  74974611  188388781  451200269  808636062  1218150890  1672119243
% europe_osm
  1133198197  3486800997  6044384851  8597003622  11244882962  13797437340
% hugebubbles-00020
  44843267  92753422  596656683  1220235909  1867018554  2541127887
% rgg_n_2_24_s0
  555162002  2616501606  5549432932  8978628618  12463993145  16130803467
% road_usa
  522011825  1365419895  2277568860  3195010855  4139558145  5091150465
% stanford
  137534315  360896010  592848920  840915102  1070456743  1312152016
];
runtime_pagerank_2_ell = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  9403218  14013673  105444934  282937075  468858949  667963941
% delaunay_n22
  62005535  152825458  397501783  755593307  1170464825  1624318808
% europe_osm
  984228541  3054692588  5416944622  7994031960  10625615297  13186276067
% hugebubbles-00020
  33122738  67613095  516264610  1142616335  1791259522  2470139953
% rgg_n_2_24_s0
  398662972  1909236642  4389013264  7818423198  11193117106  14757057509
% road_usa
  348621042  910626297  1669558765  2578272953  3518002936  4461392931
% stanford
  107633789  278865905  488556833  736122376  970313430  1209410241
];
runtime_pagerank_4_ell = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  7852550  21870894  252158456  507795427  758727958  968436875
% delaunay_n22
  92266449  388611684  1455096399  2835272041  4225891572  4734535571
% europe_osm
  1786495546  5982872125  10885588778  16491467270  22135313228  24663389101
% hugebubbles-00020
  78071016  122788989  1148080895  2353113637  3573709965  4256106585
% rgg_n_2_24_s0
  761933509  6054903112  13314967068  21895625539  30258613600  33971229233
% road_usa
  702703968  1989856605  3579174443  5358541193  7204224005  8143005295
% stanford
  258908159  677134898  1156441839  1687976644  2188553495  2424228741
];
runtime_pagerank_4_csr = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
 % adaptive
  7597385  21189869  262431382  525443556  785862471  1028572710
 % delaunay_n22
  19522177  75332450  278712135  532985664  797300479  992152688
 % europe_osm
  516091143  1753347500  3226813541  4901940297  6622111941  7685800743
 % hugebubbles-00020
  86840656  132262313  1187599946  2430357186  3690411671  4463750023
 % rgg_n_2_24_s0
  84895474  619476522  1353621049  2221828563  3108813721  4093182770
 % road_usa
  204880753  576901863  1037612471  1564793271  2128935294  2661853355
 % stanford
  7607842  20226099  35164280  51876599  68921643  82902743
 % wb-edu
  81833225  281825993  528213702  807157713  1104460961  1397573746
 % web-BerkStan
  12042934  39079041  71420270  105303386  139365324  173313106
 % web-Google
  12244648  41449474  80060690  129911898  182545948  216303050
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01; Sparse format: ELL; Mode: Using completely separate approach;
runtime_pagerank_4_ell = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  8456429  22855763  279978544  558328220  833497026  1053801481
% delaunay_n22
  105571781  410941316  1496551307  2851519406  4255691162  4739542435
% europe_osm
  2144953965  7300435986  13051066315  19278838728  25514292046  27959330986
% hugebubbles-00020
  104358731  154531709  1328604607  2695120343  4083391062  4755578543
% rgg_n_2_24_s0
  584078253  4637523981  10187399195  16796280297  23250490465  26987577723
% road_usa
  725940340  2063549251  3823139298  5890625543  8005240842  8928001553
% stanford
  196324683  537621241  917443225  1328206721  1724737999  1960391426
];