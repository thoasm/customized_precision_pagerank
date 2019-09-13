% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01; Sparse format: CSR

runtime_pagerank_1_csr = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  1729471  2666710  15745264  38260395  61680687  86970781
% delaunay_n22
  5336093  13511178  32573370  58454403  88385958  121048917
% europe_osm
  56168139  175482239  304231410  432986884  566462488  695187555
% hugebubbles-00020
  9934102  20339493  129306408  264208627  404269367  549538780
% rgg_n_2_24_s0
  71935029  343287336  723123301  1175340639  1627475906  2097746630
% road_usa
  39423997  103244986  172367404  241492863  313273178  385020462
% stanford
  2591410  6750068  11263877  15934032  20629925  25317340
% wb-edu
  47060828  147139001  256614491  370268359  483979321  601679452
% web-BerkStan
  86314225  264333572  457931128  657148554  858802832  1059373980
% web-Google
  11924749  36134219  63196965  90866031  118955337  147856562
];
runtime_pagerank_2_csr = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  1219301  1865619  14200851  36686086  60076239  85365494
% delaunay_n22
  2407775  6089061  20684030  46539777  76494866  109165508
% europe_osm
  40737275  126590883  240433048  369174059  502659716  631400341
% hugebubbles-00020
  8641759  17484038  121308319  256168906  396234138  541528635
% rgg_n_2_24_s0
  50353817  240030180  566127335  1018294230  1470405664  1940795448
% road_usa
  24071048  62948043  120106158  189193616  260985670  332768761
% stanford
  2130155  5535320  9722521  14395161  19057885  23747470
% wb-edu
  37279145  118919399  218593923  332462903  446194629  563834207
% web-BerkStan
  88003363  271495553  467239232  667648291  867351634  1066526387
% web-Google
  8849728  28203533  52955793  80457212  109074767  137427249
];
runtime_pagerank_4_csr = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
 % adaptive
  868864  2385092  28221160  56030568  83793301  111864115
 % delaunay_n22
  1983232  7686483  30624823  61439768  94361025  131098487
 % europe_osm
  43502671  143245478  267312825  413929280  569481074  698252951
 % hugebubbles-00020
  18368854  30047835  298451177  613504407  929751429  1074949849
 % rgg_n_2_24_s0
  37250159  455540441  1030676648  1710381283  2388309352  2876784316
 % road_usa
  23808899  63618431  115843672  177213301  243732821  315487759
 % stanford
  2411117  7005533  11546283  17565967  23593840  28261301
 % wb-edu
  38527661  125053722  230175434  352079990  479586018  593129987
 % web-BerkStan
  82108433  255307339  445526609  646467563  848123389  1046533784
 % web-Google
  11145489  37619761  72237410  111468095  151954566  180281318
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01; Sparse format: ELL
runtime_pagerank_1_ell = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  1664188  2560470  15156325  36721103  59211223  83476130
% delaunay_n22
  6781502  17084241  41121378  73748733  111533089  152729085
% europe_osm
  149964272  464319140  803823219  1143325054  1495418686  1834910636
% hugebubbles-00020
  9190044  18925472  120827527  246936634  377903673  513687777
% rgg_n_2_24_s0
  40958281  195721187  412377382  670298022  928244188  1196481423
% road_usa
  72496356  189205500  315657779  442066650  573394065  704675338
% stanford
  18051406  46919193  78250922  110674352  143243056  175717607
];
runtime_pagerank_2_ell = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  1140839  1750287  13547833  35128123  57611595  81888815
% delaunay_n22
  4502031  11331249  33484268  66085635  103857983  145051815
% europe_osm
  101704020  314574348  609186875  948704443  1300782490  1640307853
% hugebubbles-00020
  8165283  16535529  114211215  240287343  371196858  506999584
% rgg_n_2_24_s0
  27665763  132012427  325967083  583936232  841878251  1110102947
% road_usa
  45479294  118619981  225298265  351727338  483060769  614340053
% stanford
  11794318  30724795  57691621  90202642  122726782  155137809
];
runtime_pagerank_4_ell = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  792389  2224995  25538682  50656261  75785444  102762689
% delaunay_n22
  3690884  12902636  45246989  85069540  128763622  175110448
% europe_osm
  99695102  331220362  608441428  928510832  1277691490  1617185380
% hugebubbles-00020
  16231957  26780846  269158151  553662241  839902603  975530404
% rgg_n_2_24_s0
  23194315  170516017  373091861  612471949  867167119  1145758509
% road_usa
  42560581  114531881  205346719  309752975  424649148  555917063
% stanford
  13872741  37972829  67917243  103528982  140434188  172857388
];
