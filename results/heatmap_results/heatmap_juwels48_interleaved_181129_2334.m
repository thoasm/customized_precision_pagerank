% Optimizations: vectorized own_csr own_ell single_block
% all PageRank runtimes averaged over 10 executions (and measured in ns); PageRank executed with damping factor of 8.500000e-01; Sparse format: CSR; Mode: interleaved with 8192 Bytes; Using thread-buffer;

runtime_pagerank_4_csr = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
 % adaptive
  9349244  23541389  271187387  541071084  804995044  1045299301
 % delaunay_n22
  22311043  81276761  293029630  553038108  822673281  1016918431
 % europe_osm
  466519057  1559997011  2832025069  4261246398  5747753071  6780566904
 % hugebubbles-00020
  127699051  177747205  1314279983  2641939082  3997356042  4745453890
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01; Sparse format: ELL; Mode: interleaved with 8192 Bytes; Using thread-buffer;
runtime_pagerank_4_ell = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  10243126  25759946  325810947  650861372  960526917  1183883958
% delaunay_n22
  131564130  493599895  1737741858  3208725036  4747124389  5246303721
% europe_osm
  93412917204  93339178481  93363364389  93358769033  93323139445  93338426058
% hugebubbles-00020
  127165074  179971744  1417374118  2863177889  4323324330  4975081753
];
