% Optimizations: vectorized own_csr own_ell single_block
% all PageRank runtimes averaged over 10 executions (and measured in ns); PageRank executed with damping factor of 8.500000e-01; Sparse format: CSR; Mode: interleaved with 128 Bytes; Using thread-buffer;

runtime_pagerank_4_csr = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
 % adaptive
  14083314  34620312  308914058  605859806  898361663  1139964611
 % delaunay_n22
  35437332  123235990  396118296  701892708  1013615297  1210746723
 % europe_osm
  563992543  1851509526  3286337184  4813177722  6410081080  7438228650
 % hugebubbles-00020
  144991551  196477054  1356548779  2717030746  4099803905  4860508932
];
% all PageRank runtimes averaged over 10 executions (and measured in ns), PageRank executed with damping factor of 8.500000e-01; Sparse format: ELL; Mode: interleaved with 128 Bytes; Using thread-buffer;
runtime_pagerank_4_ell = [
%  1.000000e-02  1.000000e-04  1.000000e-06  1.000000e-08  1.000000e-10  1.000000e-12
% adaptive
  16070580  39020170  365720704  717845510  1062986718  1284648302
% delaunay_n22
  141266992  515191768  1766400657  3261101440  4802941015  5321872936
% europe_osm
  92954111650  93022367585  93129521110  93286645595  93317016966  93279929957
% hugebubbles-00020
  142090831  196600173  1434136622  2882000818  4340567057  5005666065
];
