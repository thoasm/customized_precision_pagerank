load('adaptive.mat')
mat=tril(Problem.A);sz=size(mat);adaptive_row_nnz = mat * ones(sz(2), 1);
load('delaunay_n22.mat')
mat=tril(Problem.A);sz=size(mat);delaunay_n22_row_nnz = mat * ones(sz(2), 1);
load('europe_osm.mat')
mat=tril(Problem.A);sz=size(mat);europe_osm_row_nnz = mat * ones(sz(2), 1);
load('hugebubbles-00020.mat')
mat=tril(Problem.A);sz=size(mat);hugebubbles_row_nnz = mat * ones(sz(2), 1);
load('rgg_n_2_24_s0.mat')
mat=tril(Problem.A);sz=size(mat);rgg_n_2_24_s0_row_nnz = mat * ones(sz(2), 1);
load('road_usa.mat')
mat=tril(Problem.A);sz=size(mat);road_usa_row_nnz = mat * ones(sz(2), 1);
load('Stanford.mat')
mat=Problem.A;sz=size(mat);Stanford_row_nnz = mat * ones(sz(2), 1);
load('wb-edu.mat')
mat=Problem.A;sz=size(mat);edu_row_nnz = mat * ones(sz(2), 1);
load('web-BerkStan.mat')
mat=Problem.A;sz=size(mat);BerkStan_row_nnz = mat * ones(sz(2), 1);
load('web-Google.mat')
mat=Problem.A;sz=size(mat);Google_row_nnz = mat * ones(sz(2), 1);

row_nnz = {adaptive_row_nnz; delaunay_n22_row_nnz; europe_osm_row_nnz; hugebubbles_row_nnz; rgg_n_2_24_s0_row_nnz; road_usa_row_nnz; Stanford_row_nnz; edu_row_nnz; BerkStan_row_nnz; Google_row_nnz};
full_name = {"adaptive"; "delaunay_n22"; "europe_osm"; "hugebubbles-00020"; "rgg_n_2_24_s0"; "road_usa"; "Stanford"; "wb-edu"; "web-BerkStan"; "web-Google"};
abbr_name = {"Ada"; "Del"; "Eur"; "Bub"; "Rgg"; "USA"; "Std"; "edu"; "Brk"; "Ggl"};

save("row_nnz_ours.mat", "abbr_name", "full_name", "row_nnz");
%for i=1:size(row_nnz); nzr=0; for j=1:size(row_nnz{i}); if (row_nnz{i}(j) == 0); nzr = nzr + 1; end; end; fprintf("%s: %d\n", abbr_name{i}, nzr);end;
