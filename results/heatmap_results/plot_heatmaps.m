function [ data ] = plot_heatmaps( save )

%
% loads the data specified in ''filename'' into the data array
% the size of data should be outdata(1,16)
% data stored for every test matrix is:
%  1--5:   size   (m x n)     ||   nonzeros (nnz)   ||   nnz/m   ||   stored nnz
%  6--10:  iter   ||   residual-nrm2    ||   runtime    ||   SpMV-count   ||   info
% 11--12:  precond setup t || precond application t 
% 13:      overall runtime
% 14:       energy precond setup
% 15:       energy solver
% 16:       overall energy
%

if (nargin < 1)
    save = 0;
end

% plot defaults
LW1 = 'linewidth'; lw1 = 1;
LW2 = 'linewidth'; lw2 = 2;
LW3 = 'linewidth'; lw3 = 2;
FS = 'fontsize'; fs = 20;
FS2 = 'fontsize'; fs2 = 20;
MS = 'markersize'; ms = 8;
MC = 'markerfacecolor'; mc = 'auto';
set(0,'DefaultFigureColormap',feval('parula'));

font_size = 12;
% Format for printing floating point values: here, 2 decimal places
cell_label_format = '%.2f';

close all;

% set different color scheme
mycolors=[
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0.2500    0.2500    0.2500
];

myblue    = mycolors(1,:);
myorange  = mycolors(2,:);
myyellow  = mycolors(3,:);
mymagenta = mycolors(4,:);
mygreen   = mycolors(5,:);
mycyan    = mycolors(6,:);
myred     = mycolors(7,:);
myblack   = mycolors(8,:);
% example: plot(x, y, '-o', 'color', myblue);

mtx_names_csr = { 'Ada'; 'Del'; 'Eur'; 'Bub'; 'Rgg'; 'USA'; 'Std';'edu'; 'Brk'; 'Ggl'};
mtx_names_ell = { 'Ada'; 'Del'; 'Eur'; 'Bub'; 'Rgg'; 'USA'; 'Std' };
epsilon_list = {'1e-2','1e-4','1e-6','1e-8','1e-10', '1e-12'};

xLabel = 'Matrices';
yLabel = 'Stopping threshold';

%eval('heatmap_juwels48_split_181128_0552'); fileSuffix = '_cpu_split';
%eval('heatmap_juwels48_interleaved_128B_181128_0547'); fileSuffix = '_cpu_interleaved_128B';
%eval('heatmap_juwels48_interleaved_181127_1223'); fileSuffix = '_cpu_interleaved_8KiB';
%eval('heatmap_juwels48_split_181128_0552'); fileSuffix = '_cpu_split';
eval('heatmap_v100_181128_0235'); fileSuffix = '_gpu_v100';



%colorbarVisible='on';
colorbarVisible='off';

z=runtime_pagerank_1_csr./runtime_pagerank_2_csr;
dd = z.';

figure;

h = heatmap(mtx_names_csr, epsilon_list, dd);
h.XLabel = xLabel;
h.YLabel = yLabel;
caxis([0 1.365]);

h.ColorbarVisible = colorbarVisible;
h.FontSize= font_size;
h.CellLabelFormat = cell_label_format;

plotname = strcat('heatmap_csr_2seg', fileSuffix);

if( save == 1 )
    %%% %pause();
    set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...    tl   % figure handle
            plotname,...     % name of output file without extension
        '-painters', ...      % renderer
        '-pdf', ...           % file format
        '-r72' );             % resolution in dpi
else
    %pause()
end

z=runtime_pagerank_1_csr./runtime_pagerank_4_csr;
dd = z.';

figure;

h = heatmap(mtx_names_csr, epsilon_list, dd);
h.XLabel = xLabel;
h.YLabel = yLabel;
caxis([0 1.365]);

h.ColorbarVisible = colorbarVisible;
h.FontSize= font_size;
h.CellLabelFormat = cell_label_format;

plotname = strcat('heatmap_csr_4seg', fileSuffix);

if( save == 1 )
    %%% %pause();
    set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...    tl   % figure handle
            plotname,...     % name of output file without extension
        '-painters', ...      % renderer
        '-pdf', ...           % file format
        '-r72' );             % resolution in dpi
else
    %pause()
end


z=runtime_pagerank_1_ell./runtime_pagerank_2_ell;
dd = z.';

figure;

h = heatmap(mtx_names_ell, epsilon_list, dd);
h.XLabel = xLabel;
h.YLabel = yLabel;
caxis([0 1.365]);

h.ColorbarVisible = colorbarVisible;
h.FontSize= font_size;
h.CellLabelFormat = cell_label_format;

plotname = strcat('heatmap_ell_2seg', fileSuffix);

if( save == 1 )
    %%% %pause();
    set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...    tl   % figure handle
            plotname,...     % name of output file without extension
        '-painters', ...      % renderer
        '-pdf', ...           % file format
        '-r72' );             % resolution in dpi
else
    %pause()
end



z=runtime_pagerank_1_ell./runtime_pagerank_4_ell;
dd = z.';

figure;

h = heatmap(mtx_names_ell, epsilon_list, dd);
h.XLabel = xLabel;
h.YLabel = yLabel;
caxis([0 1.365]);

h.ColorbarVisible = colorbarVisible;
h.FontSize= font_size;
h.CellLabelFormat = cell_label_format;

plotname = strcat('heatmap_ell_4seg', fileSuffix);

if( save == 1 )
    %%% %pause();
    set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...    tl   % figure handle
            plotname,...     % name of output file without extension
        '-painters', ...      % renderer
        '-pdf', ...           % file format
        '-r72' );             % resolution in dpi
else
    %pause()
end

