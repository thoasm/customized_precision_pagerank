function [ data ] = plot_runtime( save )

if (nargin < 1)
    save = false;
end

close all;

% plot defaults
LW1 = 'linewidth'; lw1 = 1;
LW2 = 'linewidth'; lw2 = 2;
LW3 = 'linewidth'; lw3 = 2;
FS = 'fontsize'; fs = 20;
FS2 = 'fontsize'; fs2 = 20;
MS = 'markersize'; ms = 8;
MS2 = 'markersize'; ms2 = 14;
MC = 'markerfacecolor'; mc = 'auto';

%file_prefix = 'runtime_plots/';

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

mtx_names_csr = { 'Ada'; 'Del'; 'Eur'; 'Bub'; 'Rgg'; 'USA'; 'Std'; 'edu'; 'Brk'; 'Ggl'};
epsilon_list = {'1e-2', '1e-4', '1e-6', '1e-8', '1e-10', '1e-12'};
%xAxis = {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12};
x_lin = linspace(1, 6, 6);

xLabel = 'Relative stopping criteria';
yLabel = 'runtime [ms]';

runtime_factor = 1e-6;

%eval('heatmap_juwels48_split_181128_0552'); fileSuffix = '_cpu_split';
%eval('heatmap_juwels48_interleaved_128B_181128_0547'); fileSuffix = '_cpu_interleaved_128B';
%eval('heatmap_juwels48_interleaved_181127_1223'); fileSuffix = '_cpu_interleaved_8KiB';
%eval('heatmap_juwels48_split_181128_0552'); fileSuffix = '_cpu_split';
eval('heatmap_v100_181128_0235'); fileSuffix = '_gpu_v100';

size_plot = [600 500];
start_pos = [20 50];

for matrix_index = 1:10

    fig = figure(matrix_index);
    set(fig, 'Position', [start_pos(1) start_pos(2) start_pos(1)+size_plot(1) start_pos(2)+size_plot(2)])
    h = plot(x_lin, runtime_factor * runtime_pagerank_1_csr(matrix_index, :), '*', 'color', myred, LW3, lw3, MC, mc, MS, ms);
    hold on;
    plot(x_lin, runtime_factor * runtime_pagerank_2_csr(matrix_index, :), 'o', 'color', myblue, LW3, lw3, MC, mc, MS, ms);
    hold on;
    plot(x_lin, runtime_factor * runtime_pagerank_4_csr(matrix_index, :), 'd', 'color', mygreen, LW3, lw3, MC, mc, MS, ms);
    
    xlabel(xLabel);
    ylabel(yLabel);
    title_string = strcat({'Runtime for matrix '''}, mtx_names_csr(matrix_index), '''');
    title(title_string);

    set(gca,FS,fs);grid on;
    leg = legend('IEEE Double', '2-segment CPMS', '4-segment CPMS');
    set(leg, 'Location', 'northwest');
    set(gca, 'xtick', [1:1:6], 'xticklabel', epsilon_list); 
    plotname = string(strcat('runtime_plot_csr', fileSuffix, '_', mtx_names_csr(matrix_index)));
    hold off;
    
    if( save )
        %%% %pause();
        set(gcf, 'Color', 'white'); % white bckgr
        export_fig( gcf, ...    tl   % figure handlCe
                plotname,...     % name of output file without extension
            '-painters', ...      % renderer
            '-pdf', ...           % file format
            '-r72' );             % resolution in dpi
    else
        %pause()
    end
end