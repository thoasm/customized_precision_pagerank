function [ data ] = plot_data_csr( save )

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
MS2 = 'markersize'; ms2 = 14;
MC = 'markerfacecolor'; mc = 'auto';

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


% ##############################################################################
%
% a plot comparing the flop count GJE-inversion vs. GJE-inversion+cond number
%
% ##############################################################################

%eval('csr_juwels48_interleaved_181127_1223'); fileSuffix = '_cpu_interleaved_8KiB';
%eval('csr_juwels48_interleaved_128B_181128_0547'); fileSuffix = '_cpu_interleaved_128B';
eval('csr_v100_181128_0251'); fileSuffix = '_gpu_v100';
%eval('csr_juwels48_split_181128_0552'); fileSuffix = '_cpu_split';


mat = linspace(1,10,10);

% plot for the matrices the relative runtime
data = data_pagerank;

h=figure(1);
set(h, 'Position', [0 0 900 600 ])
plot(mat(:),data(:,2)./data(:,2), '-', 'color', myred, LW3, lw3, MC, mc, MS, ms);
hold on;
plot(mat(:),data(:,2)./data(:,5), 'o', 'color', myblue, LW3, lw3, MC, mc, MS2, ms2);
hold on;
plot(mat(:),data(:,2)./data(:,10), 'd', 'color', mygreen, LW3, lw3, MC, mc, MS2, ms2);
hold on;
set(gca,FS,fs),grid on
ylabel('PageRank speedup',FS,fs);
%xlabel('Matrix',FS,fs);
h_legend = legend( 'IEEE Double', '2-segment CPMS', '4-segment CPMS');
set(h_legend,FS,fs,'Location','Northwest');
set(gca,FS,fs),grid on

ylim([.2 1.6]);
xlim([0 11]);
names = { 'Ada'; 'Del'; 'Eur'; 'Bub'; 'Rgg'; 'USA'; 'Std'; 'edu'; 'Brk'; 'Ggl'};
set(gca,'xtick',[1:1:10],'xticklabel',names); 

%xlim=get(gca,'xlim');
hold on
%plot(xlim,[456 456], '-', 'color', myblack, ...
%            LW3, lw3, MC, mc);
%txt1 = 'Total number of test matrices';
%text(.1,440,txt1,FS,fs)
hold off;

if( save == 1 )
    plotname = strcat('plots/rel_runtime_csr', fileSuffix);
    set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...       % figure handle
            plotname,...     % name of output file without extension
        '-painters', ...      % renderer
        '-pdf', ...           % file format
        '-r72' );             % resolution in dpi
else
    %pause()
end
hold off;



h=figure(2);
set(h, 'Position', [0 0 2000 600 ]);
% try to verify what is happening....

t = zeros(50,7);%*NaN;
o=0;
for(m = 1:10)
    
    %data_pagerank:
        %1-> iterations double; 2-> runtime double
        %3->iterations 2s; 4-> switching 2s to 64; 5-> runtime 2s
        %6->iterations 4s; 7-> switching 4s to 32; 8-> switching 4s to 48;
        %9-> switching 4s to 64; 10-> runtime 4s
    %data_oneiteration:
        %1-> runtime / iter double
        %2-> runtime / iter 2s 32 bit; 3-> runtime / iter 2s 64 bit
        %4-> runtime / iter 4s 16 bit; 5-> runtime / iter 4s 32 bit
        %6-> runtime / iter 4s 48 bit; 7-> runtime / iter 4s 64 bit
    %data_normalization_conversion:
        %1-> normalization 2s from 32 to 64; 2-> conversion 2s to double
        %3-> normalization 4s 16->32; 4-> normaliz. 4s 32->48;
        %5-> normaliz. 4s 48->64; 6->conversion 4s to double
    
    % double
    t(o+3,7) = data_pagerank(m,2);
    
    %2-segment
    t(o+2,3) = data_pagerank(m,4) * data_oneiteration(m, 2);
    t(o+2,6) = data_normalization_conversion(m,2); % precision change
    t(o+2,7) = data_pagerank(m,5) - sum(t(o+2,:));
    
    %4-segment
    t(o+1,1) = data_oneiteration(m,4)* data_pagerank(m,7);
    t(o+1,2) = data_normalization_conversion(m,3); % precision change
    t(o+1,3) = data_oneiteration(m,5)*(data_pagerank(m,8)-data_pagerank(m,7));
    t(o+1,4) = data_normalization_conversion(m,4); % precision change
    t(o+1,5) = data_oneiteration(m,6)*(data_pagerank(m,9)-data_pagerank(m,8));
    t(o+1,6) = data_normalization_conversion(m,6) - sum(data_normalization_conversion(m,3:4)); % precision change
    t(o+1,7) = data_pagerank(m,10) - sum(t(o+1,:));
    
    k = t(o+3,7);   %k ^= runtime with double precision
    t(o+1:o+3,:)= t(o+1:o+3,:)/k;

    o=o+5;
end

barh(t,'stacked'); % stacks values in each row together
h_legend = legend( '16 bit access', 'precision change', '32 bit access', 'precision change', '48 bit access', 'precision change', '64 bit access');
set(h_legend,FS,fs,'Location','Southeast');
xlabel('Normalized PageRank execution time',FS,fs);
set(gca,FS,fs),grid on
%names = { 'Adapt'; 'Delny'; 'Euro'; 'Bubble'; 'Rgg'; 'Rd\_USA';'Stanford'; 'wb\_edu'; 'BerkStan'; 'Google'};
names = { 'Ada'; 'Del'; 'Eur'; 'Bub'; 'Rgg'; 'USA'; 'Std'; 'edu'; 'Brk'; 'Ggl'};
set(gca,'ytick',[2:5:50],'yticklabel',names');


%xlim([0 2.0]);

if( save == 1 )
    plotname = strcat('plots/trace_csr', fileSuffix);
    set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...       % figure handle
            plotname,...     % name of output file without extension
        '-painters', ...      % renderer
        '-pdf', ...           % file format
        '-r72' );             % resolution in dpi
else
    %pause()
end
hold off;

data = data_basicblocks;


h=figure(3);
set(h, 'Position', [0 0 1400 600 ]);
plot(mat(:),data(:,2)./data(:,2), '-', 'color', myred, LW3, lw3, MC, mc, MS, ms);
hold on;
plot(mat(:),data(:,2)./data(:,7), 'v', 'color', myblue, LW3, lw3, MC, mc, MS2, ms2);
hold on;
plot(mat(:),data(:,2)./data(:,8), '^', 'color', myblue, LW3, lw3, MC, mc, MS2, ms2);
hold on;
plot(mat(:),data(:,2)./data(:,16), 'd', 'color', mygreen, LW3, lw3, MC, mc, MS2, ms2);
hold on;
plot(mat(:),data(:,2)./data(:,17), '*', 'color', mygreen, LW3, lw3, MC, mc, MS2, ms2);
hold on;
plot(mat(:),data(:,2)./data(:,18), 'x', 'color', mygreen, LW3, lw3, MC, mc, MS2, ms2);
hold on;
plot(mat(:),data(:,2)./data(:,19), '+', 'color', mygreen, LW3, lw3, MC, mc, MS2, ms2);
hold on;
set(gca,FS,fs),grid on
ylabel('SpMV speedup',FS,fs);
%xlabel('Matrix',FS,fs);
h_legend = legend( 'IEEE Double', '32-bit 2-seg', '64-bit 2-seg', '16-bit 4-seg', '32-bit 4-seg', '48-bit 4-seg', '64-bit 4-seg');
set(h_legend,FS,fs,'Location','westoutside');
set(gca,FS,fs),grid on

ylim([0 3]);
xlim([0 11]);
names = { 'Ada'; 'Del'; 'Eur'; 'Bub'; 'Rgg'; 'USA'; 'Std'; 'edu'; 'Brk'; 'Ggl'};
set(gca,'xtick',[1:1:10],'xticklabel',names); 

%xlim=get(gca,'xlim');
hold on
%plot(xlim,[456 456], '-', 'color', myblack, ...
%            LW3, lw3, MC, mc);
%txt1 = 'Total number of test matrices';
%text(.1,440,txt1,FS,fs)
hold off;

if( save == 1 )
    plotname = strcat('plots/rel_runtime_spmv_csr', fileSuffix);
    set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...       % figure handle
            plotname,...     % name of output file without extension
        '-painters', ...      % renderer
        '-pdf', ...           % file format
        '-r72' );             % resolution in dpi
else
    %pause()
end
hold off;


return;

% The following portion is no longer needed for the CPU plot.


eval('result_csr_v100_180813_0600');

mat = linspace(1,10,10);

% plot for the matrices the relative runtime
data = data_pagerank;

h=figure(1);
set(h, 'Position', [0 0 900 600 ])
plot(mat(:),data(:,2)./data(:,2), '-', 'color', myred, LW3, lw3, MC, mc, MS, ms);
hold on;
plot(mat(:),data(:,2)./data(:,5), 'o', 'color', myblue, LW3, lw3, MC, mc, MS2, ms2);
hold on;
plot(mat(:),data(:,2)./data(:,10), 'd', 'color', mygreen, LW3, lw3, MC, mc, MS2, ms2);
hold on;
set(gca,FS,fs),grid on
ylabel('PageRank speedup',FS,fs);
%xlabel('Matrix',FS,fs);
h_legend = legend( 'IEEE Double', '2-segment CPMS', '4-segment CPMS');
set(h_legend,FS,fs,'Location','Northwest');
set(gca,FS,fs),grid on

ylim([.6 1.6]);
xlim([0 11]);
names = { 'Ada'; 'Del'; 'Eur'; 'Bub'; 'Rgg'; 'USA'; 'Std'; 'edu'; 'Brk'; 'Ggl'};
set(gca,'xtick',[1:1:10],'xticklabel',names); 

%xlim=get(gca,'xlim');
hold on
%plot(xlim,[456 456], '-', 'color', myblack, ...
%            LW3, lw3, MC, mc);
%txt1 = 'Total number of test matrices';
%text(.1,440,txt1,FS,fs)
hold off;

if( save == 1 )
    plotname = strcat('plots/rel_runtime_csr_v100');
    set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...       % figure handle
            plotname,...     % name of output file without extension
        '-painters', ...      % renderer
        '-pdf', ...           % file format
        '-r72' );             % resolution in dpi
else
    %pause()
end
hold off;

h=figure(2);
set(h, 'Position', [0 0 2000 600 ]);
% try to verify what is happening....

t = zeros(50,7);%*NaN;
o=0;
for(m = 1:10)
    
    t(o+3,7) = data_pagerank(m,2);
    
    %2-segment
    t(o+2,3) = data_pagerank(m,4) * data_oneiteration(m, 2);
    t(o+2,6) = data_normalization_conversion(m,2); % precision change
    t(o+2,7) = data_pagerank(m,5) - sum(t(o+2,:));
    
    %4-segment
    t(o+1,1) = data_oneiteration(m,4)* data_pagerank(m,7);
    t(o+1,2) = data_normalization_conversion(m,3); % precision change
    t(o+1,3) = data_oneiteration(m,5)*(data_pagerank(m,8)-data_pagerank(m,7));
    t(o+1,4) = data_normalization_conversion(m,4); % precision change
    t(o+1,5) = data_oneiteration(m,6)*(data_pagerank(m,9)-data_pagerank(m,8));
    t(o+1,6) = data_normalization_conversion(m,6) - sum(data_normalization_conversion(m,3:4)); % precision change
    t(o+1,7) = data_pagerank(m,10) - sum(t(o+1,:));
    
    k = t(o+3,7);   %k ^= runtime with double precision
    t(o+1:o+3,:)= t(o+1:o+3,:)/k;

    o=o+5;
end

barh(t,'stacked'); % stacks values in each row together
h_legend = legend( '16 bit access', 'precision change', '32 bit access', 'precision change', '48 bit access', 'precision change', '64 bit access');
set(h_legend,FS,fs,'Location','Southeast');
xlabel('Normalized PageRank execution time',FS,fs);
set(gca,FS,fs),grid on
%names = { 'Std'; 'Ada'; 'Delny'; 'Euro'; 'Bubble'; 'Rgg'; 'Rd\_USA';'Stanford'; 'wb\_edu'; 'BerkStan'; 'Google'};
names = { 'Ada'; 'Del'; 'Eur'; 'Bub'; 'Rgg'; 'USA'; 'Std'; 'edu'; 'Brk'; 'Ggl'};
set(gca,'ytick',[2:5:50],'yticklabel',names')



if( save == 1 )
    plotname = strcat('plots/trace_csr_v100');
    set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...       % figure handle
            plotname,...     % name of output file without extension
        '-painters', ...      % renderer
        '-pdf', ...           % file format
        '-r72' );             % resolution in dpi
else
    %pause()
end
hold off;

data = data_basicblocks;


h=figure(3);
set(h, 'Position', [0 0 1400 600 ])
plot(mat(:),data(:,2)./data(:,2), '-', 'color', myred, LW3, lw3, MC, mc, MS, ms);
hold on;
plot(mat(:),data(:,2)./data(:,7), 'v', 'color', myblue, LW3, lw3, MC, mc, MS2, ms2);
hold on;
plot(mat(:),data(:,2)./data(:,8), '^', 'color', myblue, LW3, lw3, MC, mc, MS2, ms2);
hold on;
plot(mat(:),data(:,2)./data(:,16), 'd', 'color', mygreen, LW3, lw3, MC, mc, MS2, ms2);
hold on;
plot(mat(:),data(:,2)./data(:,17), '*', 'color', mygreen, LW3, lw3, MC, mc, MS2, ms2);
hold on;
plot(mat(:),data(:,2)./data(:,18), 'x', 'color', mygreen, LW3, lw3, MC, mc, MS2, ms2);
hold on;
plot(mat(:),data(:,2)./data(:,19), '+', 'color', mygreen, LW3, lw3, MC, mc, MS2, ms2);
hold on;
set(gca,FS,fs),grid on
ylabel('SpMV speedup',FS,fs);
%xlabel('Matrix',FS,fs);
h_legend = legend( 'IEEE Double', '32-bit 2-seg', '64-bit 2-seg', '16-bit 4-seg', '32-bit 4-seg', '48-bit 4-seg', '64-bit 4-seg');
set(h_legend,FS,fs,'Location','westoutside');
set(gca,FS,fs),grid on

ylim([0 3]);
xlim([0 11]);
names = { 'Ada'; 'Del'; 'Eur'; 'Bub'; 'Rgg'; 'USA'; 'Std'; 'edu'; 'Brk'; 'Ggl'};
set(gca,'xtick',[1:1:10],'xticklabel',names); 

%xlim=get(gca,'xlim');
hold on
%plot(xlim,[456 456], '-', 'color', myblack, ...
%            LW3, lw3, MC, mc);
%txt1 = 'Total number of test matrices';
%text(.1,440,txt1,FS,fs)
hold off;

if( save == 1 )
    plotname = strcat('plots/rel_runtime_spmv_csr_v100');
    set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...       % figure handle
            plotname,...     % name of output file without extension
        '-painters', ...      % renderer
        '-pdf', ...           % file format
        '-r72' );             % resolution in dpi
else
    %pause()
end
hold off;


