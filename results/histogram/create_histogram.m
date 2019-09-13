function [] = create_histogram(save)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if (nargin < 1)
    save = 0;
end
loaded = load('row_nnz_ours.mat', 'abbr_name', 'full_name', 'row_nnz');

% set different color scheme
mycolors=[
          0    0.4470    0.7410
     0.8400    0.3250    0.0980
     0.9290    0.6940    0.1250
     0.4940    0.1840    0.5560
     0.4660    0.6740    0.1880
     0.3010    0.7450    0.9330
     0.6350    0.0780    0.1840
     0.2400    0.2400    0.2400
];


myblue    = mycolors(1,:);
myorange  = mycolors(2,:);
myyellow  = mycolors(3,:);
mymagenta = mycolors(4,:);
mygreen   = mycolors(5,:);
mycyan    = mycolors(6,:);
myred     = mycolors(7,:);
myblack   = mycolors(8,:);


meanColor = myred;
medianColor = myblue;
maxColor = mygreen;
correct_ticks=[1,2,4];

% 9 numbers up
% 

histogramWidth = 350;
histogramHeight = 200;

sz = size(loaded.row_nnz);
output_folder='histogram';
for i = 1:sz(1)
    h = figure(i);
    clf(h);
    hold on;
    
    xPos = 100 + floor((i-1) / 3) * (histogramWidth + 10);
    yPos = 100 + mod(i-1, 3) * (histogramHeight + 100);
    
    set(h, 'Position', [xPos yPos histogramWidth histogramHeight ]);
    hist = histogram(loaded.row_nnz{i}, 'FaceAlpha',1, 'FaceColor',[0 0 0],'BinMethod','scott');
    xlabel('Non-zeros per row');
    ylabel('Count of rows');
    curYLim = ylim;
    medianCorrection = 1;
    if (i == 9)
        medianCorrection = 1.3;
    elseif(i == 4)
        curYLim(2) = curYLim(2) * 0.9;
    end
    %set(gca, 'XScale', 'log')
    set(gca, 'XScale', 'lin');
    set(gca, 'YScale', 'log');
    %line([2,2], [0 10^8], "Color", [1,0,0]);
    curMean = mean(loaded.row_nnz{i});
    curMedian = median(loaded.row_nnz{i});
    curMax = max(loaded.row_nnz{i});
    
    if (i == 4)
        alignMean = 'right';
        alignMedian = 'right';
        alignMax = 'right';
    else
        alignMean = 'left';
        alignMedian = 'right';
        alignMax = 'right';
    end
    
    line([curMean, curMean], ylim, 'LineStyle', '--', 'LineWidth', 1, 'Color', meanColor);
    text(curMean, curYLim(2), [' ' num2str(curMean) ' '], 'Color', meanColor, 'HorizontalAlignment', alignMean);
    line([curMedian, curMedian], ylim, 'LineStyle', '-.', 'LineWidth', 1, 'Color', medianColor);
    text(curMedian, curYLim(2) * medianCorrection, [' ' int2str(curMedian) ' '], 'Color', medianColor, 'HorizontalAlignment', alignMedian);
    line([curMax, curMax], ylim, 'LineStyle', ':', 'LineWidth', 1, 'Color', maxColor);
    text(curMax, curYLim(2), [' ' int2str(curMax) ' '], 'Color', maxColor, 'HorizontalAlignment', alignMax);
    
    legend('count', 'mean', 'median', 'max', 'Location', 'East');
    if (ismember(i, correct_ticks))
        %hist.BinCounts
        %hist.Values
        curXLim = xlim;
        xticks(0:1:fix(curXLim(2)));
    else
        set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto');
    end
    %curXLim = xlim;
    %curXLim = [0, fix(curXLim(2))];
    %xlim(curXLim);
    %xticks(curXLim(1):1:curXLim(2));
    hold off;
    if( save )
        plotname = strcat(output_folder, '/histogram_', loaded.abbr_name{i});
        set(gcf, 'Color', 'white'); % white bckgr
        export_fig( gcf, ...       % figure handle
                plotname,...     % name of output file without extension
                '-painters', ...      % renderer
                '-pdf', ...           % file format
                '-r72' );             % resolution in dpi
    end
end

end

