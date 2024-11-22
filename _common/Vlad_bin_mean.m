function [bin_center_x,mean_y,ci_x,ci_y] = Vlad_bin_mean(x,y,numBins,spacing_type,alpha)
% computes the bin mean of the quantity y in the bins of x
% e.g. x - curvature, y - spinning rate

%%% INPUT:
% x - vector of x
% y - vector of y
% numBins - number of bins
% spacing_type: 'lin' - linear spacing
%               'log' - logarithmicaly spaced
% alpha - significance level, e.g. alpha=0.05 -> confidence level 95%

%%% OUTPUT:
% bin_center_x - centers of the bins in x
% mean_y - mean value of y in the bins of x
% ci_x - confidence interval for x at significance level alpha
% ci_y - confidence interval for y at significance level alpha

    % find the bins in x
    switch spacing_type
        case 'lin'
            [~, edges] = histcounts(x, numBins);
        case 'log'
            edges = 10.^linspace(log10(min(x)),log10(max(x)),numBins+1);
    end

    binIndices = discretize(x, edges);
    
    % init
    binned_y = cell(1, numBins);
    mean_y=zeros(numBins,1);
    ci_y=zeros(numBins,2);
    ci_x=zeros(numBins,2);
    
    % compute the average for each bin
    for i = 1:numBins
        binned_y{i} = y(binIndices == i);
        mean_y(i)=mean(binned_y{i}, 'omitnan');
        ci_x(i,:) = mean_confidence_interval(x(binIndices == i),alpha);
        ci_y(i,:) = mean_confidence_interval(binned_y{i},alpha);
    end
    
    % output
    bin_center_x = (edges(1:end-1) + edges(2:end)) / 2;



end