function [centers,count]=compute_pdf(data,method,value,normalization)
%%%%%% computes an estimate of the pdf of data
%%% inputs
% data - a vector containing all the data dimension: 1xN
% method - 'BinWidth' - bin data based on the specified bin width 
%        - 'BinNumber' - bin data based on the specified number of bins
% value - if method is 'BinWidth', then this is the width of the bin
%       - if method is 'BinNumber' then this is the number of bins
% normalization - 'probability', 'pdf', 'cdf' etc. see manual of histcounts
%%% outputs
% centers - the centers of the bins
% count - the value for each bin

% e.g.
% [centers,count]=compute_pdf(data,'BinWidth',0.1,'pdf')
% plot(centers,count)

data(isnan(data))=[];

switch method
    case 'BinWidth'
        [count,edges]=histcounts(data,'BinWidth',value,'Normalization',normalization);
    case 'BinNumber'
        [count,edges]=histcounts(data,value,'Normalization',normalization);
end

centers=(edges(1:end-1)+edges(2:end))/2;
end

