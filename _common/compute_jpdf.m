function [counts,XX,YY]=compute_jpdf(x,y,BinNumberX,BinNumberY)
%%%%% computes an estimate of the joint probability function between 'x' and'y'
%%% inputs
% x - vector of abcissa values
% y - vector of ordinate values
% BinNumberX - the number of bins in 'x'
% BinNumberY - the number of bins in 'x'
%%% outpus
% counts - computed values of the jpdf
% XX - x meshgrid for contourf plot
% YY - y meshgrid for contourf plot

% e.g.
% [counts,XX,YY]=compute_jpdf(x,y,10,10);
% contourf(YY,XX,counts)

[~, edges_x]=histcounts(x,BinNumberX); 
centers_x = (edges_x(1:end-1)+edges_x(2:end))/2;

[~, edges_y]=histcounts(y,BinNumberY); 
centers_y = (edges_y(1:end-1)+edges_y(2:end))/2;

[YY,XX]=meshgrid(centers_y,centers_x);

% compute the probability
counts=histcounts2(x,y,edges_x,edges_y,'Normalization','probability');

end