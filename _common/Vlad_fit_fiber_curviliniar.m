function [px,fitresult_x,gof_x,py,fitresult_y,gof_y,pz,fitresult_z,gof_z] = Vlad_fit_fiber_curviliniar(a,b,c,d,s_resolution,use_weights,fitting_type,Bounding_Box)
d1=Bounding_Box(1,2)-Bounding_Box(1,1);
d2=Bounding_Box(2,2)-Bounding_Box(2,1);
d3=Bounding_Box(3,2)-Bounding_Box(3,1);

if d1>d2 && d1>d3
    [a,idx]=sort(a);
    b=b(idx);
    c=c(idx);
    d=d(idx);
elseif d2>d1 && d2>d3
    [b,idx]=sort(b);
    a=a(idx);
    c=c(idx);
    d=d(idx);
elseif d3>d1 && d3>d2
    [c,idx]=sort(c);
    a=a(idx);
    b=b(idx);
    d=d(idx); 
end

s=zeros(size(a));
for i = 2:length(a)
    s(i) = s(i-1) + sqrt( (b(i)-b(i-1))^2 + (a(i)-a(i-1))^2 + (c(i)-c(i-1))^2 );
end

% Set up fittype and options.
ft = fittype( fitting_type );
opts = fitoptions( 'Method', 'LinearLeastSquares' );

%% x direction
% Fit: 'untitled fit 1'.
[xData, yData, weights] = prepareCurveData( s, a, d );
if use_weights
    opts.Weights = weights;
end

% Fit model to data.
[fitresult_x, gof_x] = fit( xData, yData, ft, opts );

px=fitresult_x(linspace(0,s(end),s_resolution));
%px=fitresults_x(s);

%% y direction
% Fit: 'untitled fit 1'.
[xData, yData, weights] = prepareCurveData( s, b, d );
if use_weights
    opts.Weights = weights;
end

% Fit model to data.
[fitresult_y, gof_y] = fit( xData, yData, ft, opts );

py=fitresult_y(linspace(0,s(end),s_resolution));
%py=fitresults_y(s);

%% z direction
% Fit: 'untitled fit 1'.
[xData, yData, weights] = prepareCurveData( s, c, d );
if use_weights
    opts.Weights = weights;
end

% Fit model to data.
[fitresult_z, gof_z] = fit( xData, yData, ft, opts );

pz=fitresult_z(linspace(0,s(end),s_resolution));
%pz=fitresults_z(s);

end