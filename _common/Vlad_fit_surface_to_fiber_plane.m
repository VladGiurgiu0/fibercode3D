function [fitresult, gof] = Vlad_fit_surface_to_fiber_plane(a, b, c, d)
%% Fit: 'untitled fit 1'.
[xData, yData, zData, weights] = prepareSurfaceData( a, b, c, d );

% Set up fittype and options.
ft = fittype( 'poly11' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Weights = weights;

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, [xData, yData], zData);
% legend( h, 'untitled fit 1', 'c vs. a, b with d', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'a', 'Interpreter', 'none' );
% ylabel( 'b', 'Interpreter', 'none' );
% zlabel( 'c', 'Interpreter', 'none' );
% grid on
% view( 163.0, 66.2 );


