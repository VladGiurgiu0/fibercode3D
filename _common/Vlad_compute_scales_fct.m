function [Re]=Vlad_compute_scales_fct(Reynolds,temperature)
%%%% input
% Reynolds - shear Reynolds number [-]
% temperature - water temperature [deg C]

%%% output
% structure Re contains: 
% y - wall-normal location [mm]
% t_kol - interpolated Kolmogorov time-scale [s]
% l_kol - interpolated Komogorov length-scale [m]

plot_check=0; % plot to check if the function is working correctly
%Re.nr = 736;
%temperature=15.7;

% add path to common functions
% addpath('/Users/vlad/Owncloud/Research/_codes/___common/')
% load DNS data
load('DNS.mat')

height=81.4e-3;                                   % channel height [m]
[~,~,~,viscosity,~]=water_properties(temperature);           % [m^2/s]
Re.nr = Reynolds;
Re.u_tau = Re.nr /(height/2)*viscosity;       % [m/s]
Re.viscous_length = height/2 /Re.nr;           % [m]
Re.viscous_time = viscosity/Re.u_tau^2;         % [s] 

% compute data from the DNS
Re180.u_tau =  180 /(height/2)*viscosity;       % [m/s]
Re180.viscous_length = height/2 /180;           % [m]
Re180.viscous_time = viscosity/Re180.u_tau^2;   % [s] 
Re180.t_kol=((-chan180_kbal.dissip/Re180.viscous_time^2).^(-1)).^0.5;      % [s]
Re180.l_kol=((-chan180_kbal.dissip/Re180.viscous_length^4).^(-1)).^0.25;   % [m]
Re180.y_plus = chan180_kbal.y1;                 % [-]
Re180.u_plus = chan180.Umean;                   % [-]
Re180.travel_length = chan180.Umean*Re180.u_tau.*Re180.t_kol; % travel length during one Kolmogorov time [m]
Re180.u =chan180.Umean*Re180.u_tau;             % [m/s]

Re395.u_tau =  395 /(height/2)*viscosity;       % [m/s]
Re395.viscous_length = height/2 /395;           % [m]
Re395.viscous_time = viscosity/Re395.u_tau^2;   % [s]
Re395.t_kol=((-chan395_kbal.dissip/Re395.viscous_time^2).^(-1)).^0.5;     % [s]
Re395.l_kol=((-chan395_kbal.dissip/Re395.viscous_length^4).^(-1)).^0.25;  % [m]
Re395.y_plus = chan395_kbal.y1;                 % [-]
Re395.u_plus = chan395.Umean;                   % [-]
Re395.travel_length = chan395.Umean*Re395.u_tau.*Re395.t_kol; % travel length during one Kolmogorov time [m]
Re395.u =chan395.Umean*Re395.u_tau;             % [m/s]

Re590.u_tau =  590 /(height/2)*viscosity;       % [m/s]
Re590.viscous_length = height/2 /590;           % [m]
Re590.viscous_time = viscosity/Re590.u_tau^2;   % [s] 
Re590.t_kol=((-chan590_kbal.dissip/Re590.viscous_time^2).^(-1)).^0.5;      % [s]
Re590.l_kol=((-chan590_kbal.dissip/Re590.viscous_length^4).^(-1)).^0.25;   % [m]
Re590.y_plus = chan590_kbal.y;                  % [-]
Re590.u_plus = chan590.Umean;                   % [-]
Re590.travel_length = chan590.Umean*Re590.u_tau.*Re590.t_kol; % travel length during one Kolmogorov time [m]
Re590.u =chan590.Umean*Re590.u_tau;             % [m/s]

Re950.u_tau =  950 /(height/2)*viscosity;       % [m/s]
Re950.viscous_length = height/2 /950;           % [m]
Re950.viscous_time = viscosity/Re950.u_tau^2;   % [s] 
Re950.t_kol=((-chan950_kbal.dissip/Re950.viscous_time^2).^(-1)).^0.5;      % [s]
Re950.l_kol=((-chan950_kbal.dissip/Re950.viscous_length^4).^(-1)).^0.25;   % [m]
Re950.y_plus = chan950_kbal.y;                  % [-]
Re950.u_plus = chan950.U;                       % [-]
Re950.travel_length = chan950.U*Re950.u_tau.*Re950.t_kol; % travel length during one Kolmogorov time [m]
Re950.u =chan950.U*Re950.u_tau;             % [m/s]

clearvars chan950_kbal chan950 chan590_kbal chan590 chan395_kbal chan395 chan180_kbal chan180

%% interpolate grid of all Re nrs on the grid of the highest Re nr.
Re180.interp.t_kol = interp1(Re180.y_plus*Re180.viscous_length,Re180.t_kol,Re950.y_plus*Re950.viscous_length);
Re180.interp.l_kol = interp1(Re180.y_plus*Re180.viscous_length,Re180.l_kol,Re950.y_plus*Re950.viscous_length);
Re180.interp.u =interp1(Re180.y_plus*Re180.viscous_length,Re180.u,Re950.y_plus*Re950.viscous_length);

Re395.interp.t_kol = interp1(Re395.y_plus*Re395.viscous_length,Re395.t_kol,Re950.y_plus*Re950.viscous_length);
Re395.interp.l_kol = interp1(Re395.y_plus*Re395.viscous_length,Re395.l_kol,Re950.y_plus*Re950.viscous_length);
Re395.interp.u =interp1(Re395.y_plus*Re395.viscous_length,Re395.u,Re950.y_plus*Re950.viscous_length);

Re590.interp.t_kol = interp1(Re590.y_plus*Re590.viscous_length,Re590.t_kol,Re950.y_plus*Re950.viscous_length);
Re590.interp.l_kol = interp1(Re590.y_plus*Re590.viscous_length,Re590.l_kol,Re950.y_plus*Re950.viscous_length);
Re590.interp.u =interp1(Re590.y_plus*Re590.viscous_length,Re590.u,Re950.y_plus*Re950.viscous_length);

%% fit
% fit t_kol and l_kol for the desired Re nr.
Re_for_fit = [180 395 590 950]';

% fitting options
ft = fittype( 'exp2' );
opts_t_kol = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_t_kol.Display = 'Off';
opts_t_kol.StartPoint = [0.662678485794532 -0.0115614694365102 0.038491266480233 -0.00263084151976273];

opts_l_kol = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_l_kol.Display = 'Off';
opts_l_kol.StartPoint = [0.00080704612114824 -0.00838017315656764 0.000237299355569943 -0.00090864875416824];

opts_u = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_u.Display = 'Off';
opts_u.StartPoint = [0.247310334293544 0.00104636552119348 -0.263908098606303 -0.00136281763731484];

for i=1:size(Re950.y_plus,1)
    t_kol_to_fit(1,i) = Re180.interp.t_kol(i);
    t_kol_to_fit(2,i) = Re395.interp.t_kol(i);
    t_kol_to_fit(3,i) = Re590.interp.t_kol(i);
    t_kol_to_fit(4,i) = Re950.t_kol(i);

    l_kol_to_fit(1,i) = Re180.interp.l_kol(i);
    l_kol_to_fit(2,i) = Re395.interp.l_kol(i);
    l_kol_to_fit(3,i) = Re590.interp.l_kol(i);
    l_kol_to_fit(4,i) = Re950.l_kol(i);

    u_to_fit(1,i) = Re180.interp.u(i);
    u_to_fit(2,i) = Re395.interp.u(i);
    u_to_fit(3,i) = Re590.interp.u(i);
    u_to_fit(4,i) = Re950.u(i);

    % fit t_kol
    [xData, yData] = prepareCurveData( Re_for_fit, t_kol_to_fit(:,i) );
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts_t_kol );

    Re.t_kol(i,1) = feval(fitresult,Re.nr);

    % fit l_kol
    [xData, yData] = prepareCurveData( Re_for_fit, l_kol_to_fit(:,i) );
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts_l_kol );

    Re.l_kol(i,1) = feval(fitresult,Re.nr);

    % fit u
    [xData, yData] = prepareCurveData( Re_for_fit, u_to_fit(:,i) );
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts_u );

    Re.u(i,1) = feval(fitresult,Re.nr);

end


%% mirror the data
Re180.interp.y = Re950.y_plus*Re950.viscous_length;
Re180.interp.y = [Re180.interp.y; abs(-Re180.interp.y(end)+flip(Re180.interp.y(1:end-1))) + Re180.interp.y(end)];
Re180.interp.l_kol= [Re180.interp.l_kol; flip(Re180.interp.l_kol(1:end-1))];
Re180.interp.t_kol= [Re180.interp.t_kol; flip(Re180.interp.t_kol(1:end-1))];

Re395.interp.y = Re950.y_plus*Re950.viscous_length;
Re395.interp.y = [Re395.interp.y; abs(-Re395.interp.y(end)+flip(Re395.interp.y(1:end-1))) + Re395.interp.y(end)];
Re395.interp.l_kol= [Re395.interp.l_kol; flip(Re395.interp.l_kol(1:end-1))];
Re395.interp.t_kol= [Re395.interp.t_kol; flip(Re395.interp.t_kol(1:end-1))];

Re590.interp.y = Re950.y_plus*Re950.viscous_length;
Re590.interp.y = [Re590.interp.y; abs(-Re590.interp.y(end)+flip(Re590.interp.y(1:end-1))) + Re590.interp.y(end)];
Re590.interp.l_kol= [Re590.interp.l_kol; flip(Re590.interp.l_kol(1:end-1))];
Re590.interp.t_kol= [Re590.interp.t_kol; flip(Re590.interp.t_kol(1:end-1))];

Re950.y = Re950.y_plus*Re950.viscous_length;
Re950.y = [Re950.y; abs(-Re950.y(end)+flip(Re950.y(1:end-1))) + Re950.y(end)];
Re950.l_kol= [Re950.l_kol; flip(Re950.l_kol(1:end-1))];
Re950.t_kol= [Re950.t_kol; flip(Re950.t_kol(1:end-1))];

Re.y = Re950.y * 1e3; % [mm]
Re.l_kol = [Re.l_kol; flip(Re.l_kol(1:end-1))];
Re.t_kol = [Re.t_kol; flip(Re.t_kol(1:end-1))];
Re.u = [Re.u; flip(Re.u(1:end-1))];

%% plot to check
if plot_check==1
figure(1);
subplot(1,2,1); hold on;
plot(Re180.interp.y*1e3,Re180.interp.t_kol,'Linewidth',2)
plot(Re395.interp.y*1e3,Re395.interp.t_kol,'Linewidth',2)
plot(Re590.interp.y*1e3,Re590.interp.t_kol,'Linewidth',2)
plot(Re950.y*1e3,Re950.t_kol,'Linewidth',2)
plot(Re.y,Re.t_kol,'Linewidth',2)
box on; grid minor;
xlabel('$y$ [mm]')
ylabel('$\tau_{\eta}$ [s]')
set(gca,'FontSize',15,'YScale','log')

legend('$Re_{\tau}=180$ DNS Moser et al. 1999', ...
    '$Re_{\tau}=395$ DNS Moser et al. 1999', ...
    '$Re_{\tau}=590$ DNS Moser et al. 1999', ...
    '$Re_{\tau}=950$ DNS Del Alamo et al. 2004', ...
    strcat('$Re_{\tau}=',num2str(Re.nr),'$ interp.'),...
    'location','northoutside')

subplot(1,2,2); hold on;
plot(Re180.interp.y*1e3,Re180.interp.l_kol*1e3,'Linewidth',2)
plot(Re395.interp.y*1e3,Re395.interp.l_kol*1e3,'Linewidth',2)
plot(Re590.interp.y*1e3,Re590.interp.l_kol*1e3,'Linewidth',2)
plot(Re950.y*1e3,Re950.l_kol*1e3,'Linewidth',2)
plot(Re.y,Re.l_kol*1e3,'Linewidth',2)
box on; grid minor;
xlabel('$y$ [mm]')
ylabel('$\eta$ [mm]')
set(gca,'FontSize',15,'YScale','log')

end
end