clc;clear;
addpath('common/')

height=80e-3;               % channel height in m
temperature=20;             % temperature in C
[density,cp,lambda,viscosity,Pr]=water_properties(temperature);           % in m^2/s

load('/Users/vlad/Owncloud/Research/_data_literature/Moser_DNS_mean.mat')

fig1=figure(1);
hold all

for Re_tau=[180,395,590]
    switch Re_tau
        case 180
            u_tau=Re_tau/(height/2)*viscosity;          % estimation from shear Re from Pope
            viscous_length=height/2 /Re_tau;            % in m
            viscous_time=viscosity/u_tau^2;             % in s
            load('/Users/vlad/Documents/Research/data/literature/chandata/chan180/balances/chan180_kbal.mat')
            % kolmogorov
            t_kol=((-chan180_bal.dissip/viscous_time^2).^(-1)).^0.5;
            l_kol=((-chan180_bal.dissip/viscous_length^4).^(-1)).^0.25;

            length_travel=avg_DNS_180.Umean*u_tau.*t_kol;
            plot(chan180_bal.y1*viscous_length*1e3,length_travel*1e3,'Linewidth',3)

        case 395
            u_tau=Re_tau/(height/2)*viscosity;          % estimation from shear Re from Pope
            viscous_length=height/2 /Re_tau;            % in m
            viscous_time=viscosity/u_tau^2;             % in s
            load('/Users/vlad/Documents/Research/data/literature/chandata/chan395/balances/chan395_kbal.mat')
            % kolmogorov
            t_kol=((-chan395_kbal.dissip/viscous_time^2).^(-1)).^0.5;
            l_kol=((-chan395_kbal.dissip/viscous_length^4).^(-1)).^0.25;

            length_travel=avg_DNS_395.Umean*u_tau.*t_kol;
            plot(chan395_kbal.y1*viscous_length*1e3,length_travel*1e3,'Linewidth',3)

        case 590
            u_tau=Re_tau/(height/2)*viscosity;          % estimation from shear Re from Pope
            viscous_length=height/2 /Re_tau;            % in m
            viscous_time=viscosity/u_tau^2;             % in s
            load('/Users/vlad/Documents/Research/data/literature/chandata/chan590/balances/chan590_kbal.mat')
            % kolmogorov
            t_kol=((-chan590_kbal.dissip/viscous_time^2).^(-1)).^0.5;
            l_kol=((-chan590_kbal.dissip/viscous_length^4).^(-1)).^0.25;

            length_travel=avg_DNS_590.Umean*u_tau.*t_kol;
            plot(chan590_kbal.y*viscous_length*1e3,length_travel*1e3,'Linewidth',3)

    end
box on
grid on

set(gca,'TickLabelInterpreter','latex','Fontsize',20,'Linewidth',2)
legend('$Re_{\tau}=180$','$Re_{\tau}=395$','$Re_{\tau}=590$','Interpreter','Latex','Location','best')
xlabel('$z$ in mm','Interpreter','latex')
ylabel('$l$ in mm','Interpreter','latex')
title('Travel length during one Kolmogorov-sized eddy turnover time','Interpreter','latex')
end