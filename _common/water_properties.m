function [density,cp,lambda,viscosity,Pr]=water_properties(temperature)
% water properties at 1bar
    % temperature in °C
    
    % VDI_temperature in °C
    VDI_temperature=[10:1:20 22 24 25 26 28 30 32 34 36 38 40 42 44 46 48 50 55 60 65 70 75 80 85 90 95]';
    
    % VDI density in kg/m^3
    VDI_density=[999.70 999.61 999.5 999.38 999.25 999.1 998.94 998.78 998.6 998.41 998.21 997.77...
        997.3 997.05 996.79 996.24 995.65 995.03 994.38 993.69 992.97 992.22 991.44 990.64 989.8...
        988.94 988.05 985.71 983.21 980.57 977.78 974.86 971.80 968.62 965.32 961.89]';
    
    % VDI specific heat at constant pressure in kJ/kgK
    VDI_cp=[4.195 4.194 4.193 4.191 4.19 4.189 4.188 4.187 4.186 4.186 4.185 4.183 4.182 4.182...
        4.181 4.181 4.180 4.180 4.179 4.179 4.179 4.179 4.179 4.179 4.179 4.179 4.18 4.181 4.183...
        4.185 4.188 4.192 4.196 4.2 4.205 4.211]'*1e3;
    
    % VDI heat diffusion coefficient in W/mK
    VDI_lambda=[578.78 570.85 582.89 584.89 586.86 588.8 590.70 592.57 594.42 596.23 598.01 601.49...
        604.87 606.52 608.14 611.31 614.39 617.38 620.29 623.1 625.84 628.49 631.07 633.57...
        636.00 638.35 640.64 646.04 651.02 655.59 659.78 663.58 667.01 670.08 672.8 675.17]'*1e-3;
    
    % VDI kinematic viscosity in m^2/s
    VDI_viscosity=[1.306 1.27 1.235 1.201 1.169 1.139 1.109 1.081 1.054 1.028 1.003 0.9565 0.9131...
        0.8927 0.8729 0.8355 0.8007 0.7682 0.7379 0.7095 0.6828 0.6578 0.6343 0.6122 0.5914...
        0.5717 0.5531 0.5109 0.4740 0.4415 0.4127 0.3872 0.3643 0.3439 0.3255 0.3089]'*1e-6;
    
    % VDI Prandtl number
    VDI_Pr=[9.466 9.164 8.876 8.603 8.342 8.093 7.856 7.63 7.414 7.207 7.009 6.638 6.297 6.137 5.983 5.692...
        5.424 5.175 4.943 4.728 4.527 4.34 4.164 4 3.846 3.702 3.566 3.259 2.994 2.764 2.562...
        2.384 2.227 2.088 1.964 1.853]';
    
    density=interp1(VDI_temperature,VDI_density,temperature);
    cp=interp1(VDI_temperature,VDI_cp,temperature);
    lambda=interp1(VDI_temperature,VDI_lambda,temperature);
    viscosity=interp1(VDI_temperature,VDI_viscosity,temperature);
    Pr=interp1(VDI_temperature,VDI_Pr,temperature);
end