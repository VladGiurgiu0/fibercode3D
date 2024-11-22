function [u_tau,viscous_length,viscous_time]=compute_wall_scales(Re_tau,temperature)
        height=80e-3;               % channel height in m
        [~,~,~,viscosity,~]=water_properties(temperature);           % in m^2/s
        u_tau=Re_tau/(height/2)*viscosity;          % estimation from shear Re from Pope
        viscous_length=height/2 /Re_tau;            % in m
        viscous_time=viscosity/u_tau^2;             % in s
end