clc; clear;

%% input
% water
desired_nr_fibers = 200; % [-]
measurement_volume = 40 * 24 * 18 * 1e-9; % [m^3]
volume_water_total = 2.5; % [m^3]

% fibers
density = 1150; % [kg/m^3]
length = 1.2e-3; % [m]
diameter = 10e-6; % [m]
volume_fiber = pi*(diameter/2)^2 * length;
weight_fiber = density * volume_fiber;

%% computation
desired_fiber_concentration = desired_nr_fibers / measurement_volume; % [1/m^3]
total_nr_fibers = desired_fiber_concentration * volume_water_total; % [-]

weight_fiber_total_g = weight_fiber * total_nr_fibers * 1000; % [kg]