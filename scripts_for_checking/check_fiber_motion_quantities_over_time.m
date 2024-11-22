clc; clear; close all;
%clearvars -except AllFibers p 
close all
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

addpath("_common\")

%% load data
load('G:\__PRL\processed_data_close_wall\svd_rlowess_ts_31_fk_45\_all_data.mat')

%% delete data at the edges of the volume
avg_z_fitted=mean(z_fitted,2,"omitnan");

fibre_length(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
fibre_curvature(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
omega_spinning(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
omega_tumbling(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
omega_s_x(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
omega_s_y(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
omega_s_z(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
omega_b_x(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
omega_b_y(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
omega_b_z(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
x(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
y(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
z(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
x_fitted(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
y_fitted(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
z_fitted(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
x_dot(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
y_dot(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
z_dot(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
red(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:,:)=NaN;
green(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:,:)=NaN;
blue(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:,:)=NaN;
red_fitted(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:,:)=NaN;
green_fitted(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:,:)=NaN;
blue_fitted(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:,:)=NaN;
L1(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
L2(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;
L3(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:)=NaN;

px(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:,:)=NaN;
py(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:,:)=NaN;
pz(avg_z_fitted<bottom_threshold | avg_z_fitted>top_threshold,:,:)=NaN;

%% delete data based on length and curvature
avg_fibre_length = mean(fibre_length,2,"omitnan");
avg_fibre_curvature = mean(fibre_curvature,2,"omitnan");

fibre_length(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
fibre_curvature(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
omega_spinning(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
omega_tumbling(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
omega_s_x(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
omega_s_y(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
omega_s_z(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
omega_b_x(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
omega_b_y(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
omega_b_z(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
x(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
y(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
z(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
x_fitted(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
y_fitted(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
z_fitted(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
x_dot(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
y_dot(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
z_dot(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
red(avg_fibre_length < min_length | avg_fibre_length > max_length,:,:)=NaN;
green(avg_fibre_length < min_length | avg_fibre_length > max_length,:,:)=NaN;
blue(avg_fibre_length < min_length | avg_fibre_length > max_length,:,:)=NaN;
red_fitted(avg_fibre_length < min_length | avg_fibre_length > max_length,:,:)=NaN;
green_fitted(avg_fibre_length < min_length | avg_fibre_length > max_length,:,:)=NaN;
blue_fitted(avg_fibre_length < min_length | avg_fibre_length > max_length,:,:)=NaN;
L1(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
L2(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;
L3(avg_fibre_length < min_length | avg_fibre_length > max_length,:)=NaN;

px(avg_fibre_length < min_length | avg_fibre_length > max_length,:,:)=NaN;
py(avg_fibre_length < min_length | avg_fibre_length > max_length,:,:)=NaN;
pz(avg_fibre_length < min_length | avg_fibre_length > max_length,:,:)=NaN;


fibre_length(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
fibre_curvature(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
omega_spinning(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
omega_tumbling(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
omega_s_x(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
omega_s_y(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
omega_s_z(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
omega_b_x(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
omega_b_y(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
omega_b_z(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
x(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
y(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
z(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
x_fitted(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
y_fitted(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
z_fitted(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
x_dot(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
y_dot(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
z_dot(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
red(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:,:)=NaN;
green(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:,:)=NaN;
blue(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:,:)=NaN;
red_fitted(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:,:)=NaN;
green_fitted(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:,:)=NaN;
blue_fitted(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:,:)=NaN;
L1(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
L2(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;
L3(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:)=NaN;

px(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:,:)=NaN;
py(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:,:)=NaN;
pz(avg_fibre_curvature < min_curv | avg_fibre_curvature > max_curv,:,:)=NaN;

%% delete data based on std of curvature within each track
std_fibre_curvature = std(fibre_curvature,0,2,'omitnan');

fibre_length(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
fibre_curvature(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
omega_spinning(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
omega_tumbling(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
omega_s_x(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
omega_s_y(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
omega_s_z(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
omega_b_x(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
omega_b_y(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
omega_b_z(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
x(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
y(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
z(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
x_fitted(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
y_fitted(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
z_fitted(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
x_dot(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
y_dot(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
z_dot(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
red(std_fibre_curvature>std_fibre_curvature_threshold,:,:)=NaN;
green(std_fibre_curvature>std_fibre_curvature_threshold,:,:)=NaN;
blue(std_fibre_curvature>std_fibre_curvature_threshold,:,:)=NaN;
red_fitted(std_fibre_curvature>std_fibre_curvature_threshold,:,:)=NaN;
green_fitted(std_fibre_curvature>std_fibre_curvature_threshold,:,:)=NaN;
blue_fitted(std_fibre_curvature>std_fibre_curvature_threshold,:,:)=NaN;
L1(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
L2(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
L3(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;

px(std_fibre_curvature>std_fibre_curvature_threshold,:,:)=NaN;
py(std_fibre_curvature>std_fibre_curvature_threshold,:,:)=NaN;
pz(std_fibre_curvature>std_fibre_curvature_threshold,:,:)=NaN;

%% plot

for ij=1:size(fibre_length,1)
    %index_time=find(~isnan(fibre_length(ij,:)));
    index_time=find(~cellfun('isempty',AllFibers.Centroid(ij,:)')==1);
    timesteps=index_time-index_time(1)+1;
    time=(timesteps-1)*p.dt;

%%% length and curvature
    figure(1); set(gcf,"Position",[1 50 400 800]); clf
    subplot(2,1,1); hold on; grid on; box on;
    plot(fibre_length(ij,~isnan(fibre_length(ij,:))),'.--')
    plot([1 timesteps(end)],[mean(fibre_length(ij,:),'omitnan'),mean(fibre_length(ij,:),'omitnan')],'k-')
    ylabel('$L_f$ [mm]')
    xlabel('time-step')
    legend('Data','Mean','location','best')
%
    subplot(2,1,2); hold on; grid on; box on;
    plot(fibre_curvature(ij,~isnan(fibre_curvature(ij,:))),'.--')
    plot([1 timesteps(end)],[mean(fibre_curvature(ij,:),'omitnan'),mean(fibre_curvature(ij,:),'omitnan')],'k-')
    ylabel('$\kappa$ [1]')
    xlabel('time-step')
    legend('Data','Mean','location','best')

%%% positions, velocities
    figure(2); set(gcf,"Position",[400 50 1000 800]); clf
    subplot(2,3,1); hold on; grid on; box on;
    plot(x(ij,~isnan(fibre_curvature(ij,:))),'.')
    plot(x_fitted(ij,~isnan(fibre_curvature(ij,:))),'-')
    ylabel('$x$ [mm]')
    xlabel('time-step')
    legend('raw','fitted','location','best')

    subplot(2,3,2); hold on; grid on; box on;
    plot(y(ij,~isnan(fibre_curvature(ij,:))),'.')
    plot(y_fitted(ij,~isnan(fibre_curvature(ij,:))),'-')
    ylabel('$y$ [mm]')
    xlabel('time-step')
    legend('raw','fitted','location','best')

    subplot(2,3,3); hold on; grid on; box on;
    plot(z(ij,~isnan(fibre_curvature(ij,:))),'.')
    plot(z_fitted(ij,~isnan(fibre_curvature(ij,:))),'-')
    ylabel('$z$ [mm]')
    xlabel('time-step')
    legend('raw','fitted','location','best')
%
    subplot(2,3,4); hold on; grid on; box on;
    plot(x_dot(ij,~isnan(fibre_curvature(ij,:))),'.-')
    ylabel('$\dot{x}$ [mm/s]')
    xlabel('time-step')

    subplot(2,3,5); hold on; grid on; box on;
    plot(y_dot(ij,~isnan(fibre_curvature(ij,:))),'.-')
    ylabel('$\dot{y}$ [mm/s]')
    xlabel('time-step')

    subplot(2,3,6); hold on; grid on; box on;
    plot(z_dot(ij,~isnan(fibre_curvature(ij,:))),'.-')
    ylabel('$\dot{z}$ [mm/s]')
    xlabel('time-step')

%%% vectors, rotation rates
    figure(3); set(gcf,"Position",[800 50 1000 800]); clf
    subplot(2,3,1); hold on; grid on; box on;
    plot(red(ij,~isnan(fibre_curvature(ij,:)),1),'r.')
    plot(red(ij,~isnan(fibre_curvature(ij,:)),2),'g.')
    plot(red(ij,~isnan(fibre_curvature(ij,:)),3),'b.')
    plot(red_fitted(ij,~isnan(fibre_curvature(ij,:)),1),'r-')
    plot(red_fitted(ij,~isnan(fibre_curvature(ij,:)),2),'g-')
    plot(red_fitted(ij,~isnan(fibre_curvature(ij,:)),3),'b-')
    ylabel('$e_{11},\ e_{12},\ e_{13}$')
    xlabel('time-step')
    p1=plot([],[],'k.');
    p2=plot([],[],'k-');
    legend([p1,p2],{'raw','fitted'},'location','best')

    subplot(2,3,2); hold on; grid on; box on;
    plot(green(ij,~isnan(fibre_curvature(ij,:)),1),'r.')
    plot(green(ij,~isnan(fibre_curvature(ij,:)),2),'g.')
    plot(green(ij,~isnan(fibre_curvature(ij,:)),3),'b.')
    plot(green_fitted(ij,~isnan(fibre_curvature(ij,:)),1),'r-')
    plot(green_fitted(ij,~isnan(fibre_curvature(ij,:)),2),'g-')
    plot(green_fitted(ij,~isnan(fibre_curvature(ij,:)),3),'b-')
    ylabel('$e_{21},\ e_{22},\ e_{23}$')
    xlabel('time-step')
    p1=plot([],[],'k.');
    p2=plot([],[],'k-');
    legend([p1,p2],{'raw','fitted'},'location','best')

    subplot(2,3,3); hold on; grid on; box on;
    plot(blue(ij,~isnan(fibre_curvature(ij,:)),1),'r.')
    plot(blue(ij,~isnan(fibre_curvature(ij,:)),2),'g.')
    plot(blue(ij,~isnan(fibre_curvature(ij,:)),3),'b.')
    plot(blue_fitted(ij,~isnan(fibre_curvature(ij,:)),1),'r-')
    plot(blue_fitted(ij,~isnan(fibre_curvature(ij,:)),2),'g-')
    plot(blue_fitted(ij,~isnan(fibre_curvature(ij,:)),3),'b-')
    ylabel('$e_{31},\ e_{32},\ e_{33}$')
    xlabel('time-step')
    p1=plot([],[],'k.');
    p2=plot([],[],'k-');
    legend([p1,p2],{'raw','fitted'},'location','best')

    subplot(2,3,4); hold on; grid on; box on;
    plot(omega_spinning(ij,~isnan(fibre_curvature(ij,:))),'.-')
    hold on
    plot([1 timesteps(end)],[mean(omega_spinning(ij,:),'omitnan'),mean(omega_spinning(ij,:),'omitnan')],'k-')
    ylabel('$\Omega_s$ [deg/frame]')
    xlabel('time-step')
    legend('Data','Mean','location','best')

    subplot(2,3,5); hold on; grid on; box on;
    plot(omega_tumbling(ij,~isnan(fibre_curvature(ij,:))),'.-')
    hold on
    plot([1 timesteps(end)],[mean(omega_tumbling(ij,:),'omitnan'),mean(omega_tumbling(ij,:),'omitnan')],'k-')
    ylabel('$\Omega_t$ [deg/frame]')
    xlabel('time-step')
    legend('Data','Mean','location','best')

    
    pause
end