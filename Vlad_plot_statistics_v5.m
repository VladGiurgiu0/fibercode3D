clc; clear; close all;
%clearvars -except AllFibers p 
close all
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

addpath("_common\")

which_data='near_wall';

switch which_data
    case 'near_wall'
        
        fname.main='D:\__PRL\processed_data_close_wall\rlowess_ts_15_fk_45\';
        fname.subfolder_root = 'Loop=';
        fname.ending = '\4_Quantities_Refined_fibers\AllFibers_Only_data.mat';
        
        reload_data_from_loops=0;
        starting_loop_number = 1;
        ending_loop_number=75;

        if reload_data_from_loops==0
            % load data
            load(strcat(fname.main,'_all_data.mat'))
            reload_data_from_loops=0;
        end

        k0=2.62; % reference curvature for a 1.2mm long fibre [1/mm]
        l_nominal = 1.2; % nominal length of fibers [mm]
        
        %%% the wall-normal has to be pointing up in the MART
        ref_wall = 'top'; % which wall should be used as reference: 'bottom' or 'top'
        %z_pos_1_mart = 34.8; % distance between the first bottom plane of MART and the ref_wall [mm]
        z_pos_1_mart = 16.15; % should be 16, because that is where the wall is, but in the PRL paper it was shown with 16.15
         
        
        % discard data based on curvature fluctuations within the track
        std_fibre_curvature_threshold=0.05;  % threshold for the std of fibre curvature within each track
        
        % discard data based on proximity to the edges of the volume in wall-normal
        % direction
        top_threshold=20; % [mm]
        bottom_threshold=1; % [mm]
        
        % discard data based on track average fibre length and curvature
        min_length = 0.7; % [mm]
        max_length = 1.4; % [mm]
        
        min_curv = 0.2; % [-]
        max_curv = 0.6; % [-]
        
        % number of bins in wall-normal direction
        numBins_z=12;
        alpha=0.01; % significance level for confidence interval computation
        which_bins = [7,9,12]; % for which three bins to show convergence
        
        % bin width in wall-normal direction to measure concentration
        bin_width_concentration=1; % [mm]

    case 'center'

        fname.main = 'F:\Fiber23May2023b\_processed_data_for_PRL\svd_rlowess_ts_30_fk_45\';
        fname.subfolder_root = 'Loop=';
        fname.ending = '\4_Quantities_Refined_fibers\AllFibers_Only_data.mat';

        reload_data_from_loops=0;
        starting_loop_number = 1;
        ending_loop_number=133;

        if reload_data_from_loops==0
            % load data
            load(strcat(fname.main,'_all_data.mat'))
            reload_data_from_loops=0;
        end

        k0=2.62; % reference curvature for a 1.2mm long fibre [1/mm]
        l_nominal = 1.2; % nominal length of fibers [mm]
        
        %%% the wall-normal has to be pointing up in the MART
        ref_wall = 'bottom'; % which wall should be used as reference: 'bottom' or 'top'
        z_pos_1_mart = 34.8; % distance between the first bottom plane of MART and the ref_wall [mm]
        %z_pos_1_mart = 16.3;
        
        % discard data based on curvature fluctuations within the track
        std_fibre_curvature_threshold=0.05;  % threshold for the std of fibre curvature within each track
        
        % discard data based on proximity to the edges of the volume in wall-normal
        % direction
        top_threshold=10; % [mm]
        bottom_threshold=1; % [mm]
        
        % discard data based on track average fibre length and curvature
        min_length = 0.7; % [mm]
        max_length = 1.4; % [mm]
        
        min_curv = 0.2; % [-]
        max_curv = 0.6; % [-]
        
        % number of bins in wall-normal direction
        numBins_z=3;
        alpha=0.01; % significance level for confidence interval computation
        which_bins = [1,2,3]; % for which three bins to show convergence
        
        % bin width in wall-normal direction to measure concentration
        bin_width_concentration=1; % [mm]

end


%% extract data
if reload_data_from_loops==1
for ll = starting_loop_number:ending_loop_number
    load(strcat(fname.main,fname.subfolder_root,num2str(ll),fname.ending))
    [~,~,tau]=compute_wall_scales(p.Re_tau,p.Temperature);

    disp(ll)
    
    if ll==starting_loop_number
        nr_fibres_old=0;
    else
        nr_fibres_old = nr_fibres_old +nr_fibres_new;
    end

    nr_fibres_new = size(AllFibers.Centroid(:,1),1);

    for ij=1:nr_fibres_new
        ik = ij+nr_fibres_old;
        it=find(~cellfun('isempty',AllFibers.Centroid(ij,:)))';
    
        track_length(ik)=numel(find(~cellfun('isempty',AllFibers.Centroid(ij,:)))');
        
        fibre_length(ik,it)=cell2mat(AllFibers.Length(ij,it))*p.dx;             % [mm]
        fibre_curvature(ik,it)=cell2mat(AllFibers.Curvature(ij,it))/p.dx / k0; % dimensionless curvature [1]
    
        omega_spinning(ik,it)=cell2mat(AllFibers.Spinning_rate_RotationMatrix(ij,it)); % spinning rate [deg/frame]
        omega_tumbling(ik,it)=cell2mat(AllFibers.Tumbling_rate_RotationMatrix(ij,it)); % tumbling rate [deg/frame]

        omega_s_x(ik,it)=cell2mat(AllFibers.Omega_s_x_RotationMatrix(ij,it)); % rotation rate around laboratory x [deg/frame]
        omega_s_y(ik,it)=cell2mat(AllFibers.Omega_s_y_RotationMatrix(ij,it)); % rotation rate around laboratory y [deg/frame]
        omega_s_z(ik,it)=cell2mat(AllFibers.Omega_s_z_RotationMatrix(ij,it)); % rotation rate around laboratory z [deg/frame]

        omega_b_x(ik,it)=cell2mat(AllFibers.Omega_b_x_RotationMatrix(ij,it)); % rotation rate around fibre x - spinning [deg/frame]
        omega_b_y(ik,it)=cell2mat(AllFibers.Omega_b_y_RotationMatrix(ij,it)); % rotation rate around fibre y - tumbling 2 [deg/frame]
        omega_b_z(ik,it)=cell2mat(AllFibers.Omega_b_z_RotationMatrix(ij,it)); % rotation rate around fibre z - tumbling 3[deg/frame]
    
        x(ik,it)=cell2mat(AllFibers.x(ij,it)); % [mm]
        y(ik,it)=cell2mat(AllFibers.y(ij,it)); % [mm]
        z(ik,it)=cell2mat(AllFibers.z(ij,it)); % [mm]
    
        x_fitted(ik,it)=cell2mat(AllFibers.x_fitted(ij,it)); % [mm]
        y_fitted(ik,it)=cell2mat(AllFibers.y_fitted(ij,it)); % [mm]
        z_fitted(ik,it)=cell2mat(AllFibers.z_fitted(ij,it)); % [mm]
    
        x_dot(ik,it)=cell2mat(AllFibers.x_dot(ij,it)); % [mm/s]
        y_dot(ik,it)=cell2mat(AllFibers.y_dot(ij,it)); % [mm/s]
        z_dot(ik,it)=cell2mat(AllFibers.z_dot(ij,it)); % [mm/s]

        x_dot_dot(ik,it)=cell2mat(AllFibers.x_dot_dot(ij,it)); % [mm/s]
        y_dot_dot(ik,it)=cell2mat(AllFibers.y_dot_dot(ij,it)); % [mm/s]
        z_dot_dot(ik,it)=cell2mat(AllFibers.z_dot_dot(ij,it)); % [mm/s]

        for tt=it'
            red(ik,tt,:)=cell2mat(AllFibers.red_tensor(ij,tt));         % red vector - lowest moment of inertia
            green(ik,tt,:)=cell2mat(AllFibers.green_tensor(ij,tt));     % green vector - medium moment of inertia
            blue(ik,tt,:)=cell2mat(AllFibers.blue_tensor(ij,tt));       % blue vector - highest moment of inertia

            red_fitted(ik,tt,:)=cell2mat(AllFibers.red_tensor_fitted(ij,tt));         % red vector - lowest moment of inertia
            green_fitted(ik,tt,:)=cell2mat(AllFibers.green_tensor_fitted(ij,tt));     % green vector - medium moment of inertia
            blue_fitted(ik,tt,:)=cell2mat(AllFibers.blue_tensor_fitted(ij,tt));       % blue vector - highest moment of inertia

            px(ik,tt,:)=cell2mat(AllFibers.px_Tensor(ij,tt));               % x component of the fitted polynomial
            py(ik,tt,:)=cell2mat(AllFibers.py_Tensor(ij,tt));               % x component of the fitted polynomial
            pz(ik,tt,:)=cell2mat(AllFibers.pz_Tensor(ij,tt));               % x component of the fitted polynomial

        
            if ll>22
            % lengths
                L1(ik,tt)=AllFibers.PrincipalAxisLength{ij,tt}(1);  % longest axis length
                L2(ik,tt)=AllFibers.PrincipalAxisLength{ij,tt}(2);  % medium axis length
                L3(ik,tt)=AllFibers.PrincipalAxisLength{ij,tt}(3);  % shortest axis length
            end
        end

    end

end

%% change 0 to NaN in data
fibre_length(fibre_length==0)=NaN;
fibre_curvature(fibre_curvature==0)=NaN;
omega_spinning(omega_spinning==0)=NaN;
omega_tumbling(omega_tumbling==0)=NaN; 
omega_s_x(omega_s_x==0)=NaN;
omega_s_y(omega_s_y==0)=NaN;
omega_s_z(omega_s_z==0)=NaN;
omega_b_x(omega_b_x==0)=NaN;
omega_b_y(omega_b_y==0)=NaN;
omega_b_z(omega_b_z==0)=NaN;
x(x==0)=NaN;
y(y==0)=NaN;
z(z==0)=NaN;
x_fitted(x_fitted==0)=NaN;
y_fitted(y_fitted==0)=NaN;
z_fitted(z_fitted==0)=NaN;
x_dot(x_dot==0)=NaN;
y_dot(y_dot==0)=NaN;
z_dot(z_dot==0)=NaN;
x_dot_dot(x_dot_dot==0)=NaN;
y_dot_dot(y_dot_dot==0)=NaN;
z_dot_dot(z_dot_dot==0)=NaN;

red(red==0)=NaN;
green(green==0)=NaN;
blue(blue==0)=NaN;
red_fitted(red_fitted==0)=NaN;
green_fitted(green_fitted==0)=NaN;
blue_fitted(blue_fitted==0)=NaN;

L1(L1==0)=NaN;
L2(L2==0)=NaN;
L3(L3==0)=NaN;

px(px==0)=NaN;
py(py==0)=NaN;
pz(pz==0)=NaN;

%% save the data
save(strcat(fname.main,'_all_data.mat'))

% else
%     %% load data
%     load(strcat(fname.main,'_all_data.mat'))
%     reload_data_from_loops=0;

end

%% delete data at the edges of the volume

fibre_length(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
fibre_curvature(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;

omega_spinning(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
omega_tumbling(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
omega_s_x(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
omega_s_y(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
omega_s_z(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
omega_b_x(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
omega_b_y(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
omega_b_z(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;

x(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
y(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
z(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
x_fitted(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
y_fitted(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
z_fitted(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
x_dot(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
y_dot(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
z_dot(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
x_dot_dot(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
y_dot_dot(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
z_dot_dot(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;

red(z_fitted<bottom_threshold | z_fitted>top_threshold,:)=NaN;
green(z_fitted<bottom_threshold | z_fitted>top_threshold,:)=NaN;
blue(z_fitted<bottom_threshold | z_fitted>top_threshold,:)=NaN;
red_fitted(z_fitted<bottom_threshold | z_fitted>top_threshold,:)=NaN;
green_fitted(z_fitted<bottom_threshold | z_fitted>top_threshold,:)=NaN;
blue_fitted(z_fitted<bottom_threshold | z_fitted>top_threshold,:)=NaN;

L1(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
L2(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;
L3(z_fitted<bottom_threshold | z_fitted>top_threshold)=NaN;

px(z_fitted<bottom_threshold | z_fitted>top_threshold,:)=NaN;
py(z_fitted<bottom_threshold | z_fitted>top_threshold,:)=NaN;
pz(z_fitted<bottom_threshold | z_fitted>top_threshold,:)=NaN;

avg_z_fitted=mean(z_fitted,2,"omitnan");

%% delete data based on length and curvature
avg_fibre_curvature = mean(fibre_curvature,2,"omitnan");

fibre_length(fibre_length < min_length | fibre_length > max_length)=NaN;
fibre_curvature(fibre_length < min_length | fibre_length > max_length)=NaN;

omega_spinning(fibre_length < min_length | fibre_length > max_length)=NaN;
omega_tumbling(fibre_length < min_length | fibre_length > max_length)=NaN;
omega_s_x(fibre_length < min_length | fibre_length > max_length)=NaN;
omega_s_y(fibre_length < min_length | fibre_length > max_length)=NaN;
omega_s_z(fibre_length < min_length | fibre_length > max_length)=NaN;
omega_b_x(fibre_length < min_length | fibre_length > max_length)=NaN;
omega_b_y(fibre_length < min_length | fibre_length > max_length)=NaN;
omega_b_z(fibre_length < min_length | fibre_length > max_length)=NaN;

x(fibre_length < min_length | fibre_length > max_length)=NaN;
y(fibre_length < min_length | fibre_length > max_length)=NaN;
z(fibre_length < min_length | fibre_length > max_length)=NaN;
x_fitted(fibre_length < min_length | fibre_length > max_length)=NaN;
y_fitted(fibre_length < min_length | fibre_length > max_length)=NaN;
z_fitted(fibre_length < min_length | fibre_length > max_length)=NaN;
x_dot(fibre_length < min_length | fibre_length > max_length)=NaN;
y_dot(fibre_length < min_length | fibre_length > max_length)=NaN;
z_dot(fibre_length < min_length | fibre_length > max_length)=NaN;
x_dot_dot(fibre_length < min_length | fibre_length > max_length)=NaN;
y_dot_dot(fibre_length < min_length | fibre_length > max_length)=NaN;
z_dot_dot(fibre_length < min_length | fibre_length > max_length)=NaN;
red(fibre_length < min_length | fibre_length > max_length,:)=NaN;
green(fibre_length < min_length | fibre_length > max_length,:)=NaN;
blue(fibre_length < min_length | fibre_length > max_length,:)=NaN;
red_fitted(fibre_length < min_length | fibre_length > max_length,:)=NaN;
green_fitted(fibre_length < min_length | fibre_length > max_length,:)=NaN;
blue_fitted(fibre_length < min_length | fibre_length > max_length,:)=NaN;
L1(fibre_length < min_length | fibre_length > max_length)=NaN;
L2(fibre_length < min_length | fibre_length > max_length)=NaN;
L3(fibre_length < min_length | fibre_length > max_length)=NaN;

px(fibre_length < min_length | fibre_length > max_length,:)=NaN;
py(fibre_length < min_length | fibre_length > max_length,:)=NaN;
pz(fibre_length < min_length | fibre_length > max_length,:)=NaN;


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

avg_fibre_length = mean(fibre_length,2,"omitnan");


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
x_dot_dot(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
y_dot_dot(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;
z_dot_dot(std_fibre_curvature>std_fibre_curvature_threshold,:)=NaN;

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

%% pdf avg, std of length and curvature

avg_fibre_length = mean(fibre_length,2,"omitnan");
std_fibre_length = std(fibre_length,0,2,'omitnan');

avg_fibre_curvature = mean(fibre_curvature,2,"omitnan");
std_fibre_curvature = std(fibre_curvature,0,2,'omitnan');

% plot histograms of avg curvature and length
figure(); hold all;

subplot(1,2,1); hold on
nbins=10;
[centers_len,count_len]=compute_pdf(avg_fibre_length,'BinNumber',nbins,'pdf');
plot(centers_len,count_len,'b.-','MarkerSize',15)
set(gca,'FontSize',20)
xlabel('$L_f$ [mm]')
ylabel('p.d.f')
xlim([0 2])
hold on; box on; grid minor;

subplot(1,2,2); hold on;
nbins=10;
[centers_k,count_k]=compute_pdf(avg_fibre_curvature,'BinNumber',nbins,'probability');
plot(centers_k,count_k,'b.-','MarkerSize',15)
set(gca,'FontSize',20)
xlabel('$\kappa^*$ [-]')
ylabel('p.d.f')
xlim([0 1])
hold on; box on; grid minor;

% plot histogram of std curvature
figure(); hold on; grid on; box on;
title('Histogram of std. of fibre curvatures within each track')
[centers_std_k,count_std_k]=compute_pdf(std_fibre_curvature,'BinWidth',0.005,'probability');
plot(centers_std_k,count_std_k,'k.-','MarkerSize',20); 
set(gca,'FontSize',15)
ylabel('pdf')
xlim([0 0.3])
%ylim([0 0.3])

% plot histogram of track lengths
figure(); hold on; grid on; box on;
title('Histogram of fibre track lengths in timesteps [-]')
histogram(track_length);



%% compute wall-distance and Kolmogorov scales based on wall-distance
[scales]=Vlad_compute_scales_fct(p.Re_tau,p.Temperature);

switch ref_wall
    case 'bottom'
        mean_wall_distance_fibre=mean(z_fitted,2,"omitnan") + z_pos_1_mart;
        wall_distance_fibre = z_fitted + z_pos_1_mart;
    case 'top'
        mean_wall_distance_fibre=-mean(z_fitted,2,"omitnan") + z_pos_1_mart;
        wall_distance_fibre = -z_fitted + z_pos_1_mart;
end


local_l_kol_fibre = scales.l_kol(knnsearch(scales.y,mean_wall_distance_fibre));
local_t_kol_fibre = scales.t_kol(knnsearch(scales.y,mean_wall_distance_fibre));

% plot to check
figure(); hold on; box on; grid on;
[~,idx_sort]=sort(mean_wall_distance_fibre);
subplot(1,2,1); grid on;
plot(mean_wall_distance_fibre(idx_sort),1e3*local_t_kol_fibre(idx_sort),'r-')
ylabel('$\tau_{\eta}$ [ms]')
xlabel('wall-distance [mm]')
subplot(1,2,2); grid on;
plot(mean_wall_distance_fibre(idx_sort),1e3*local_l_kol_fibre(idx_sort),'b-')
xlabel('wall-distance [mm]')
ylabel('$\eta$ [mm]')

figure(); hold on; box on; grid on;
histogram(mean_wall_distance_fibre,100);
title('Mean distance between fibre and wall')
xlabel('Distance between fibre and wall [mm]')

%% concentration
% linearly spaced bins
wall_dist=wall_distance_fibre(:);

[concentration_counts,concentration_edges]=histcounts(wall_dist,'BinWidth',bin_width_concentration,'Normalization','count');
concentration_centers=(concentration_edges(1:end-1)+concentration_edges(2:end))/2; clearvars concentration_edges
figure();
plot(concentration_centers,concentration_counts,'k.-','MarkerSize',20)
title('Concentration');xlabel('$y$ [mm]'); ylabel('Nr. fibres in each bin [-]'); grid on;



%% rotation rates 
omega_tumbling_a = deg2rad(abs(omega_tumbling))/p.dt; % [rad/s]
omega_spinning_a = deg2rad(abs(omega_spinning))/p.dt; % [rad/s]

omega_s_x = deg2rad(abs(omega_s_x))/p.dt; % [rad/s]
omega_s_y = deg2rad(abs(omega_s_y))/p.dt; % [rad/s]
omega_s_z = deg2rad(abs(omega_s_z))/p.dt; % [rad/s]

omega_b_x = deg2rad(abs(omega_b_x))/p.dt; % [rad/s]
omega_b_y = deg2rad(abs(omega_b_y))/p.dt; % [rad/s]
omega_b_z = deg2rad(abs(omega_b_z))/p.dt; % [rad/s]

omega_tumbling_kol = omega_tumbling_a .* local_t_kol_fibre;
omega_spinning_kol = omega_spinning_a .* local_t_kol_fibre;

omega_tumbling_visc = omega_tumbling_a * scales.viscous_time;
omega_spinning_visc = omega_spinning_a * scales.viscous_time;

mean_tumbling_kol = mean(reshape(omega_tumbling_kol,1,[]),'omitnan');
mean_spinning_kol = mean(reshape(omega_spinning_kol,1,[]),'omitnan');

mean_tumbling_visc = mean(reshape(omega_tumbling_visc,1,[]),'omitnan');
mean_spinning_visc = mean(reshape(omega_spinning_visc,1,[]),'omitnan');

mean_squared_tumbling_kol = mean(reshape(omega_tumbling_kol,1,[]).^2,'omitnan');
mean_squared_spinning_kol = mean(reshape(omega_spinning_kol,1,[]).^2,'omitnan');

mean_squared_tumbling_visc = mean(reshape(omega_tumbling_visc,1,[]).^2,'omitnan');
mean_squared_spinning_visc = mean(reshape(omega_spinning_visc,1,[]).^2,'omitnan');

var_tumbling_kol=mean(reshape(omega_tumbling_kol-mean_tumbling_kol,1,[]).^2,'omitnan');
var_spinning_kol=mean(reshape(omega_spinning_kol-mean_spinning_kol,1,[]).^2,'omitnan');

%% local l/eta over y+
% replenish the data
omega_tumbling_visc = omega_tumbling_a * scales.viscous_time;

l_eta = l_nominal *1e-3 * local_l_kol_fibre.^(-1);
l_eta = repmat(l_eta,1,size(omega_tumbling_kol,2));
l_eta(isnan(omega_tumbling_visc))=NaN;
l_eta= reshape(l_eta,1,[]);

wall_dist = wall_distance_fibre;
wall_dist= reshape(wall_dist,1,[]);

[bin_z_Centers_rot_visc,mean_l_eta,ci_z_Centers_rot_visc,ci_mean_l_eta]=Vlad_bin_mean(wall_dist,l_eta,numBins_z,'log',alpha);


%% plot mean velocity
wall_dist = wall_distance_fibre;
wall_dist(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];

x_dot_b = x_dot;
y_dot_b = y_dot;
z_dot_b = z_dot;

x_dot_b(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];
y_dot_b(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];
z_dot_b(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];

% bin velocities
%numBins_z=10;
wall_dist= reshape(wall_dist,1,[]);
x_velocities = reshape(x_dot_b,1,[]);
y_velocities = reshape(y_dot_b,1,[]);
z_velocities = reshape(z_dot_b,1,[]);

[~, edges] = histcounts(wall_dist, numBins_z);
binIndices = discretize(wall_dist, edges);

binned_x_velocities = cell(1, numBins_z);
binned_y_velocities = cell(1, numBins_z);
binned_z_velocities = cell(1, numBins_z);

mean_x_velocities =[]; mean_y_velocities =[]; mean_z_velocities =[];

for i = 1:numBins_z
    binned_x_velocities{i} = x_velocities(binIndices == i);
    binned_y_velocities{i} = y_velocities(binIndices == i);
    binned_z_velocities{i} = z_velocities(binIndices == i);

    mean_x_velocities(i)=mean(binned_x_velocities{i},'omitnan');
    mean_y_velocities(i)=mean(binned_y_velocities{i},'omitnan');
    mean_z_velocities(i)=mean(binned_z_velocities{i},'omitnan');
end

bin_z_Centers_vel = (edges(1:end-1) + edges(2:end)) / 2;
% plot
figure();
subplot(1,3,1); grid on; box on; hold on;
plot(bin_z_Centers_vel,mean_x_velocities,'.-')
xlabel('$z$ [mm]')
ylabel('$\dot{x}$ [mm/s]')
title('Stream-wise velocity')

subplot(1,3,2); grid on; box on; hold on;
plot(bin_z_Centers_vel,mean_y_velocities,'.-')
xlabel('$z$ [mm]')
ylabel('$\dot{y}$ [mm/s]')
title('Span-wise velocity')

subplot(1,3,3); grid on; box on; hold on;
plot(bin_z_Centers_vel,mean_z_velocities,'.-')
xlabel('$z$ [mm]')
ylabel('$\dot{z}$ [mm/s]')
title('Wall-normal velocity')

%% mean velocities
% replenish the data
wall_dist = wall_distance_fibre;

% plot mean square tumbling and spinning over wall-distance
%numBins_z=12;
wall_dist= reshape(wall_dist,1,[]);

[bin_z_Centers_velocity,mean_velocity_x,ci_z_Centers_velocity,ci_mean_velocity_x]=Vlad_bin_mean(wall_dist,x_dot(:),numBins_z,'log',alpha);
[bin_z_Centers_velocity,mean_velocity_y,ci_z_Centers_velocity,ci_mean_velocity_y]=Vlad_bin_mean(wall_dist,y_dot(:),numBins_z,'log',alpha);
[bin_z_Centers_velocity,mean_velocity_z,ci_z_Centers_velocity,ci_mean_velocity_z]=Vlad_bin_mean(wall_dist,z_dot(:),numBins_z,'log',alpha);


%% mean rotation rates (not squared) over y+
omega_1_visc = omega_spinning_a * scales.viscous_time;
omega_2_visc = omega_b_y * scales.viscous_time;
omega_3_visc = omega_b_z * scales.viscous_time;

omega_x_visc = omega_s_x * scales.viscous_time;
omega_y_visc = omega_s_y * scales.viscous_time;
omega_z_visc = omega_s_z * scales.viscous_time;


wall_dist = wall_distance_fibre;
wall_dist= reshape(wall_dist,1,[]);

% lab reference frame
[bin_z_Centers_simple_omega_visc,mean_simple_omega_x_visc,ci_z_Centers_simple_omega_visc,ci_mean_simple_omega_x_visc]=...
            Vlad_bin_mean(wall_dist,reshape(omega_x_visc,1,[]),numBins_z,'log',alpha);
[bin_z_Centers_simple_omega_visc,mean_simple_omega_y_visc,ci_z_Centers_simple_omega_visc,ci_mean_simple_omega_y_visc]=...
            Vlad_bin_mean(wall_dist,reshape(omega_y_visc,1,[]),numBins_z,'log',alpha);
[bin_z_Centers_simple_omega_visc,mean_simple_omega_z_visc,ci_z_Centers_simple_omega_visc,ci_mean_simple_omega_z_visc]=...
            Vlad_bin_mean(wall_dist,reshape(omega_z_visc,1,[]),numBins_z,'log',alpha);

% fibre reference frame
[bin_z_Centers_simple_omega_visc,mean_simple_omega_1_visc,ci_z_Centers_simple_omega_visc,ci_mean_simple_omega_1_visc]=...
            Vlad_bin_mean(wall_dist,reshape(omega_1_visc,1,[]),numBins_z,'log',alpha);
[bin_z_Centers_simple_omega_visc,mean_simple_omega_2_visc,ci_z_Centers_simple_omega_visc,ci_mean_simple_omega_2_visc]=...
            Vlad_bin_mean(wall_dist,reshape(omega_2_visc,1,[]),numBins_z,'log',alpha);
[bin_z_Centers_simple_omega_visc,mean_simple_omega_3_visc,ci_z_Centers_simple_omega_visc,ci_mean_simple_omega_3_visc]=...
            Vlad_bin_mean(wall_dist,reshape(omega_3_visc,1,[]),numBins_z,'log',alpha);


%% tumbling spinning scaled by Kolmogorov time over y+
% replenish the data
omega_tumbling_kol = omega_tumbling_a .* local_t_kol_fibre;
omega_spinning_kol = omega_spinning_a .* local_t_kol_fibre;
omega_s_x_kol = omega_s_x .* local_t_kol_fibre;
omega_s_y_kol = omega_s_y .* local_t_kol_fibre;
omega_s_z_kol = omega_s_z .* local_t_kol_fibre;

omega_b_x_kol = omega_b_x .* local_t_kol_fibre;
omega_b_y_kol = omega_b_y .* local_t_kol_fibre;
omega_b_z_kol = omega_b_z .* local_t_kol_fibre;

omega_full_s_kol  = sqrt(omega_s_x_kol.^2 + omega_s_y_kol.^2 + omega_s_z_kol.^2);
omega_full_b_kol  = sqrt(omega_b_x_kol.^2 + omega_b_y_kol.^2 + omega_b_z_kol.^2);


wall_dist = wall_distance_fibre;

% plot mean square tumbling and spinning over wall-distance
%numBins_z=10;
wall_dist= reshape(wall_dist,1,[]);
square_omega_tumbling_kol = reshape(omega_tumbling_kol,1,[]).^2;
square_omega_spinning_kol = reshape(omega_spinning_kol,1,[]).^2;
%square_omega_full_body_rot_rate_kol = reshape(sqrt(omega_tumbling_kol.^2 + omega_spinning_kol.^2),1,[]);
square_omega_full_body_rot_rate_s_kol = reshape(omega_full_s_kol,1,[]).^2;
square_omega_full_body_rot_rate_b_kol = reshape(omega_full_b_kol,1,[]).^2;

[~, edges] = histcounts(wall_dist, numBins_z);
binIndices = discretize(wall_dist, edges);

binned_tumbl = cell(1, numBins_z);
binned_spinn = cell(1, numBins_z);
binned_full_rot_rate = cell(1, numBins_z);
mean_tumbl=[];
mean_spinn=[];
mean_full_rot_rate = [];


for i = 1:numBins_z
    binned_tumbl{i} = square_omega_tumbling_kol(binIndices == i);
    binned_spinn{i} = square_omega_spinning_kol(binIndices == i);
    binned_full_rot_rate{i} = square_omega_full_body_rot_rate_s_kol(binIndices == i);

    mean_tumbl(i)=mean(binned_tumbl{i}, 'omitnan');
    mean_spinn(i)=mean(binned_spinn{i}, 'omitnan');
    mean_full_rot_rate(i)=mean(binned_full_rot_rate{i}, 'omitnan');
end

bin_z_Centers_rot = (edges(1:end-1) + edges(2:end)) / 2;
% plot
figure(); grid on; box on; hold on;
plot(bin_z_Centers_rot/1e3/scales.viscous_length,mean_tumbl,'.-')
plot(bin_z_Centers_rot/1e3/scales.viscous_length,mean_spinn,'.-')
xlabel('$y^+$ [-]')
ylabel('$\langle \omega_i^2 \rangle \tau_{\eta}^2$ [-]')
title('Mean square rotation rates over wall-normal distance')
legend('$\langle \omega_t^2 \rangle \tau_{\eta}$','$\langle \omega_s^2 \rangle \tau_{\eta}$')
set(gca,'Yscale','log','Xscale','log')

%% tumbling and spinning scaled by viscous time over y+
% replenish the data
omega_tumbling_visc = omega_tumbling_a * scales.viscous_time;
omega_spinning_visc = omega_spinning_a * scales.viscous_time;
wall_dist = wall_distance_fibre;

% plot mean square tumbling and spinning over wall-distance
%numBins_z=12;
wall_dist= reshape(wall_dist,1,[]);
square_omega_tumbling_visc = reshape(omega_tumbling_visc,1,[]).^2;
square_omega_spinning_visc = reshape(omega_spinning_visc,1,[]).^2;

mean_squared_tumbling_visc = mean(reshape(omega_tumbling_visc,1,[]).^2,'omitnan');
mean_squared_spinning_visc = mean(reshape(omega_spinning_visc,1,[]).^2,'omitnan');

[bin_z_Centers_rot_visc,mean_tumbl_visc,ci_z_Centers_rot_visc,ci_mean_tumbl_visc]=Vlad_bin_mean(wall_dist,square_omega_tumbling_visc,numBins_z,'log',alpha);
[bin_z_Centers_rot_visc,mean_spinn_visc,ci_z_Centers_rot_visc,ci_mean_spinn_visc]=Vlad_bin_mean(wall_dist,square_omega_spinning_visc,numBins_z,'log',alpha);

% plot
figure(); grid on; box on; hold on;
%plot(bin_z_Centers_rot_visc/1e3/scales.viscous_length,mean_tumbl_visc,'.-')
errorbar(bin_z_Centers_rot_visc/1e3/scales.viscous_length,mean_tumbl_visc,...
    mean_tumbl_visc-ci_mean_tumbl_visc(:,1),mean_tumbl_visc-ci_mean_tumbl_visc(:,1),...
    (bin_z_Centers_rot_visc' - ci_z_Centers_rot_visc(:,1))/1e3/scales.viscous_length,...
    (bin_z_Centers_rot_visc' - ci_z_Centers_rot_visc(:,1))/1e3/scales.viscous_length);
%plot(bin_z_Centers_rot_visc/1e3/scales.viscous_length,mean_spinn_visc,'.-')
errorbar(bin_z_Centers_rot_visc/1e3/scales.viscous_length,mean_spinn_visc,...
    mean_spinn_visc-ci_mean_spinn_visc(:,1),mean_spinn_visc-ci_mean_spinn_visc(:,1),...
    (bin_z_Centers_rot_visc' - ci_z_Centers_rot_visc(:,1))/1e3/scales.viscous_length,...
    (bin_z_Centers_rot_visc' - ci_z_Centers_rot_visc(:,1))/1e3/scales.viscous_length);
xlabel('$y^+$ [-]')
ylabel('$\langle \omega_i^2 \rangle \tau^2$ [-]')
title('Mean square rotation rates over wall-normal distance')
legend('$\langle \omega_t^2 \rangle \tau$','$\langle \omega_s^2 \rangle \tau$')
set(gca,'Yscale','lin','Xscale','log')

%% convergence of tumbl and spinn (visc) over y plus within 3 bins
%nr_samples_bins = round(logspace(log10(10),log10(size(omega_spinning,1))-2,30)); % logarithmically spaced
nr_samples_bins = round(linspace(100,size(omega_spinning,1)-2,30));                  % linearly spaced
% replenish the data
%omega_tumbling_visc = omega_tumbling_a * scales.viscous_time;
%omega_spinning_visc = omega_spinning_a * scales.viscous_time;
wall_dist = wall_distance_fibre;

mean_tumbl_visc_conv_bins=[];
mean_spinn_visc_conv_bins=[];
nr_traj_used_conv_bins=[];

lll = 0;
for ij=nr_samples_bins

    lll=lll+1; % counter

    wall_dist_temp= reshape(wall_dist(1:ij,:),1,[]);

    omega_tumbling_visc = omega_tumbling_a * scales.viscous_time;
    omega_spinning_visc = omega_spinning_a * scales.viscous_time;

    omega_tumbling_visc(ij+1:end,:)=NaN;
    omega_spinning_visc(ij+1:end,:)=NaN;

    square_omega_tumbling_visc_temp = reshape(omega_tumbling_visc,1,[]).^2;
    square_omega_spinning_visc_temp = reshape(omega_spinning_visc,1,[]).^2;

    [~,mean_tumbl_visc_conv_bins(lll,:),~,~]=Vlad_bin_mean(wall_dist(:),square_omega_tumbling_visc_temp,numBins_z,'log',alpha);
    [~,mean_spinn_visc_conv_bins(lll,:),~,~]=Vlad_bin_mean(wall_dist(:),square_omega_spinning_visc_temp,numBins_z,'log',alpha);

    % how many trajectories have the mean wall distance within each bin
    edges = 10.^linspace(log10(min(wall_dist(:))),log10(max(wall_dist(:))),numBins_z+1);
    binIndices = discretize(mean(wall_dist(1:ij,:),2,'omitnan'), edges);

    for gg=1:numBins_z
        nr_traj_used_conv_bins(lll,gg) = sum(binIndices==gg);
    end

    %clearvars edges binIndices
end

% plot
figure(); hold on; grid on; box on; set(gcf,'Position',[50 50 900 500]);
set(gca,'FontSize',20,'XScale','lin','Yscale','log');
pbaspect([10 5 1])

p0 = plot(nr_traj_used_conv_bins(:,which_bins(1)),mean_tumbl_visc_conv_bins(:,which_bins(1)),'b>-','MarkerSize',5,'MarkerFaceColor','b');
p1 = plot(nr_traj_used_conv_bins(:,which_bins(2)),mean_tumbl_visc_conv_bins(:,which_bins(2)),'bo-','MarkerSize',5,'MarkerFaceColor','b');
p2 = plot(nr_traj_used_conv_bins(:,which_bins(3)),mean_tumbl_visc_conv_bins(:,which_bins(3)),'bs-','MarkerSize',6,'MarkerFaceColor','b');

p3 = plot(nr_traj_used_conv_bins(:,which_bins(1)),mean_spinn_visc_conv_bins(:,which_bins(1)),'r>-','MarkerSize',5,'MarkerFaceColor','r');
p4 = plot(nr_traj_used_conv_bins(:,which_bins(2)),mean_spinn_visc_conv_bins(:,which_bins(2)),'ro-','MarkerSize',5,'MarkerFaceColor','r');
p5 = plot(nr_traj_used_conv_bins(:,which_bins(3)),mean_spinn_visc_conv_bins(:,which_bins(3)),'rs-','MarkerSize',6,'MarkerFaceColor','r');

ylabel('$\langle \omega_i^2 \rangle \tau^2$')
xlabel('Nr. of trajectories')

legend([p0,p1,p2,p3,p4,p5],{['$i=t$, ', '$y^+ = $', num2str(bin_z_Centers_rot_visc(which_bins(1))*1e-3/scales.viscous_length,'%.0f')],...
                            ['$i=t$, ', '$y^+ = $', num2str(bin_z_Centers_rot_visc(which_bins(2))*1e-3/scales.viscous_length,'%.0f')],...
                            ['$i=t$, ', '$y^+ = $', num2str(bin_z_Centers_rot_visc(which_bins(3))*1e-3/scales.viscous_length,'%.0f')],...
                            ['$i=s$, ', '$y^+ = $', num2str(bin_z_Centers_rot_visc(which_bins(1))*1e-3/scales.viscous_length,'%.0f')],...
                            ['$i=s$, ', '$y^+ = $', num2str(bin_z_Centers_rot_visc(which_bins(2))*1e-3/scales.viscous_length,'%.0f')],...
                            ['$i=s$, ', '$y^+ = $', num2str(bin_z_Centers_rot_visc(which_bins(3))*1e-3/scales.viscous_length,'%.0f')]},...
                            'location','eastoutside')


%% tumbling components (2 and 3) scaled by viscous time over y+
% replenish the data
omega_tumbling_2_visc = omega_b_y * scales.viscous_time;
omega_tumbling_3_visc = omega_b_z * scales.viscous_time;
wall_dist = wall_distance_fibre;

omega_spinning_1_visc = omega_spinning_a * scales.viscous_time;

% plot mean square tumbling and spinning over wall-distance
%numBins_z=12;
wall_dist= reshape(wall_dist,1,[]);
square_omega_tumbling_2_visc = reshape(omega_tumbling_2_visc,1,[]).^2;
square_omega_tumbling_3_visc = reshape(omega_tumbling_3_visc,1,[]).^2;

square_omega_spinning_1_visc = reshape(omega_spinning_1_visc,1,[]).^2;

square_omega_body_tbl_spin_visc = square_omega_spinning_1_visc + square_omega_tumbling_2_visc + square_omega_tumbling_3_visc;

mean_squared_tumbling_2_visc = mean(reshape(omega_tumbling_2_visc,1,[]).^2,'omitnan');
mean_squared_tumbling_3_visc = mean(reshape(omega_tumbling_3_visc,1,[]).^2,'omitnan');

[bin_z_Centers_rot_tumbl_2_3_visc,mean_tumbl_2_visc,ci_z_Centers_rot_tumbl_2_3_visc,ci_mean_tumbl_2_visc]=Vlad_bin_mean(wall_dist,square_omega_tumbling_2_visc,numBins_z,'log',alpha);
[bin_z_Centers_rot_tumbl_2_3_visc,mean_tumbl_3_visc,ci_z_Centers_rot_tumbl_2_3_visc,ci_mean_tumbl_3_visc]=Vlad_bin_mean(wall_dist,square_omega_tumbling_3_visc,numBins_z,'log',alpha);

[bin_z_Centers_rot_tumbl_2_3_visc,mean_omega_body_tbl_spin_visc,ci_z_Centers_rot_tumbl_2_3_visc,ci_mean_omega_body_tbl_spin_visc]=Vlad_bin_mean(wall_dist,square_omega_body_tbl_spin_visc,numBins_z,'log',alpha);


% plot
figure(); grid on; box on; hold on;
plot(bin_z_Centers_rot_tumbl_2_3_visc/1e3/scales.viscous_length,mean_tumbl_2_visc,'.-')
plot(bin_z_Centers_rot_tumbl_2_3_visc/1e3/scales.viscous_length,mean_tumbl_3_visc,'.-')
xlabel('$y^+$ [-]')
ylabel('$\langle \omega_i^2 \rangle \tau^2$ [-]')
title('Mean square rotation rates over wall-normal distance')
legend('$\langle \omega_2^2 \rangle \tau$','$\langle \omega_3^2 \rangle \tau$')
set(gca,'Yscale','log','Xscale','log')
xlim([1 1e3])

%% tumbling and spinning components (1, 2 and 3) projected onto x y z scaled by viscous time over y+
% replenish the data
%close all
omega_spinning_1_visc = omega_b_x * scales.viscous_time;
omega_tumbling_2_visc = omega_b_y * scales.viscous_time;
omega_tumbling_3_visc = omega_b_z * scales.viscous_time;
wall_dist = wall_distance_fibre;

proj_x(:,:,1) = ones(size(green,[1 2]));
proj_x(:,:,2) = zeros(size(green,[1 2]));
proj_x(:,:,3) = zeros(size(green,[1 2]));

proj_y(:,:,1) = zeros(size(green,[1 2]));
proj_y(:,:,2) = ones(size(green,[1 2]));
proj_y(:,:,3) = zeros(size(green,[1 2]));

proj_z(:,:,1) = zeros(size(green,[1 2]));
proj_z(:,:,2) = zeros(size(green,[1 2]));
proj_z(:,:,3) = ones(size(green,[1 2]));
% project fitted vectors onto lab reference frame, i.e. (1,0,0),(0,1,0),(0,0,1)
% scaling by the norm of the vectors is not necessary because all vectors
% have the norm=1
proj_spinning_1_x_visc = omega_spinning_1_visc.*sum(red_fitted.*proj_x,3);
proj_spinning_1_y_visc = omega_spinning_1_visc.*sum(red_fitted.*proj_y,3);
proj_spinning_1_z_visc = omega_spinning_1_visc.*sum(red_fitted.*proj_z,3);

proj_tumbling_2_x_visc = omega_tumbling_2_visc.*sum(green_fitted.*proj_x,3);
proj_tumbling_2_y_visc = omega_tumbling_2_visc.*sum(green_fitted.*proj_y,3);
proj_tumbling_2_z_visc = omega_tumbling_2_visc.*sum(green_fitted.*proj_z,3);

proj_tumbling_3_x_visc = omega_tumbling_3_visc.*sum(blue_fitted.*proj_x,3);
proj_tumbling_3_y_visc = omega_tumbling_3_visc.*sum(blue_fitted.*proj_y,3);
proj_tumbling_3_z_visc = omega_tumbling_3_visc.*sum(blue_fitted.*proj_z,3);


% plot mean square tumbling and spinning over wall-distance
%numBins_z=12;
wall_dist= reshape(wall_dist,1,[]);

square_proj_spinning_1_x_visc = reshape(proj_spinning_1_x_visc,1,[]).^2;
square_proj_spinning_1_y_visc = reshape(proj_spinning_1_y_visc,1,[]).^2;
square_proj_spinning_1_z_visc = reshape(proj_spinning_1_z_visc,1,[]).^2;

square_proj_tumbling_2_x_visc = reshape(proj_tumbling_2_x_visc,1,[]).^2;
square_proj_tumbling_2_y_visc = reshape(proj_tumbling_2_y_visc,1,[]).^2;
square_proj_tumbling_2_z_visc = reshape(proj_tumbling_2_z_visc,1,[]).^2;

square_proj_tumbling_3_x_visc = reshape(proj_tumbling_3_x_visc,1,[]).^2;
square_proj_tumbling_3_y_visc = reshape(proj_tumbling_3_y_visc,1,[]).^2;
square_proj_tumbling_3_z_visc = reshape(proj_tumbling_3_z_visc,1,[]).^2;

[bin_z_Centers_rot_proj_spinn_1_visc,mean_proj_spinn_1_x_visc,ci_z_Centers_rot_proj_spinn_1_visc,ci_mean_proj_spinn_1_x_visc]=Vlad_bin_mean(wall_dist,square_proj_spinning_1_x_visc,numBins_z,'log',alpha);
[bin_z_Centers_rot_proj_spinn_1_visc,mean_proj_spinn_1_y_visc,ci_z_Centers_rot_proj_spinn_1_visc,ci_mean_proj_spinn_1_y_visc]=Vlad_bin_mean(wall_dist,square_proj_spinning_1_y_visc,numBins_z,'log',alpha);
[bin_z_Centers_rot_proj_spinn_1_visc,mean_proj_spinn_1_z_visc,ci_z_Centers_rot_proj_spinn_1_visc,ci_mean_proj_spinn_1_z_visc]=Vlad_bin_mean(wall_dist,square_proj_spinning_1_z_visc,numBins_z,'log',alpha);

[bin_z_Centers_rot_proj_tumbl_2_3_visc,mean_proj_tumbl_2_x_visc,ci_z_Centers_rot_proj_tumbl_2_3_visc,ci_mean_proj_tumbl_2_x_visc]=Vlad_bin_mean(wall_dist,square_proj_tumbling_2_x_visc,numBins_z,'log',alpha);
[bin_z_Centers_rot_proj_tumbl_2_3_visc,mean_proj_tumbl_2_y_visc,ci_z_Centers_rot_proj_tumbl_2_3_visc,ci_mean_proj_tumbl_2_y_visc]=Vlad_bin_mean(wall_dist,square_proj_tumbling_2_y_visc,numBins_z,'log',alpha);
[bin_z_Centers_rot_proj_tumbl_2_3_visc,mean_proj_tumbl_2_z_visc,ci_z_Centers_rot_proj_tumbl_2_3_visc,ci_mean_proj_tumbl_2_z_visc]=Vlad_bin_mean(wall_dist,square_proj_tumbling_2_z_visc,numBins_z,'log',alpha);

[bin_z_Centers_rot_proj_tumbl_2_3_visc,mean_proj_tumbl_3_x_visc,ci_z_Centers_rot_proj_tumbl_2_3_visc,ci_mean_proj_tumbl_3_x_visc]=Vlad_bin_mean(wall_dist,square_proj_tumbling_3_x_visc,numBins_z,'log',alpha);
[bin_z_Centers_rot_proj_tumbl_2_3_visc,mean_proj_tumbl_3_y_visc,ci_z_Centers_rot_proj_tumbl_2_3_visc,ci_mean_proj_tumbl_3_y_visc]=Vlad_bin_mean(wall_dist,square_proj_tumbling_3_y_visc,numBins_z,'log',alpha);
[bin_z_Centers_rot_proj_tumbl_2_3_visc,mean_proj_tumbl_3_z_visc,ci_z_Centers_rot_proj_tumbl_2_3_visc,ci_mean_proj_tumbl_3_z_visc]=Vlad_bin_mean(wall_dist,square_proj_tumbling_3_z_visc,numBins_z,'log',alpha);

% plot
figure(); grid on; box on; hold on;

plot(bin_z_Centers_rot_proj_tumbl_2_3_visc/1e3/scales.viscous_length,mean_proj_spinn_1_x_visc,'.-')
plot(bin_z_Centers_rot_proj_tumbl_2_3_visc/1e3/scales.viscous_length,mean_proj_spinn_1_y_visc,'.-')
plot(bin_z_Centers_rot_proj_tumbl_2_3_visc/1e3/scales.viscous_length,mean_proj_spinn_1_z_visc,'.-')

plot(bin_z_Centers_rot_proj_tumbl_2_3_visc/1e3/scales.viscous_length,mean_proj_tumbl_2_x_visc,'.-')
plot(bin_z_Centers_rot_proj_tumbl_2_3_visc/1e3/scales.viscous_length,mean_proj_tumbl_2_y_visc,'.-')
plot(bin_z_Centers_rot_proj_tumbl_2_3_visc/1e3/scales.viscous_length,mean_proj_tumbl_2_z_visc,'.-')

plot(bin_z_Centers_rot_proj_tumbl_2_3_visc/1e3/scales.viscous_length,mean_proj_tumbl_3_x_visc,'.-')
plot(bin_z_Centers_rot_proj_tumbl_2_3_visc/1e3/scales.viscous_length,mean_proj_tumbl_3_y_visc,'.-')
plot(bin_z_Centers_rot_proj_tumbl_2_3_visc/1e3/scales.viscous_length,mean_proj_tumbl_3_z_visc,'.-')

xlabel('$y^+$ [-]')
ylabel('$\langle \omega_i^2 \rangle \tau^2$ [-]')
title('Mean square rotation rates over wall-normal distance')
legend('$\langle \omega_2^2\  \vec{e_2} \cdot \vec{x} \rangle \tau^2$',...
    '$\langle \omega_2^2\  \vec{e_2} \cdot \vec{y} \rangle \tau^2$',...
    '$\langle \omega_2^2\  \vec{e_2} \cdot \vec{z} \rangle \tau^2$', ...
    '$\langle \omega_3^2\  \vec{e_3} \cdot \vec{x} \rangle \tau^2$',...
    '$\langle \omega_3^2\  \vec{e_3} \cdot \vec{y} \rangle \tau^2$',...
    '$\langle \omega_3^2\  \vec{e_3} \cdot \vec{z} \rangle \tau^2$','location','best')
set(gca,'Yscale','log','Xscale','log','FontSize',15)
set(gca,'FontSize',15)
%xlim([1 1e3])

%% rotation rate components (x,y,z) scaled by viscous time over y+
% replenish the data
omega_x_visc = omega_s_x * scales.viscous_time;
omega_y_visc = omega_s_y * scales.viscous_time;
omega_z_visc = omega_s_z * scales.viscous_time;
wall_dist = wall_distance_fibre;

% plot mean square tumbling and spinning over wall-distance
%numBins_z=12;
wall_dist= reshape(wall_dist,1,[]);
square_omega_x_visc = reshape(omega_x_visc,1,[]).^2;
square_omega_y_visc = reshape(omega_y_visc,1,[]).^2;
square_omega_z_visc = reshape(omega_z_visc,1,[]).^2;

% solid body rotation
square_omega_body_visc = square_omega_x_visc + square_omega_y_visc + square_omega_z_visc;

mean_squared_x_visc = mean(reshape(omega_x_visc,1,[]).^2,'omitnan');
mean_squared_y_visc = mean(reshape(omega_y_visc,1,[]).^2,'omitnan');
mean_squared_z_visc = mean(reshape(omega_z_visc,1,[]).^2,'omitnan');

[bin_z_Centers_rot_x_y_z_visc,mean_omega_x_visc,ci_z_Centers_rot_x_y_z_visc,ci_mean_omega_x_visc]=Vlad_bin_mean(wall_dist,square_omega_x_visc,numBins_z,'log',alpha);
[bin_z_Centers_rot_x_y_z_visc,mean_omega_y_visc,ci_z_Centers_rot_x_y_z_visc,ci_mean_omega_y_visc]=Vlad_bin_mean(wall_dist,square_omega_y_visc,numBins_z,'log',alpha);
[bin_z_Centers_rot_x_y_z_visc,mean_omega_z_visc,ci_z_Centers_rot_x_y_z_visc,ci_mean_omega_z_visc]=Vlad_bin_mean(wall_dist,square_omega_z_visc,numBins_z,'log',alpha);

[bin_z_Centers_rot_x_y_z_visc,mean_omega_body_visc,ci_z_Centers_rot_x_y_z_visc,ci_mean_omega_body_visc]=Vlad_bin_mean(wall_dist,square_omega_body_visc,numBins_z,'log',alpha);


% plot
figure(); grid on; box on; hold on;
plot(bin_z_Centers_rot_x_y_z_visc/1e3/scales.viscous_length,mean_omega_x_visc,'.-')
plot(bin_z_Centers_rot_x_y_z_visc/1e3/scales.viscous_length,mean_omega_y_visc,'.-')
plot(bin_z_Centers_rot_x_y_z_visc/1e3/scales.viscous_length,mean_omega_z_visc,'.-')
xlabel('$y^+$ [-]')
ylabel('$\langle \omega_i^2 \rangle \tau^2$ [-]')
title('Mean square rotation rates over wall-normal distance')
legend('$\langle \omega_x^2 \rangle \tau$','$\langle \omega_y^2 \rangle \tau$','$\langle \omega_z^2 \rangle \tau$')
set(gca,'Yscale','log','Xscale','log')
xlim([1 1e3])

%% tumbling and spinning scaled by viscous time over l/eta
% replenish the data
omega_tumbling_visc = omega_tumbling_a * scales.viscous_time;
omega_spinning_visc = omega_spinning_a * scales.viscous_time;

l_eta = l_nominal *1e-3 * local_l_kol_fibre.^(-1);
l_eta = repmat(l_eta,1,size(omega_tumbling_kol,2));

l_eta(isnan(omega_tumbling_visc))=NaN;

wall_dist = wall_distance_fibre;
wall_dist= reshape(wall_dist,1,[]);

% plot mean square tumbling and spinning over wall-distance
numBins=13;
l_eta= reshape(l_eta,1,[]);
square_omega_tumbling_visc = reshape(omega_tumbling_visc,1,[]).^2;
square_omega_spinning_visc = reshape(omega_spinning_visc,1,[]).^2;

[bin_l_eta_Centers_visc,mean_tumbl_l_eta_visc,ci_l_eta_Centers_visc,ci_mean_tumbl_l_eta_visc]=Vlad_bin_mean(l_eta,square_omega_tumbling_visc,numBins,'lin',alpha);
[bin_l_eta_Centers_visc,mean_spinn_l_eta_visc,ci_l_eta_Centers_visc,ci_mean_spinn_l_eta_visc]=Vlad_bin_mean(l_eta,square_omega_spinning_visc,numBins,'lin',alpha);

% plot
figure(); grid on; box on; hold on;
plot(bin_l_eta_Centers_visc,mean_tumbl_l_eta_visc,'.-')
plot(bin_l_eta_Centers_visc,mean_spinn_l_eta_visc,'.-')
xlabel('$l/\eta$ [-]')
ylabel('$\langle \omega_i^2 \rangle \tau^2$ [-]')
title('Mean square rotation rates over l over eta')
legend('$\langle \omega_t^2 \rangle \tau_{\eta}$','$\langle \omega_s^2 \rangle \tau_{\eta}$')
set(gca,'Yscale','lin','Xscale','lin')

%% tumbling and spinning scaled by Kolmogorov time over l/eta
% replenish the data
omega_tumbling_kol = omega_tumbling_a .* local_t_kol_fibre;
omega_spinning_kol = omega_spinning_a .* local_t_kol_fibre;
l_eta = l_nominal *1e-3 * local_l_kol_fibre.^(-1);
%l_eta = avg_fibre_length*1e-3 .* local_l_kol_fibre.^(-1);
l_eta = repmat(l_eta,1,size(omega_tumbling_kol,2));

l_eta(isnan(omega_tumbling_kol))=NaN;

% plot mean square tumbling and spinning over wall-distance
numBins=10;
square_omega_tumbling_kol = reshape(omega_tumbling_kol,1,[]).^2;
square_omega_spinning_kol = reshape(omega_spinning_kol,1,[]).^2;

[~, edges] = histcounts(l_eta, numBins);
binIndices = discretize(l_eta, edges);

binned_tumbl_l_eta = cell(1, numBins);
binned_spinn_l_eta = cell(1, numBins);
mean_tumbl_l_eta=[];
mean_spinn_l_eta=[];


for i = 1:numBins
    binned_tumbl_l_eta{i} = square_omega_tumbling_kol(binIndices == i);
    binned_spinn_l_eta{i} = square_omega_spinning_kol(binIndices == i);

    mean_tumbl_l_eta(i)=mean(binned_tumbl_l_eta{i}, 'omitnan');
    mean_spinn_l_eta(i)=mean(binned_spinn_l_eta{i}, 'omitnan');
end

bin_l_eta_Centers = (edges(1:end-1) + edges(2:end)) / 2;

% plot
figure('Position',[100 100 800 700]); grid on; box on; hold on;
plot(bin_l_eta_Centers,mean_tumbl_l_eta,'.-')
plot(bin_l_eta_Centers,mean_spinn_l_eta,'.-')
xlabel('$l_f/{\eta}$ [-]')
ylabel('$\langle \omega_i^2 \rangle \tau_{\eta}^2$ [-]')
title('Mean square rotation rates over dimensionless fiber length')

set(gca,'Yscale','log','Xscale','log','FontSize',20)

l_eta_scaling = [6,13];
scaling = l_eta_scaling.^(-4/3);
plot(l_eta_scaling,scaling,'k-')

legend('$\langle \omega_t^2 \rangle \tau^2_{\eta}$',...
    '$\langle \omega_s^2 \rangle \tau^2_{\eta}$',...
    '$\langle \omega^2 \rangle \tau^2_{\eta}\propto(l/\eta)^{-4/3}$ Parsa\& Voth (2014)','FontSize',20,'location','best')

%% tumbling and spinning scaled by Kolmogorov over k^*
% replenish the data
omega_tumbling_kol = omega_tumbling_a .* local_t_kol_fibre;
omega_spinning_kol = omega_spinning_a .* local_t_kol_fibre;
k_star = repmat(avg_fibre_curvature,1,size(omega_tumbling_kol,2));

k_star(isnan(omega_tumbling_kol))=NaN;

% plot mean square tumbling and spinning over wall-distance
numBins=5;
square_omega_tumbling_kol = reshape(omega_tumbling_kol,1,[]).^2;
square_omega_spinning_kol = reshape(omega_spinning_kol,1,[]).^2;

[~, edges] = histcounts(k_star, numBins);
binIndices = discretize(k_star, edges);

binned_tumbl_k_star = cell(1, numBins);
binned_spinn_k_star = cell(1, numBins);
mean_tumbl_k_star=[];
mean_spinn_k_star=[];

for i = 1:numBins
    binned_tumbl_k_star{i} = square_omega_tumbling_kol(binIndices == i);
    binned_spinn_k_star{i} = square_omega_spinning_kol(binIndices == i);

    mean_tumbl_k_star(i)=mean(binned_tumbl_k_star{i}, 'omitnan');
    mean_spinn_k_star(i)=mean(binned_spinn_k_star{i}, 'omitnan');
end

bin_k_star_Centers = (edges(1:end-1) + edges(2:end)) / 2;

% plot
figure(); grid on; box on; hold on;
plot(bin_k_star_Centers,mean_tumbl_k_star,'.-')
plot(bin_k_star_Centers,mean_spinn_k_star,'.-')
xlabel('$\kappa^*$ [-]')
ylabel('$\langle \omega_i^2 \rangle \tau_{\eta}^2$ [-]')
title('Mean square rotation rates over dimensionless fiber curvature')
legend('$\langle \omega_t^2 \rangle \tau_{\eta}$','$\langle \omega_s^2 \rangle \tau_{\eta}$')
set(gca,'Yscale','log','Xscale','lin')

%% tumbling and spinning scaled by viscous time scale over k^*
% replenish the data
omega_tumbling_kol = omega_tumbling_a .* scales.viscous_time;
omega_spinning_kol = omega_spinning_a .* scales.viscous_time;
k_star = repmat(avg_fibre_curvature,1,size(omega_tumbling_kol,2));

k_star(isnan(omega_tumbling_kol))=NaN;

% plot mean square tumbling and spinning over wall-distance
numBins=3;
k_star = reshape(k_star,1,[]);
square_omega_tumbling_kol = reshape(omega_tumbling_kol,1,[]).^2;
square_omega_spinning_kol = reshape(omega_spinning_kol,1,[]).^2;

[bin_k_star_Centers_visc,mean_tumbl_k_star_visc,ci_k_star_Centers_visc,ci_mean_tumbl_k_star_visc]=Vlad_bin_mean(k_star,square_omega_tumbling_kol,numBins,'lin',alpha);
[bin_k_star_Centers_visc,mean_spinn_k_star_visc,ci_k_star_Centers_visc,ci_mean_spinn_k_star_visc]=Vlad_bin_mean(k_star,square_omega_spinning_kol,numBins,'lin',alpha);

% plot
figure();grid on; box on; hold on;
%plot(bin_k_star_Centers_visc,mean_tumbl_k_star_visc,'.-')
errorbar(bin_k_star_Centers_visc,mean_tumbl_k_star_visc,...
        mean_tumbl_k_star_visc-ci_mean_tumbl_k_star_visc(:,1),...
        mean_tumbl_k_star_visc-ci_mean_tumbl_k_star_visc(:,1),...
        bin_k_star_Centers_visc' - ci_k_star_Centers_visc(:,1),...
        bin_k_star_Centers_visc' - ci_k_star_Centers_visc(:,1));
%plot(bin_k_star_Centers_visc,mean_spinn_k_star_visc,'.-')
errorbar(bin_k_star_Centers_visc,mean_spinn_k_star_visc,...
        mean_spinn_k_star_visc-ci_mean_spinn_k_star_visc(:,1),...
        mean_spinn_k_star_visc-ci_mean_spinn_k_star_visc(:,1),...
        bin_k_star_Centers_visc' - ci_k_star_Centers_visc(:,1),...
        bin_k_star_Centers_visc' - ci_k_star_Centers_visc(:,1));

xlabel('$\kappa^*$ [-]')
ylabel('$\langle \omega_i^2 \rangle \tau^2$ [-]')
title('Mean square rotation rates over dimensionless fiber curvature')
legend('$\langle \omega_t^2 \rangle \tau_{\eta}$','$\langle \omega_s^2 \rangle \tau_{\eta}$','location','best')
set(gca,'Yscale','log','Xscale','lin')

%% check convergence
figure(); hold on; grid on; box on; set(gcf,'Position',[1 50 1800 900])
nr_samples = round(logspace(log10(10),log10(size(omega_spinning,1)),50));
% replenish the data
omega_tumbling_kol = omega_tumbling_a .* local_t_kol_fibre;
omega_spinning_kol = omega_spinning_a .* local_t_kol_fibre;
for ij=nr_samples
    omega_tumbling_kol_conv = omega_tumbling_kol(1:ij,:);
    omega_spinning_kol_conv = omega_spinning_kol(1:ij,:);

    mean_tumbling_kol_conv = mean(reshape(omega_tumbling_kol_conv,1,[]),'omitnan');
    mean_spinning_kol_conv = mean(reshape(omega_spinning_kol_conv,1,[]),'omitnan');

    subplot(1,3,1); hold on; grid on; box on; set(gca,'Fontsize',20,'XScale','log')
    p1=plot(ij,mean_tumbling_kol_conv,'b.','MarkerSize',20);
    p2=plot(ij,mean_spinning_kol_conv,'r.','MarkerSize',20);
    title('Convergence mean')
    legend([p1,p2],'$\langle\Omega_t\rangle\tau_{\eta}$','$\langle\Omega_s\rangle\tau_{\eta}$')
    xlabel('Nr. Tracks')

    mean_squared_tumbling_kol_conv = mean(reshape(omega_tumbling_kol_conv,1,[]).^2,'omitnan');
    mean_squared_spinning_kol_conv = mean(reshape(omega_spinning_kol_conv,1,[]).^2,'omitnan');

    subplot(1,3,2); hold on; grid on; box on; set(gca,'Fontsize',20,'XScale','log','YScale','lin')
    p1=plot(ij,mean_squared_tumbling_kol_conv,'b.','MarkerSize',20);
    p2=plot(ij,mean_squared_spinning_kol_conv,'r.','MarkerSize',20);
    title('Convergence mean squared')
    legend([p1,p2],'$\langle\Omega_t^2\rangle \tau_{\eta}^2$','$\langle\Omega_s^2\rangle\tau_{\eta}^2$')
    xlabel('Nr. Tracks')

    var_tumbling_kol_conv=mean(reshape(omega_tumbling_kol_conv-mean_tumbling_kol_conv,1,[]).^2,'omitnan');
    var_spinning_kol_conv=mean(reshape(omega_spinning_kol_conv-mean_spinning_kol_conv,1,[]).^2,'omitnan');

    subplot(1,3,3); hold on; grid on; box on; set(gca,'Fontsize',20,'XScale','log','YScale','lin')
    p1=plot(ij,var_tumbling_kol_conv,'b.','MarkerSize',20);
    p2=plot(ij,var_spinning_kol_conv,'r.','MarkerSize',20);
    title('Convergence variance')
    legend([p1,p2],'$\langle(\Omega_t - \langle\Omega_t\rangle)^2\rangle\tau_{\eta}^2$','$\langle(\Omega_s - \langle\Omega_s\rangle)^2\rangle\tau_{\eta}^2$')
    xlabel('Nr. Tracks')
end


%% check convergence for different curvature classes
% replenish the data
omega_tumbling_kol = omega_tumbling_a .* local_t_kol_fibre;
omega_spinning_kol = omega_spinning_a .* local_t_kol_fibre;

figure(); hold on; grid on; box on; set(gcf,'Position',[1 50 1800 900])

% std(k*)<0.05
% low curvature k*<0.28
omega_tumbling_kol(avg_fibre_curvature>0.28,:)=NaN;
omega_spinning_kol(avg_fibre_curvature>0.28,:)=NaN;

nr_samples = round(logspace(log10(10),log10(size(omega_spinning_kol,1)),50));
for ij=nr_samples
    omega_tumbling_kol_conv = omega_tumbling_kol(1:ij,:);
    omega_spinning_kol_conv = omega_spinning_kol(1:ij,:);

    mean_tumbling_kol_conv = mean(reshape(omega_tumbling_kol_conv,1,[]),'omitnan');
    mean_spinning_kol_conv = mean(reshape(omega_spinning_kol_conv,1,[]),'omitnan');

    subplot(1,3,1); hold on; grid on; box on;
    plot(ij,mean_tumbling_kol_conv,'bo','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    plot(ij,mean_spinning_kol_conv,'ro','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);

    mean_squared_tumbling_kol_conv = mean(reshape(omega_tumbling_kol_conv,1,[]).^2,'omitnan');
    mean_squared_spinning_kol_conv = mean(reshape(omega_spinning_kol_conv,1,[]).^2,'omitnan');

    subplot(1,3,2); hold on; grid on; box on;
    plot(ij,mean_squared_tumbling_kol_conv,'bo','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    plot(ij,mean_squared_spinning_kol_conv,'ro','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);

    var_tumbling_kol_conv=mean(reshape(omega_tumbling_kol_conv-mean_tumbling_kol_conv,1,[]).^2,'omitnan');
    var_spinning_kol_conv=mean(reshape(omega_spinning_kol_conv-mean_spinning_kol_conv,1,[]).^2,'omitnan');

    subplot(1,3,3); hold on; grid on; box on;
    plot(ij,var_tumbling_kol_conv,'bo','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    plot(ij,var_spinning_kol_conv,'ro','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);


end

% replenish the data
omega_tumbling_kol = omega_tumbling_a .* local_t_kol_fibre;
omega_spinning_kol = omega_spinning_a .* local_t_kol_fibre;

% std(k*)<0.05
% low curvature 0.28<k*<0.42
omega_tumbling_kol(or(avg_fibre_curvature>0.42,avg_fibre_curvature<0.28),:)=NaN;
omega_spinning_kol(or(avg_fibre_curvature>0.42,avg_fibre_curvature<0.28),:)=NaN;

nr_samples = round(logspace(log10(10),log10(size(omega_spinning_kol,1)),50));
for ij=nr_samples
    omega_tumbling_kol_conv = omega_tumbling_kol(1:ij,:);
    omega_spinning_kol_conv = omega_spinning_kol(1:ij,:);

    mean_tumbling_kol_conv = mean(reshape(omega_tumbling_kol_conv,1,[]),'omitnan');
    mean_spinning_kol_conv = mean(reshape(omega_spinning_kol_conv,1,[]),'omitnan');

    subplot(1,3,1); hold on;
    plot(ij,mean_tumbling_kol_conv,'bs','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    plot(ij,mean_spinning_kol_conv,'rs','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);

    mean_squared_tumbling_kol_conv = mean(reshape(omega_tumbling_kol_conv,1,[]).^2,'omitnan');
    mean_squared_spinning_kol_conv = mean(reshape(omega_spinning_kol_conv,1,[]).^2,'omitnan');

    subplot(1,3,2); hold on;
    plot(ij,mean_squared_tumbling_kol_conv,'bs','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    plot(ij,mean_squared_spinning_kol_conv,'rs','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);

    var_tumbling_kol_conv=mean(reshape(omega_tumbling_kol_conv-mean_tumbling_kol_conv,1,[]).^2,'omitnan');
    var_spinning_kol_conv=mean(reshape(omega_spinning_kol_conv-mean_spinning_kol_conv,1,[]).^2,'omitnan');

    subplot(1,3,3); hold on;
    plot(ij,var_tumbling_kol_conv,'bs','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    plot(ij,var_spinning_kol_conv,'rs','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);

end

% replenish the data
omega_tumbling_kol = omega_tumbling_a .* local_t_kol_fibre;
omega_spinning_kol = omega_spinning_a .* local_t_kol_fibre;

% std(k*)<0.05
% low curvature k*>0.42
omega_tumbling_kol(avg_fibre_curvature<0.42,:)=NaN;
omega_spinning_kol(avg_fibre_curvature<0.42,:)=NaN;

nr_samples = round(logspace(log10(10),log10(size(omega_spinning_kol,1)),50));
for ij=nr_samples
    omega_tumbling_kol_conv = omega_tumbling_kol(1:ij,:);
    omega_spinning_kol_conv = omega_spinning_kol(1:ij,:);

    mean_tumbling_kol_conv = mean(reshape(omega_tumbling_kol_conv,1,[]),'omitnan');
    mean_spinning_kol_conv = mean(reshape(omega_spinning_kol_conv,1,[]),'omitnan');

    subplot(1,3,1); hold on; grid on; box on; set(gca,'Fontsize',20,'XScale','log')
    plot(ij,mean_tumbling_kol_conv,'b^','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    plot(ij,mean_spinning_kol_conv,'r^','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    title('Convergence mean')
    %legend([p1,p2],'$\langle\Omega_t\rangle\tau_{\eta}$','$\langle\Omega_s\rangle\tau_{\eta}$')
    xlabel('Nr. Tracks')
    ylabel('$\langle\Omega_t\rangle\tau_{\eta}$, $\langle\Omega_s\rangle\tau_{\eta}$')


    xlim([10 1e4])
    ylim([0 0.6])


    mean_squared_tumbling_kol_conv = mean(reshape(omega_tumbling_kol_conv,1,[]).^2,'omitnan');
    mean_squared_spinning_kol_conv = mean(reshape(omega_spinning_kol_conv,1,[]).^2,'omitnan');

    subplot(1,3,2); hold on; grid on; box on; set(gca,'Fontsize',20,'XScale','log','YScale','log')
    plot(ij,mean_squared_tumbling_kol_conv,'b^','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    plot(ij,mean_squared_spinning_kol_conv,'r^','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    % voth soldati 2017, AR=40
    plot([10 1e4],[0.1 0.1],'b-','LineWidth',2) % tumbling
    plot([10 1e4],[0.17 0.17],'r-','LineWidth',2) % spinning

    title('Convergence mean squared')
    %legend([p1,p2],'$\langle\Omega_t^2\rangle \tau_{\eta}^2$','$\langle\Omega_s^2\rangle\tau_{\eta}^2$')
    xlabel('Nr. Tracks')
    ylabel('$\langle\Omega_t^2\rangle \tau_{\eta}^2$, $\langle\Omega_s^2\rangle\tau_{\eta}^2$')
    xlim([10 1e4])
    ylim([1e-2 2e0])

    var_tumbling_kol_conv=mean(reshape(omega_tumbling_kol_conv-mean_tumbling_kol_conv,1,[]).^2,'omitnan');
    var_spinning_kol_conv=mean(reshape(omega_spinning_kol_conv-mean_spinning_kol_conv,1,[]).^2,'omitnan');

    subplot(1,3,3); hold on; grid on; box on; set(gca,'Fontsize',20,'XScale','log','YScale','log')
    plot(ij,var_tumbling_kol_conv,'b^','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    plot(ij,var_spinning_kol_conv,'r^','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    title('Convergence variance')
    %legend([p1,p2],'$\langle(\Omega_t - \langle\Omega_t\rangle)^2\rangle\tau_{\eta}^2$','$\langle(\Omega_s - \langle\Omega_s\rangle)^2\rangle\tau_{\eta}^2$')
    xlabel('Nr. Tracks')
    ylabel('$\langle(\Omega_t - \langle\Omega_t\rangle)^2\rangle\tau_{\eta}^2$, $\langle(\Omega_s - \langle\Omega_s\rangle)^2\rangle\tau_{\eta}^2$')
    xlim([10 1e4])
    ylim([1e-2 1e0])

        % make legend
    p1=plot(1,100,'bo','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    p2=plot(1,100,'bs','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    p3=plot(1,100,'b^','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    p4=plot(1,100,'ro','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    p5=plot(1,100,'rs','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);
    p6=plot(1,100,'r^','MarkerSize',7,'MarkerFaceColor','none','LineWidth',2);

    p7=plot(1,100,'b-','LineWidth',2);
    p8=plot(1,100,'r-','LineWidth',2);

    legend([p1,p2,p3,p4,p5,p6,p7,p8],'tumbl. $\kappa^* < 0.28$','tumbl. $0.28 <\kappa^* < 0.42$', 'tumbl. $\kappa^* > 0.42$', ...
                            'spin. $\kappa^* < 0.28$','spin. $0.28 <\kappa^* < 0.42$', 'spin. $\kappa^* > 0.42$',...
                            'tumbl. Voth \& Soldati 2017','spin. Voth \& Soldati 2017','location','southoutside')
end


%% PDF of spinning and tumbling scaled by the mean squared - for different curvature classes
figure(); hold on; grid on; box on;
set(gcf,'Position',[50 50 800 800])
set(gca,'FontSize',20,'YScale','log');
ylabel('p.d.f')
xlabel('$\omega_s^2 / \langle \omega_s^2 \rangle$, $\omega_t^2 / \langle \omega_t^2 \rangle$')

% replenish the data
omega_tumbling_kol = omega_tumbling_a;
omega_spinning_kol = omega_spinning_a;

%%%% std(k*)<0.05
% low curvature k*<0.28
omega_tumbling_kol(avg_fibre_curvature>0.28,:)=NaN;
omega_spinning_kol(avg_fibre_curvature>0.28,:)=NaN;

mean_squared_tumbling_kol = mean(reshape(omega_tumbling_kol,1,[]).^2,'omitnan');
mean_squared_spinning_kol = mean(reshape(omega_spinning_kol,1,[]).^2,'omitnan');

% spinning pdf
nbins=20;
data = reshape(omega_spinning_kol,1,[]).^2 / mean_squared_spinning_kol;
data(isnan(data))=[];
[n_spinning_low_k,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_spinning_low_k=(edges(1:end-1)+edges(2:end))/2;
plot(centers_spinning_low_k,n_spinning_low_k,'ro-','MarkerSize',15)

% tumbling pdf
nbins=20;
data = reshape(omega_tumbling_kol,1,[]).^2 / mean_squared_tumbling_kol;
data(isnan(data))=[];
[n_tumbling_low_k,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_tumbling_low_k=(edges(1:end-1)+edges(2:end))/2;
plot(centers_tumbling_low_k,n_tumbling_low_k,'bo-','MarkerSize',15)


% replenish the data
omega_tumbling_kol = omega_tumbling_a;
omega_spinning_kol = omega_spinning_a;

%%%% std(k*)<0.05
% medium curvature 0.28<k*<0.42
omega_tumbling_kol(or(avg_fibre_curvature>0.42,avg_fibre_curvature<0.28),:)=NaN;
omega_spinning_kol(or(avg_fibre_curvature>0.42,avg_fibre_curvature<0.28),:)=NaN;

% spinning pdf
nbins=20;
data = reshape(omega_spinning_kol,1,[]).^2 / mean_squared_spinning_kol;
data(isnan(data))=[];
[n_spinning_med_k,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_spinning_med_k=(edges(1:end-1)+edges(2:end))/2;
plot(centers_spinning_med_k,n_spinning_med_k,'rs-','MarkerSize',15)

% tumbling pdf
nbins=23;
data = reshape(omega_tumbling_kol,1,[]).^2 / mean_squared_tumbling_kol;
data(isnan(data))=[];
[n_tumbling_med_k,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_tumbling_med_k=(edges(1:end-1)+edges(2:end))/2;
plot(centers_tumbling_med_k,n_tumbling_med_k,'bs-','MarkerSize',15)


% replenish the data
omega_tumbling_kol = omega_tumbling_a;
omega_spinning_kol = omega_spinning_a;

%%%% std(k*)<0.05
% low curvature k*>0.42
omega_tumbling_kol(avg_fibre_curvature<0.42,:)=[];
omega_spinning_kol(avg_fibre_curvature<0.42,:)=[];

% spinning pdf
nbins=20;
data = reshape(omega_spinning_kol,1,[]).^2 / mean_squared_spinning_kol;
data(isnan(data))=[];
[n_spinning_high_k,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_spinning_high_k=(edges(1:end-1)+edges(2:end))/2;
plot(centers_spinning_high_k,n_spinning_high_k,'rd-','MarkerSize',15)

% tumbling pdf
nbins=22;
data = reshape(omega_tumbling_kol,1,[]).^2 / mean_squared_tumbling_kol;
data(isnan(data))=[];
[n_tumbling_high_k,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_tumbling_high_k=(edges(1:end-1)+edges(2:end))/2;
plot(centers_tumbling_high_k,n_tumbling_high_k,'bd-','MarkerSize',15)

%% PDF of spinning and tumbling scaled by the mean squared - all curvatures together
figure(); hold on; grid on; box on;
set(gcf,'Position',[50 50 800 800])
set(gca,'FontSize',20,'YScale','log');
ylabel('p.d.f')
xlabel('$\omega_s^2 / \langle \omega_s^2 \rangle$, $\omega_t^2 / \langle \omega_t^2 \rangle$')

% replenish the data
omega_tumbling_kol = omega_tumbling_a ;
omega_spinning_kol = omega_spinning_a ;


mean_squared_tumbling_kol = mean(reshape(omega_tumbling_kol,1,[]).^2,'omitnan');
mean_squared_spinning_kol = mean(reshape(omega_spinning_kol,1,[]).^2,'omitnan');

% spinning pdf
nbins=15;
data = reshape(omega_spinning_kol,1,[]).^2 / mean_squared_spinning_kol;
data(isnan(data))=[];

[n_spinning,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_spinning=(edges(1:end-1)+edges(2:end))/2;
plot(centers_spinning,n_spinning,'r.-','MarkerSize',15)

data=[];

% tumbling pdf
nbins=50;
data = reshape(omega_tumbling_kol,1,[]).^2 / mean_squared_tumbling_kol;
data(isnan(data))=[];

[n_tumbling,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_tumbling=(edges(1:end-1)+edges(2:end))/2;
plot(centers_tumbling,n_tumbling,'b.-','MarkerSize',15)

% remove half the data
numToRemove = round(length(data) / 2);
% Generate random indices to remove
indicesToRemove = randperm(length(data), numToRemove);
% Create a new vector without the randomly selected elements
newData = data;
newData(indicesToRemove) = [];
[n_tumbling_half,~]=histcounts(newData,edges,"Normalization","probability");

% remove 3/4th of the data
numToRemove = round(3*length(data) / 4);
% Generate random indices to remove
indicesToRemove = randperm(length(data), numToRemove);
% Create a new vector without the randomly selected elements
newData = data;
newData(indicesToRemove) = [];
[n_tumbling_quarter,~]=histcounts(newData,edges,"Normalization","probability");

% remove 90% of the data
numToRemove = round(9*length(data) / 10);
% Generate random indices to remove
indicesToRemove = randperm(length(data), numToRemove);
% Create a new vector without the randomly selected elements
newData = data;
newData(indicesToRemove) = [];
[n_tumbling_one_tenth,~]=histcounts(newData,edges,"Normalization","probability");

% remove 99% of the data
numToRemove = round(99*length(data) / 100);
% Generate random indices to remove
indicesToRemove = randperm(length(data), numToRemove);
% Create a new vector without the randomly selected elements
newData = data;
newData(indicesToRemove) = [];
[n_tumbling_one_onehundreth,~]=histcounts(newData,edges,"Normalization","probability");

ylim([1e-7 1])

legend('$\omega_s^2 / \langle \omega_s^2 \rangle$', '$\omega_t^2 / \langle \omega_t^2 \rangle$')

%% PDF spinning and tumbling to compare with Oehmke et al. 2021
figure(); hold on; grid on; box on;
set(gcf,'Position',[50 50 800 800])
set(gca,'FontSize',20,'YScale','log');
ylabel('p.d.f')
xlabel('$(\omega_s - \langle \omega_s \rangle )/$rms$(\omega_s)$','Interpreter','latex')

omega_spinning_oehmke = reshape(deg2rad(omega_spinning)/p.dt,1,[]);

rms_spinning_oehmke = rms(omega_spinning_oehmke,'omitnan');
mean_spinning_oehmke = mean(omega_spinning_oehmke,'omitnan');

data = (omega_spinning_oehmke - mean_spinning_oehmke)/rms_spinning_oehmke;

[centers_spin_oehmke,count_spin_oehmke]=compute_pdf(data,'BinNumber',21,'pdf');
plot(centers_spin_oehmke,count_spin_oehmke,'b.-','MarkerSize',15)

figure(); hold on; grid on; box on;
set(gcf,'Position',[50 50 800 800])
set(gca,'FontSize',20,'YScale','log');
ylabel('p.d.f')
xlabel('$(\omega_s - \langle \omega_s \rangle )/$rms$(\omega_s)$','Interpreter','latex')

omega_tumbling_oehmke = reshape(deg2rad(omega_tumbling)/p.dt,1,[]);

rms_tumbling_oehmke = rms(omega_tumbling_oehmke,'omitnan');
mean_tumbling_oehmke = mean(omega_tumbling_oehmke,'omitnan');

data = (omega_tumbling_oehmke - mean_tumbling_oehmke)/rms_tumbling_oehmke;

[centers_tumbl_oehmke,count_tumbl_oehmke]=compute_pdf(data,'BinNumber',21,'pdf');
plot(centers_tumbl_oehmke,count_tumbl_oehmke,'b.-','MarkerSize',15)

%% PDF square rotation rates normalised by the mean to compare to enstrophy
% omega^2/ <omega^2>
figure(); hold on; grid on; box on; set(gcf,'Position',[50 50 400 400])
ylabel('p.d.f')
xlabel('$\omega^2 / \langle \omega^2 \rangle$')
set(gca,'Yscale','log','FontSize',20)
squared_omega_space = omega_s_x(:).^2 + omega_s_y(:).^2 +  omega_s_z(:).^2;
squared_omega_body = omega_b_x(:).^2 + omega_b_y(:).^2 +  omega_b_z(:).^2;


data = squared_omega_space/mean(squared_omega_space,'omitnan');
[centers_square_omega_over_mean_space,count_square_omega_over_mean_space] = compute_pdf(data,'BinNumber',50,'pdf');
p1=plot(centers_square_omega_over_mean_space,count_square_omega_over_mean_space);

data = squared_omega_body/mean(squared_omega_body,'omitnan');
[centers_square_omega_over_mean_body,count_square_omega_over_mean_body] = compute_pdf(data,'BinNumber',50,'pdf');
p2=plot(centers_square_omega_over_mean_body,count_square_omega_over_mean_body);

data = omega_s_x(:).^2/mean(omega_s_x(:).^2,'omitnan');
[centers_square_omega_over_mean_x,count_square_omega_over_mean_x] = compute_pdf(data,'BinNumber',50,'pdf');
p3=plot(centers_square_omega_over_mean_x,count_square_omega_over_mean_x);

data = omega_s_y(:).^2/mean(omega_s_y(:).^2,'omitnan');
[centers_square_omega_over_mean_y,count_square_omega_over_mean_y] = compute_pdf(data,'BinNumber',50,'pdf');
p5=plot(centers_square_omega_over_mean_y,count_square_omega_over_mean_y);

data = omega_s_z(:).^2/mean(omega_s_z(:).^2,'omitnan');
[centers_square_omega_over_mean_z,count_square_omega_over_mean_z] = compute_pdf(data,'BinNumber',50,'pdf');
p6=plot(centers_square_omega_over_mean_z,count_square_omega_over_mean_z);

data = omega_b_x(:).^2/mean(omega_b_x(:).^2,'omitnan');
[centers_square_omega_over_mean_spinn,count_square_omega_over_mean_spinn] = compute_pdf(data,'BinNumber',50,'pdf');
p7=plot(centers_square_omega_over_mean_spinn,count_square_omega_over_mean_spinn);

data = (omega_b_y(:).^2 + omega_b_z(:).^2) /mean(omega_b_y(:).^2 + omega_b_z(:).^2,'omitnan');
[centers_square_omega_over_mean_tumbl,count_square_omega_over_mean_tumbl] = compute_pdf(data,'BinNumber',50,'pdf');
p8=plot(centers_square_omega_over_mean_tumbl,count_square_omega_over_mean_tumbl);

legend([p1 p2 p3 p5 p6 p7 p8],{'$\omega$ space frame', '$\omega$ body frame',...
                '$\omega_x$','$\omega_y$','$\omega_z$',...
                '$\omega_s$','$\omega_t$'},'FontSize',20,'location','best')

%% orientation red vector pdf
figure()
set(gcf,'Position',[50 50 1200 600])

% red streamwise component (p_x)
subplot(1,3,1)
nbins=10;
data = reshape(abs(red(:,:,1)),1,[]);
data(isnan(data))=[];

[n_red_x,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_red_x=(edges(1:end-1)+edges(2:end))/2;
plot(centers_red_x,n_red_x,'r.-','MarkerSize',15)

set(gca,'FontSize',20,'YScale','lin');
ylabel('p.d.f')
xlabel('$|p_x|$')

% red spanwise component (p_y)
subplot(1,3,2)
nbins=10;
data = reshape(abs(red(:,:,2)),1,[]);
data(isnan(data))=[];

[n_red_y,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_red_y=(edges(1:end-1)+edges(2:end))/2;
plot(centers_red_y,n_red_y,'r.-','MarkerSize',15)

set(gca,'FontSize',20,'YScale','lin');
ylabel('p.d.f')
xlabel('$|p_y|$')

% red wall-normal component (p_z)
subplot(1,3,3)
nbins=10;
data = reshape(abs(red(:,:,3)),1,[]);
data(isnan(data))=[];

[n_red_z,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_red_z=(edges(1:end-1)+edges(2:end))/2;
plot(centers_red_z,n_red_z,'r.-','MarkerSize',15)

set(gca,'FontSize',20,'YScale','lin');
ylabel('p.d.f')
xlabel('$|p_z|$')

%% orientation vectors over y+
wall_dist = reshape(wall_distance_fibre,1,[]);

%numBins_z=15;

[bin_red_x,mean_red_x,ci_red_x,ci_mean_red_x] = Vlad_bin_mean(wall_dist,reshape(abs(red(:,:,1)),1,[]),numBins_z,'lin',alpha);
[bin_red_x,mean_red_y,ci_red_x,ci_mean_red_y] = Vlad_bin_mean(wall_dist,reshape(abs(red(:,:,2)),1,[]),numBins_z,'lin',alpha);
[bin_red_x,mean_red_z,ci_red_x,ci_mean_red_z] = Vlad_bin_mean(wall_dist,reshape(abs(red(:,:,3)),1,[]),numBins_z,'lin',alpha);

[bin_red_x,mean_green_x,ci_red_x,ci_mean_green_x] = Vlad_bin_mean(wall_dist,reshape(abs(green(:,:,1)),1,[]),numBins_z,'lin',alpha);
[bin_red_x,mean_green_y,ci_red_x,ci_mean_green_y] = Vlad_bin_mean(wall_dist,reshape(abs(green(:,:,2)),1,[]),numBins_z,'lin',alpha);
[bin_red_x,mean_green_z,ci_red_x,ci_mean_green_z] = Vlad_bin_mean(wall_dist,reshape(abs(green(:,:,3)),1,[]),numBins_z,'lin',alpha);

[bin_red_x,mean_blue_x,ci_red_x,ci_mean_blue_x] = Vlad_bin_mean(wall_dist,reshape(abs(blue(:,:,1)),1,[]),numBins_z,'lin',alpha);
[bin_red_x,mean_blue_y,ci_red_x,ci_mean_blue_y] = Vlad_bin_mean(wall_dist,reshape(abs(blue(:,:,2)),1,[]),numBins_z,'lin',alpha);
[bin_red_x,mean_blue_z,ci_red_x,ci_mean_blue_z] = Vlad_bin_mean(wall_dist,reshape(abs(blue(:,:,3)),1,[]),numBins_z,'lin',alpha);

% plot
figure(); set(gcf,'Position',[900 10 600 900]);grid on; box on; hold on;
subplot(3,1,1);grid on; box on; hold on;
set(gca,'XScale','log','FontSize',10)
xlabel('$z^+$')
ylabel('$e_{1j}$')
ylim([0 1])
plot(bin_red_x/1e3/scales.viscous_length,mean_red_x,'k-','LineWidth',2)
plot(bin_red_x/1e3/scales.viscous_length,mean_red_y,'k--','LineWidth',2)
plot(bin_red_x/1e3/scales.viscous_length,mean_red_z,'k:','LineWidth',2)
legend('$\vec{e_1} \cdot \vec{x}$','$\vec{e_1} \cdot \vec{y}$','$\vec{e_1} \cdot \vec{z}$','location','best','FontSize',10)

subplot(3,1,2);grid on; box on; hold on;
set(gca,'XScale','log','FontSize',10)
xlabel('$z^+$')
ylabel('$e_{2j}$')
ylim([0 1])
plot(bin_red_x/1e3/scales.viscous_length,mean_green_x,'k-','LineWidth',2)
plot(bin_red_x/1e3/scales.viscous_length,mean_green_y,'k--','LineWidth',2)
plot(bin_red_x/1e3/scales.viscous_length,mean_green_z,'k:','LineWidth',2)
legend('$\vec{e_2} \cdot \vec{x}$','$\vec{e_2} \cdot \vec{y}$','$\vec{e_2} \cdot \vec{z}$','location','best','FontSize',10)

subplot(3,1,3);grid on; box on; hold on;
set(gca,'XScale','log','FontSize',10)
xlabel('$z^+$')
ylabel('$e_{3j}$')
ylim([0 1])
%xlim([1 1000])
plot(bin_red_x/1e3/scales.viscous_length,mean_blue_x,'k-','LineWidth',2)
plot(bin_red_x/1e3/scales.viscous_length,mean_blue_y,'k--','LineWidth',2)
plot(bin_red_x/1e3/scales.viscous_length,mean_blue_z,'k:','LineWidth',2)
legend('$\vec{e_3} \cdot \vec{x}$','$\vec{e_3} \cdot \vec{y}$','$\vec{e_3} \cdot \vec{z}$','location','best','FontSize',10)

title('Mean components of the orientation vectors over wall-normal distance','FontSize',10)
%legend('$\langle \omega_t^2 \rangle \tau$','$\langle \omega_s^2 \rangle \tau$')

%% orientation vectors over y+ using euler angles
% wall distance
wall_dist = reshape(wall_distance_fibre,1,[]);

% number of bins in y direction
%numBins_z=30;

% compute the euler angles v1 - slow
% R=NaN(size(blue,1),size(blue,2),3,3);
% eulerAngles=NaN(size(blue,1),size(blue,2),3);
% for it = 1:size(blue,2)
%     for ij = 1:size(blue,1)
%         R(ij,it,:,:) = [squeeze(red(ij,it,:)), squeeze(green(ij,it,:)), squeeze(blue(ij,it,:))];
%         eulerAngles(ij,it,:) = rotm2eul(R(ij,it,:,:));
%     end
% end

% compute the euler angles v2 - fast
R=[];eulerAngles=[];
red_x = reshape(red(:,:,1),1,[]);
red_y = reshape(red(:,:,2),1,[]);
red_z = reshape(red(:,:,3),1,[]);

green_x = reshape(green(:,:,1),1,[]);
green_y = reshape(green(:,:,2),1,[]);
green_z = reshape(green(:,:,3),1,[]);

blue_x = reshape(blue(:,:,1),1,[]);
blue_y = reshape(blue(:,:,2),1,[]);
blue_z = reshape(blue(:,:,3),1,[]);

R(1,1,:)=red_x;
R(2,1,:)=red_y;
R(3,1,:)=red_z;

R(1,2,:)=green_x;
R(2,2,:)=green_y;
R(3,2,:)=green_z;

R(1,3,:)=blue_x;
R(2,3,:)=blue_y;
R(3,3,:)=blue_z;


eulerAngles = rotm2eul(R);


% compute 
[bin_eul,mean_eul_1,ci_bin_eul,ci_mean_eul_1] = Vlad_bin_mean(wall_dist,reshape(abs(eulerAngles(:,1)),1,[]),numBins_z,'lin',alpha);
[bin_eul,mean_eul_2,ci_bin_eul,ci_mean_eul_2] = Vlad_bin_mean(wall_dist,reshape(abs(eulerAngles(:,2)),1,[]),numBins_z,'lin',alpha);
[bin_eul,mean_eul_3,ci_bin_eul,ci_mean_eul_3] = Vlad_bin_mean(wall_dist,reshape(abs(eulerAngles(:,3)),1,[]),numBins_z,'lin',alpha);
figure()
plot(bin_eul,mean_eul_1);hold on;
plot(bin_eul,mean_eul_2)
plot(bin_eul,mean_eul_3)

%% orientation vector components jpdf with y+
wall_dist = wall_distance_fibre;
wall_dist= reshape(wall_dist,1,[]);
x=wall_dist/1e3/scales.viscous_length;

bin_nr_x = 10;
bin_nr_vect = 10;

figure(); set(gcf,'Position',[500 10 900 900]);grid on; box on; hold on;
subplot(3,3,1);box on;
y=reshape(abs(red(:,:,1)),1,[]);
[counts,YY,XX]=compute_jpdf(x,y,bin_nr_x,bin_nr_vect);
contourf(YY,XX,counts,'EdgeColor','none')
xlabel('$z^+$')
ylabel('$\vec{e_1} \cdot \vec{x}$')
colormap(parula(512))
colorbar('eastoutside')
clim([0 4e-3])

subplot(3,3,2);box on;
y=reshape(abs(red(:,:,2)),1,[]);
[counts,YY,XX]=compute_jpdf(x,y,bin_nr_x,bin_nr_vect);
contourf(YY,XX,counts,'EdgeColor','none')
xlabel('$z^+$')
ylabel('$\vec{e_1} \cdot \vec{y}$')
colormap(parula(512))
%colorbar('northoutside')
clim([0 4e-3])

subplot(3,3,3);box on;
y=reshape(abs(red(:,:,3)),1,[]);
[counts,YY,XX]=compute_jpdf(x,y,bin_nr_x,bin_nr_vect);
contourf(YY,XX,counts,'EdgeColor','none')
xlabel('$z^+$')
ylabel('$\vec{e_1} \cdot \vec{z}$')
colormap(parula(512))
%colorbar('northoutside')
clim([0 4e-3])

subplot(3,3,4);box on;
y=reshape(abs(green(:,:,1)),1,[]);
[counts,YY,XX]=compute_jpdf(x,y,bin_nr_x,bin_nr_vect);
contourf(YY,XX,counts,'EdgeColor','none')
xlabel('$z^+$')
ylabel('$\vec{e_2} \cdot \vec{x}$')
colormap(parula(512))
%colorbar('northoutside')
clim([0 4e-3])

subplot(3,3,5);box on;
y=reshape(abs(green(:,:,2)),1,[]);
[counts,YY,XX]=compute_jpdf(x,y,bin_nr_x,bin_nr_vect);
contourf(YY,XX,counts,'EdgeColor','none')
xlabel('$z^+$')
ylabel('$\vec{e_2} \cdot \vec{y}$')
colormap(parula(512))
%colorbar('northoutside')
clim([0 4e-3])

subplot(3,3,6);box on;
y=reshape(abs(green(:,:,3)),1,[]);
[counts,YY,XX]=compute_jpdf(x,y,bin_nr_x,bin_nr_vect);
contourf(YY,XX,counts,'EdgeColor','none')
xlabel('$z^+$')
ylabel('$\vec{e_2} \cdot \vec{z}$')
colormap(parula(512))
%colorbar('northoutside')
clim([0 4e-3])

subplot(3,3,7);box on;
y=reshape(abs(blue(:,:,1)),1,[]);
[counts,YY,XX]=compute_jpdf(x,y,bin_nr_x,bin_nr_vect);
contourf(YY,XX,counts,'EdgeColor','none')
xlabel('$z^+$')
ylabel('$\vec{e_3} \cdot \vec{x}$')
colormap(parula(512))
%colorbar('northoutside')
clim([0 4e-3])

subplot(3,3,8);box on;
y=reshape(abs(blue(:,:,2)),1,[]);
[counts,YY,XX]=compute_jpdf(x,y,bin_nr_x,bin_nr_vect);
contourf(YY,XX,counts,'EdgeColor','none')
xlabel('$z^+$')
ylabel('$\vec{e_3} \cdot \vec{y}$')
colormap(parula(512))
%colorbar('northoutside')
clim([0 4e-3])

subplot(3,3,9);box on;
y=reshape(abs(blue(:,:,3)),1,[]);
[counts,YY,XX]=compute_jpdf(x,y,bin_nr_x,bin_nr_vect);
contourf(YY,XX,counts,'EdgeColor','none')
xlabel('$z^+$')
ylabel('$\vec{e_3} \cdot \vec{z}$')
colormap(parula(512))
%colorbar('northoutside')
clim([0 4e-3])

%% principal axis lengths

figure(); hold on; grid minor; box on;
set(gcf,'Position',[50 50 1200 500])

% longest axis length L1
nbins=20;
data = reshape(L1,1,[]);
data(isnan(data))=[];

data = data.*p.dx; % scale to milimeters

[n_L1,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_L1=(edges(1:end-1)+edges(2:end))/2;
p1=plot(centers_L1,n_L1,'s-','MarkerSize',10,'Linewidth',2,'Color','k','MarkerFaceColor','k');


% medium axis length L2

nbins=12;
data = reshape(L2,1,[]);
data(isnan(data))=[];

data = data.*p.dx; % scale to milimeters

[n_L2,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_L2=(edges(1:end-1)+edges(2:end))/2;
p2=plot(centers_L2,n_L2,'^-','MarkerSize',10,'Linewidth',2,'Color','b','MarkerFaceColor','b');

% shortest axis length L3
nbins=7;
data = reshape(L3,1,[]);
data(isnan(data))=[];

data = data.*p.dx; % scale to milimeters

[n_L3,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_L3=(edges(1:end-1)+edges(2:end))/2;
p3=plot(centers_L3,n_L3,'.-','MarkerSize',35,'Linewidth',2,'Color','r','MarkerFaceColor','r');


legend([p1,p2,p3],{'$l_1$','$l_2$','$l_3$'},'FontSize',20)
set(gca,'FontSize',20,'YScale','lin');
ylabel('p.d.f')
xlabel('$l_1$, $l_2$, $l_3$ [mm]')
xlim([0 2])


%% principal axis lengths scaled by kolmogorov


L1_kol = L1*p.dx./(local_l_kol_fibre*1e3);
L2_kol = L2*p.dx./(local_l_kol_fibre*1e3);
L3_kol = L3*p.dx./(local_l_kol_fibre*1e3);

figure(); hold on; grid minor; box on;
set(gcf,'Position',[50 50 1200 500])

% longest axis length L1
nbins=30;
data = reshape(L1_kol,1,[]);
data(isnan(data))=[];

[n_L1_kol,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_L1_kol=(edges(1:end-1)+edges(2:end))/2;
p1=plot(centers_L1_kol,n_L1_kol,'s-','MarkerSize',10,'Linewidth',2,'Color','k','MarkerFaceColor','k');

% medium axis length L2

nbins=20;
data = reshape(L2_kol,1,[]);
data(isnan(data))=[];

[n_L2_kol,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_L2_kol=(edges(1:end-1)+edges(2:end))/2;
p2=plot(centers_L2_kol,n_L2_kol,'^-','MarkerSize',10,'Linewidth',2,'Color','b','MarkerFaceColor','b');

% shortest axis length L3

nbins=10;
data = reshape(L3_kol,1,[]);
data(isnan(data))=[];

[n_L3_kol,edges]=histcounts(data,nbins,"Normalization","probability"); 
centers_L3_kol=(edges(1:end-1)+edges(2:end))/2;
p3=plot(centers_L3_kol,n_L3_kol,'.-','MarkerSize',35,'Linewidth',2,'Color','r','MarkerFaceColor','r');


xlim([0 20])

legend([p1,p2,p3],{'$l_1$','$l_2$','$l_3$'},'FontSize',20)
set(gca,'FontSize',20,'YScale','lin');
ylabel('p.d.f')
xlabel('$l_1/\eta$, $l_2/\eta$, $l_3/\eta$ [-]')


%% jpdf tumbling and spinning
% omega_tumbling_visc = abs(omega_tumbling_a * scales.viscous_time);
% omega_spinning_visc = abs(omega_spinning_a * scales.viscous_time);
% 
% omega_tumbling_visc(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];
% omega_spinning_visc(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];
% figure()
% [counts,TUMB,SPINN]=compute_jpdf(reshape(omega_tumbling_visc,1,[]),reshape(omega_spinning_visc,1,[]),50,50);
% contourf(TUMB,SPINN,counts)
% ylabel('$\omega_s\tau$')
% xlabel('$\omega_t\tau$')
% ccc=colorbar;
% ccc.Label.String='jpdf';

%% jpdf tumbling or spinning with velocity
% omega_tumbling_visc = abs(omega_tumbling_a * scales.viscous_time);
% omega_spinning_visc = abs(omega_spinning_a * scales.viscous_time);
% 
% omega_tumbling_visc(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];
% omega_spinning_visc(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];
% %x_dot(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];
% %y_dot(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];
% %z_dot(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];
% 
% figure()
% [counts,Z_DOT,SPINN]=compute_jpdf(reshape(z_dot,1,[]),reshape(omega_spinning_visc,1,[]),50,50);
% contourf(Z_DOT,SPINN,counts)
% ylabel('$\omega_s\tau$')
% xlabel('$\dot{z}$')
% ccc=colorbar;
% ccc.Label.String='jpdf';

%% which subset of fibres show extreme events in the spinning rate
omega_spinning_sub = omega_spinning_a.^2./mean(mean(omega_spinning_a.^2,2,'omitnan'),'omitnan');
omega_spinning_sub(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];

curv=fibre_curvature;
curv(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];

wld=wall_distance_fibre;
wld(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];

len=fibre_length;
len(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];

z_dot_sub = z_dot;
z_dot_sub(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];

omega_tumbling_sub = omega_tumbling_a.^2./mean(mean(omega_tumbling_a.^2,2,'omitnan'),'omitnan');
omega_tumbling_sub(std_fibre_curvature>std_fibre_curvature_threshold,:)=[];

% % check components of eigenvectors
% red_sub = red;
% red_sub(std_fibre_curvature>std_fibre_curvature_threshold,:,:)=[];
% green_sub = green;
% green_sub(std_fibre_curvature>std_fibre_curvature_threshold,:,:)=[];
% blue_sub = blue;
% blue_sub(std_fibre_curvature>std_fibre_curvature_threshold,:,:)=[];
% figure(1)
% for rrr=row'
%     clf
%     subplot(1,3,1); hold on; grid on; box on;
%     plot(red(rrr,:,1),'r.-')
%     plot(red(rrr,:,2),'g.-')
%     plot(red(rrr,:,3),'b.-')
% 
% 
%     subplot(1,3,2); hold on; grid on; box on;
%     plot(green(rrr,:,1),'r.-')
%     plot(green(rrr,:,2),'g.-')
%     plot(green(rrr,:,3),'b.-')
% 
%     subplot(1,3,3); hold on; grid on; box on;
%     plot(blue(rrr,:,1),'r.-')
%     plot(blue(rrr,:,2),'g.-')
%     plot(blue(rrr,:,3),'b.-')
%     pause
% end

figg=figure(); figg.Position=[100 100 1400 500];
subplot(1,4,1); hold on; grid on; box on;
[centers,count]=compute_pdf(reshape(curv,1,[]),'BinWidth',0.05,'pdf');
plot(centers,count,'k-')
[centers,count]=compute_pdf(reshape(curv(omega_spinning_sub>8),1,[]),'BinWidth',0.05,'pdf');
plot(centers,count,'-','Color',colors_vlad('blue','dark'))
xlabel('$\kappa^*$')
ylabel('pdf')
legend('all data','$\omega_s^2/\langle \omega_s^2 \rangle > 8$','Location','south')
set(gca,'FontSize',15)

subplot(1,4,2); hold on; grid on; box on;
[centers,count]=compute_pdf(reshape(len,1,[]),'BinWidth',0.1,'pdf');
plot(centers,count,'k-')
[centers,count]=compute_pdf(reshape(len(omega_spinning_sub>8),1,[]),'BinWidth',0.1,'pdf');
plot(centers,count,'-','Color',colors_vlad('blue','dark'))
xlabel('$L_f$ [mm]')
ylabel('pdf')
%legend('all data','$\omega_s^2/\langle \omega_s^2 \rangle > 8$')
set(gca,'FontSize',15)

subplot(1,4,3); hold on; grid on; box on;
[centers,count]=compute_pdf(reshape(wld,1,[]),'BinWidth',2,'pdf');
plot(centers,count,'k-')
[centers,count]=compute_pdf(reshape(wld(omega_spinning_sub>8),1,[]),'BinWidth',2,'pdf');
plot(centers,count,'-','Color',colors_vlad('blue','dark'))
xlabel('$y$ [mm]')
ylabel('pdf')
%legend('all data','$\omega_s^2/\langle \omega_s^2 \rangle > 8$')
set(gca,'FontSize',15)

subplot(1,4,4); hold on; grid on; box on;
[centers,count]=compute_pdf(reshape(omega_tumbling_sub,1,[]),'BinWidth',1,'pdf');
plot(centers,count,'k-')
[centers,count]=compute_pdf(reshape(omega_tumbling_sub(omega_spinning_sub>8),1,[]),'BinWidth',1,'pdf');
plot(centers,count,'-','Color',colors_vlad('blue','dark'))
xlabel('$\omega_t^2/\langle \omega_t^2 \rangle$')
ylabel('pdf')
%legend('all data','$\omega_s^2/\langle \omega_s^2 \rangle > 8$')
set(gca,'YScale','log','FontSize',15)


%% save data

save(strcat(fname.main,'statistics_',which_data,'.mat'),'centers_spinning','centers_tumbling','n_spinning','n_tumbling',...
    'bin_z_Centers_rot','scales',"mean_tumbl","mean_spinn",'bin_z_Centers_vel','mean_x_velocities',...
    'bin_k_star_Centers',"mean_tumbl_k_star","mean_spinn_k_star",'bin_z_Centers_rot_visc','mean_tumbl_visc','mean_spinn_visc',...
    'bin_l_eta_Centers','mean_tumbl_l_eta','mean_spinn_l_eta',...
    'bin_l_eta_Centers_visc','mean_tumbl_l_eta_visc','mean_spinn_l_eta_visc',...
    'mean_squared_tumbling_visc','mean_squared_spinning_visc',...
    'n_red_x','n_red_y','n_red_z','centers_red_x','centers_red_y','centers_red_z',...
    'centers_L3','n_L3','centers_L2','n_L2','centers_L1','n_L1','centers_k','count_k',...
    'centers_tumbling_high_k','n_tumbling_high_k','centers_spinning_high_k','n_spinning_high_k',...
    'centers_tumbling_med_k','n_tumbling_med_k','centers_spinning_med_k','n_spinning_med_k',...
    'centers_tumbling_low_k','n_tumbling_low_k','centers_spinning_low_k','n_spinning_low_k',...
    'bin_z_Centers_rot_tumbl_2_3_visc',...
    'mean_tumbl_2_visc','mean_tumbl_3_visc',...
    'bin_z_Centers_rot_x_y_z_visc',...
    'mean_omega_x_visc','mean_omega_y_visc','mean_omega_z_visc',...
    'bin_z_Centers_rot_proj_tumbl_2_3_visc',...
    'mean_proj_tumbl_2_x_visc','mean_proj_tumbl_2_y_visc','mean_proj_tumbl_2_z_visc', ...
    'mean_proj_tumbl_3_x_visc','mean_proj_tumbl_3_y_visc','mean_proj_tumbl_3_z_visc',...
    'bin_z_Centers_rot_proj_spinn_1_visc',...
    'mean_proj_spinn_1_x_visc','mean_proj_spinn_1_y_visc','mean_proj_spinn_1_z_visc',...
    'bin_red_x',...
    'mean_red_x','mean_red_y','mean_red_z',...
    'mean_green_x','mean_green_y','mean_green_z',...
    'mean_blue_x','mean_blue_y','mean_blue_z',...
    'centers_spin_oehmke','count_spin_oehmke',...
    'centers_tumbl_oehmke','count_tumbl_oehmke',...
    'bin_eul','mean_eul_1','mean_eul_2','mean_eul_3',...
    'red','green','blue','wall_distance_fibre',...
    'bin_k_star_Centers_visc','mean_tumbl_k_star_visc','ci_k_star_Centers_visc','ci_mean_tumbl_k_star_visc',...
    'mean_spinn_k_star_visc','ci_mean_spinn_k_star_visc',...
    'bin_z_Centers_rot_visc','mean_tumbl_visc','ci_z_Centers_rot_visc','ci_mean_tumbl_visc',...
    'mean_spinn_visc','ci_z_Centers_rot_visc','ci_mean_spinn_visc',...
    'bin_z_Centers_rot_tumbl_2_3_visc','mean_tumbl_2_visc','ci_z_Centers_rot_tumbl_2_3_visc','ci_mean_tumbl_2_visc',...
    'mean_tumbl_3_visc','ci_mean_tumbl_3_visc',...
    'alpha',...
    'mean_tumbl_visc_conv_bins','mean_spinn_visc_conv_bins','nr_traj_used_conv_bins',...
    'concentration_centers','concentration_counts',...
    'mean_omega_body_visc','ci_mean_omega_body_visc',...
    'mean_omega_body_tbl_spin_visc','ci_mean_omega_body_tbl_spin_visc',...
    'mean_l_eta','ci_mean_l_eta',...
    'bin_z_Centers_simple_omega_visc','mean_simple_omega_x_visc','mean_simple_omega_y_visc','mean_simple_omega_z_visc',...
    'mean_simple_omega_1_visc','mean_simple_omega_2_visc','mean_simple_omega_3_visc',...
    'ci_z_Centers_rot_x_y_z_visc','ci_mean_omega_x_visc','ci_mean_omega_y_visc','ci_mean_omega_z_visc',...
    'ci_z_Centers_simple_omega_visc','ci_mean_simple_omega_x_visc','ci_mean_simple_omega_y_visc','ci_mean_simple_omega_z_visc',...
    'ci_mean_simple_omega_1_visc','ci_mean_simple_omega_2_visc','ci_mean_simple_omega_3_visc',...
    'centers_square_omega_over_mean_tumbl','count_square_omega_over_mean_tumbl',...
    'centers_square_omega_over_mean_spinn','count_square_omega_over_mean_spinn',...
    'centers_square_omega_over_mean_x','count_square_omega_over_mean_x',...
    'centers_square_omega_over_mean_y','count_square_omega_over_mean_y',...
    'centers_square_omega_over_mean_z','count_square_omega_over_mean_z',...
    'centers_square_omega_over_mean_space','count_square_omega_over_mean_space',...
    'centers_square_omega_over_mean_body','count_square_omega_over_mean_body',...
    'bin_z_Centers_velocity','mean_velocity_x','mean_velocity_y','mean_velocity_z',...
    'ci_z_Centers_velocity','ci_mean_velocity_x','ci_mean_velocity_y','ci_mean_velocity_z',...
    'avg_fibre_length','avg_fibre_curvature','n_tumbling_half','n_tumbling_quarter','n_tumbling_one_tenth','n_tumbling_one_onehundreth',...
    'p')

%% copy files to remote client
% Define source and destination paths
sourcePath =strcat(fname.main,'statistics_',which_data,'.mat');

switch which_data
    case 'center'
        destinationPath = '\\tsclient\remote_temp\____PRL\svd_rlowess_ts_30_fk_45_centre\';
        
    case 'near_wall' 
        destinationPath = '\\tsclient\remote_temp\____PRL\rlowess_ts_15_fk_45_near_wall\';
end

try
    % Copy the file
    copyfile(sourcePath, destinationPath);
    
    % Display success message
    disp('File copied successfully!');
catch
    % Display error message if the copy operation fails
    disp('Error: Unable to copy the file.');
end

close all
