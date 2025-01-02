clc; clear; close all;
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

addpath('_common\')
addpath("readim7_aux\")

%% parameters
%%% experiment
p.Re_tau=720;               % mean shear Reynolds number
p.Temperature=26.2;           % mean Temperature during the experiment
p.dt=1/801.192;                 % time sepration between two timesteps in s
p.dx=1/49.3;                  % scale factor in mm/vox

%%% computation structure
p.parallel_cores=3;        % Number of cores to be used in parallel for discretisation. 
                            % Decide based on available RAM and available CPU cores (1 object of 1800x1300x700 is approx 12.5 GB of RAM).
                            % Works only if p.plot=0

%%% discretisation
p.in=1;                     % Starting timestep number
p.ki=2;                   % Ending timestep number 

p.imbin_thres=0.5;          % Threshold for binarization, enter 0 if you dont want to consider a threshold
p.l_min=40;                 % Minimum length (in voxels) of the longest principal axis of the equivalent ellipsoid (see regionprops3)

%%% tracking
p.radius=25;                % radius of search in voxel

p.min_track_length=6;      % minimum track length - the rest of the tracks are deleted
p.max_track_length=10;     % maximum track length - the rest of the track is considered new tracks 


%%% refinement
p.peak_finding_technique="Majority";                % Dilation, imregionalmax, max, Skeletonize, Majority, Majority + Skeletonize RECOMMENDED: 'Majority'
p.polynomial_finding_technique="Mobin";             % Mobin - fit polynomial to fibre object (Alipour et al. 2021) RECOMMENDED: 'Mobin'
p.use_weights=0;                                    % use the light intensity for weighting the polynomial fit RECOMMENDED: 0
p.fitting_type_fibre='poly2';                       % polynomial order to use for fiber fitting 
                                                    % if p.polynomial_finding_technique="Mobin" RECOMMENDED: 'poly2'

p.percentile_majority_skeletonize=1;                % applies to "Skeletonize", "Majority", and "Majority + Skeletonize" [percent] RECOMMENDED: 1
p.kernel_dilation=7;                                % applies to "Dilation" RECOMMENDED: 7 
p.min_branch_length_skeletonize=5;                  % applies to "Skeletonize" RECOMMENDED: 5
                                                    % the above 3 parameters can be used to generate points in space which approximate the fibre object
                                                    % if p.peak_finding_technique="Majority" and p.percentile_majority_skeletonize=1 then
                                                    % only the top 99% brightest points from the fibre object are used for fitting the fibre
                                                    

%%% quantities computation
% track fitting
p.fitting_type="rlowess";               %%% to fit a polynomial on all the positions, Euler Angles, components of eigenvectors of orientation within a track use:
                                        % "Poly1", "Poly2", "Poly3","Poly4", "Poly5", "Poly9", "Automatic"
                                        %%% to filter all positions, Euler Angles, components of eigenvectors of orientation within a track use:
                                        % "Marco Sgolay2 - 1", "Marco Smooth 2x" (Alipour et al. 2021),
                                        % 'movmean','movmedian','gaussian', 'lowess', 'loess','rlowess', 'rloess', 'sgolay','spline'
                                        % see description of smoothdata Matlab function

p.kernel_trajectory_positions=20;       % kernel length for filtering the positions in timesteps
p.kernel_trajectory_angles=10;          % kernel length for filtering the Euler Angles in timesteps

p.use_which_vectors_for_orientation='region';  % which eigenvectors to use for computation of the rotation rates
                                            % 'poly' - from the inertia tensor of the polynomial fitted to the fibre
                                            % 'region' - the Eigenvectors from regionprops3 on the fibre object

p.kernel_trajectory_vectors=20;             % kernel length for filtering the components of the eigenvectors in timesteps
p.use_fitted_vectors_for_rotation_matrix=1; % if the filtered components should be used for the rotation rate computation (1) or not (0)
p.skip_timesteps_rotation_matrix=2;         % how many timesteps to skip when computing the derivative of the rotation matrix


% derivative computation
p.derivative_stencil='5 points stencil';    % for computing positions and angles derivatives use a '2 points stencil' or '5 points stencil'
p.disable_edge_points=0;                    % 1 - when computing derivatives the edge points of the trajectory are not computed  - RECOMMENDED: 0

% angle corrections
p.correct_orientations=1;                       % if fiber fixed reference frame is not consistent from frame to frame (if it changes by more than the angle below)
                                                % fix the coordinate system and recompute euler angles and rotational velocities
p.angle=90;                                     % which angle to use for the above method [deg] RECOMMENDED: 90


%%% saving data
p.save_data_tracking=0;     % as tracking and modelling are intermediate steps in the processing
p.save_data_modelling=0;    % we can choose if this intermediate data is saved in p.save folder (large space requirement) - Both RECOMMENDED: 0
p.save_fiber_objects=0;     % if the 3D objects of the fibers are saved in the final processing step 
                            % (large space requirement) - RECOMMENDED: 0 (use 1 if you need to visualise the fibers afterwards) 

%%% processing steps
p.discretize=1;         % which processing steps to perform, useful when changing parameters linked to a particular step and 
p.track=1;              % the previous steps are not required to be recomputed
p.modelling=1;
p.quantities=1;


%%% plotting
p.plot=1;                                           % show each step of the process
p.print=0;                                          % save each plot in p.save\Figures_Processing\
p.pause_enabled=1;                                  % 1 - after each plot pause the program and wait for a key press. 0 - don't wait for key press
p.levellist=linspace(p.imbin_thres,20*p.imbin_thres,5);                            % list of isosurface values of intensity to show the fibers
%p.levellist = [2 4 6];
p.facealphalist=linspace(0.1,0.4,numel(p.levellist)); % list of alpha values for each isosurface of instensity
p.skip=30;                                           % how many time steps to skip when showing the fiber track and eigenvectors (for clarity)

p.make_movie=1;                                     % makes a movie for each timestep of each fiber tracked: camera follows the centroid
p.movie_framerate=5;                               % framerate of movie
p.limits_plot=30;                                   % limits of the plot for the movie: xlim: [x-p.limits_plot, x+p.limits_plot] analogue for y and z

%%% check code with synthetic fibres
p.virtual_data=0; % 1 - if the input data is generated virtually to check the code

%% processing
p.main_folder = 'E:\PRL2_CC_center_v2\Fibers_test\ImgPreproc\TomographicPIV\Data\';
p.subfolders = {'example_data'};
nr_of_loops = 1;

count = 1;
for i_sub = 1:numel(p.subfolders)

    for i_loops = 1:nr_of_loops
        p.load = strcat(p.main_folder,p.subfolders{i_sub},'\loop=',num2str(i_loops-1),'\Data\');
        p.save = strcat(p.main_folder,'\example_data\_processed_data\loop=',num2str(count-1),'\');

        p.load = p.main_folder;
        p.save = p.main_folder;
        % processing step
        Vlad_main_for_loops(p)
        
        count = count+1;
    end

end