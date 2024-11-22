function []=Vlad_main_for_loops(p)

%%%%%%% ------ Main program to: ------ %%%%%%%
% 1. DISCRETIZE fibers from tracers and noise (MART exported objects -> fiber objects)
% 2. TRACK fibers (fiber objects -> bundle of fiber objects ob each track)
% 3. MODEL fibers (fiber objects -> polynomial fit for each object)
% 4. COMPUTE fiber quantities (polynomial fit -> position, velocity, rotation rates, length, curvature etc.)
        
% save parameters in folder
if ~exist(p.save,'dir'); mkdir(p.save); end
save(strcat(p.save,'parameters.mat'),"p")

%% DISCRETIZE
% Rread .im7 Davis files containing planes generated with MART. Fibers are discretised based on intensity threshold and other thresholds such as length of the longest principal axis (regionprops3).
% For each timestep a struct "Fibers" is saved. It contains the properties of each fiber. Each row corresponds to a unique fiber. Columns correspond to the timestep.

if p.discretize==1

p.ik=count_files(strcat(p.load,'S00001'),'.im7');                   % Find the number of 2D planes in each timestep; each plane is stored as one .im7 file

if p.virtual_data==0
    %%% if plotting is not required for initial parameter setup, then run discretisation in parallel
    if p.plot==0
        % create parallel worker pool
        parpool(p.parallel_cores);        
    
        % create parallel for loop: each core workes on one timestep at a time
        f1=p.load;      % initial asignment for less overhead
        ik=p.ik;
        parfor kk=p.in:p.ki
            [nr_fibers,Fibers_All{kk}]=Vlad_fiber_discrimination_3D_v4(permute(read_im7(f1,kk,ik),[2 1 3]),p,kk); 
            % permute is done so x direction is along the first dimension of the
            % instensity matrix (e.g. I(:,1,1)) and y direction is along the second
            % direction of the intensity matrix (e.g. I(1,:,1))
    
            disp(strcat('Processed object number: ',num2str(kk),' containing:',num2str(nr_fibers),' Fibers')) % show progress
        end
        
        % stop the parallel worker pool
        delete(gcp('nocreate'))
    else
        f1=p.load;
        ik=p.ik;
        for kk=p.in:p.ki
            [nr_fibers,Fibers_All{kk}]=Vlad_fiber_discrimination_3D_v4(permute(read_im7(f1,kk,ik),[2 1 3]),p,kk); % discretize
            disp(strcat('Processed object number: ',num2str(kk),' containing:',num2str(nr_fibers),' Fibers')) % show progress
        end
    end
elseif p.virtual_data==1    
    %load(strcat(p.load,'fibers.mat'))
    p.dx = 1;
    p.dt = 1;
    for kk=p.in:p.ki
        load(strcat(p.load,'I_',num2str(kk),'.mat'))
        [nr_fibers,Fibers_All{kk}]=Vlad_fiber_discrimination_3D_v4(full(I(:,:,:)),p,kk); % discretize
        clear I
        disp(strcat('Processed object number: ',num2str(kk),' containing:',num2str(nr_fibers),' Fibers')) % show progress
    end
end
%%% save data
disp('Saving all fibers to disc')
% make folder for saving 
if ~exist(strcat(p.save,'1_FiberOnly\'),'dir')
    mkdir(strcat(p.save,'1_FiberOnly\'))
end
% save all fibers in each timestep separately
for kk=p.in:p.ki
    Fibers=Fibers_All{kk};
    fn=strcat(p.save,'1_FiberOnly\','fiberonly_',num2str(kk),'.mat');
    save(fn,'Fibers')
end

end

%% TRACK
if p.track==1

clc; clearvars -except p

% initialize
yxz_positions=[];
II=zeros(1,1,1);

% track
for kk=p.in:p.ki
    if kk==p.in
        %%% generate database of fibres from the first instance
        load(strcat(p.save,'1_FiberOnly\','fiberonly_',num2str(kk),'.mat'))
        AllFibers=Fibers;
        clear Fibers
    else
        %%% collect all locations of fibers from this instance
        yxz_positions=[];
        nr_fibers=sum(cellfun('isempty',AllFibers.Centroid(:,kk-1))==0);
        for ii=find(~cellfun('isempty',AllFibers.Centroid(:,kk-1)))'
            yxz_positions(ii,1)=AllFibers.Centroid{ii,kk-1}(1);
            yxz_positions(ii,2)=AllFibers.Centroid{ii,kk-1}(2);
            yxz_positions(ii,3)=AllFibers.Centroid{ii,kk-1}(3);
        end

        yxz_positions(and(yxz_positions(:,1)==0,yxz_positions(:,2)==0),1:3)=NaN;

        ptCloud=pointCloud(yxz_positions); % generate point could of positions in previous timestep
        %%% load fibres from next timestep
        load(strcat(p.save,'1_FiberOnly\','fiberonly_',num2str(kk),'.mat'))
        %%% go through list of all fibers in new timestep
        nr_fibers=size(Fibers.Centroid,1);
        for ii=1:nr_fibers
            point=[Fibers.Centroid{ii,1}(1),Fibers.Centroid{ii,1}(2),Fibers.Centroid{ii,1}(3)]; % current fiber position in new timestep
            [indices,distances]=findNeighborsInRadius(ptCloud,point,p.radius);                    % find all distances in interrogation radius
            if ~isempty(distances) %&& size(distances,1)==1                                      % if match is found and it is unique in the interrogation radius
                                                                                                % add more conditions to constrain the tracking if needed
                found=select(ptCloud,indices);
                previous_position=found.Location(1,:);
                [~,id_previous_position]=ismember(previous_position,yxz_positions,'rows');
                track_length_until_now=sum(cellfun('isempty',AllFibers.Centroid(id_previous_position,1:kk-1))==0,2);

                if track_length_until_now<=p.max_track_length
                    %%% asign properties of fiber of new timestep to the same
                    %%% fiber from the old timestep
                    %%% do this for all properties
                    AllFibers.Centroid{id_previous_position,kk}=Fibers.Centroid{ii};    
                    AllFibers.Orientation{id_previous_position,kk}=Fibers.Orientation{ii};
                    AllFibers.Object{id_previous_position,kk}=Fibers.Object{ii};
                    AllFibers.EigenVectors{id_previous_position,kk}=Fibers.EigenVectors{ii};
                    AllFibers.BoundingBox{id_previous_position,kk}=Fibers.BoundingBox{ii};
                else
                    %%% begin new track
                    AllFibers.Centroid{size(AllFibers.Centroid(:,kk-1),1)+1,kk}=Fibers.Centroid{ii};
                    AllFibers.Orientation{size(AllFibers.Orientation(:,kk-1),1)+1,kk}=Fibers.Orientation{ii};
                    AllFibers.Object{size(AllFibers.Object(:,kk-1),1)+1,kk}=Fibers.Object{ii};
                    AllFibers.EigenVectors{size(AllFibers.EigenVectors(:,kk-1),1)+1,kk}=Fibers.EigenVectors{ii};
                    AllFibers.BoundingBox{size(AllFibers.BoundingBox(:,kk-1),1)+1,kk}=Fibers.BoundingBox{ii};
                end

            else
                %%% if no match is found or more than one match is found in
                %%% the interrogation radius the consider the current fiber
                %%% in the new timestep as unique and asign it to a new row

                AllFibers.Centroid{size(AllFibers.Centroid(:,kk-1),1)+1,kk}=Fibers.Centroid{ii};
                AllFibers.Orientation{size(AllFibers.Orientation(:,kk-1),1)+1,kk}=Fibers.Orientation{ii};
                AllFibers.Object{size(AllFibers.Object(:,kk-1),1)+1,kk}=Fibers.Object{ii};
                AllFibers.EigenVectors{size(AllFibers.EigenVectors(:,kk-1),1)+1,kk}=Fibers.EigenVectors{ii};
                AllFibers.BoundingBox{size(AllFibers.BoundingBox(:,kk-1),1)+1,kk}=Fibers.BoundingBox{ii};
            end
        end
    end

    %%% delete tracks shorter than min_track_length
    index_ended_tracks = find(cellfun('isempty',AllFibers.Centroid(:,kk))==1);
    track_length_ended_tracks = sum(cellfun('isempty',AllFibers.Centroid(index_ended_tracks,1:kk))==0,2);
    index_ended_tracks_shorter_than_min = find(track_length_ended_tracks<p.min_track_length);

    index_tracks_to_delete = index_ended_tracks(index_ended_tracks_shorter_than_min);

    AllFibers.Centroid(index_tracks_to_delete,:)=[];
    AllFibers.Orientation(index_tracks_to_delete,:)=[];
    AllFibers.Object(index_tracks_to_delete,:)=[];
    AllFibers.EigenVectors(index_tracks_to_delete,:)=[];
    AllFibers.BoundingBox(index_tracks_to_delete,:)=[];


    % display progress
    disp(strcat('Processed timestep: ',num2str(kk)))
end

%%% delete tracks shorter than min_track_length
track_length=sum(cellfun('isempty',AllFibers.Centroid(:,:))==0,2);
AllFibers.Centroid(track_length<p.min_track_length,:)=[];
AllFibers.Orientation(track_length<p.min_track_length,:)=[];
AllFibers.Object(track_length<p.min_track_length,:)=[];
AllFibers.EigenVectors(track_length<p.min_track_length,:)=[];
AllFibers.BoundingBox(track_length<p.min_track_length,:)=[];

disp(strcat('Tracks found: ',num2str(size(AllFibers.Centroid,1))))

%%% save data
% make folder for saving 
if ~exist(strcat(p.save,'2_TrackedFibers\'),'dir')
    mkdir(strcat(p.save,'2_TrackedFibers\'))
end

if p.save_data_tracking==1
    % save as one big file
    save(strcat(p.save,'2_TrackedFibers\',"AllFibers.mat"),"AllFibers","p")
end

%%% plot figures of tracks
if p.plot==1
    plot_tracks(p,AllFibers)
end

end

%% MODELLING
if p.modelling==1

if p.track==0
    load(strcat(p.save,'2_TrackedFibers\',"AllFibers.mat"),"AllFibers")
end

clc; clearvars -except p AllFibers
% initilize
nr_fibers=size(AllFibers.Centroid,1);
v=[];
% refine
for ii=1:nr_fibers
    if and(p.make_movie==1,p.plot==1)
        if ~exist(strcat(p.save,'Figures_Movies_Processing\4_Modelling\'),'dir')
            mkdir(strcat(p.save,'Figures_Movies_Processing\4_Modelling\'))
        end
        v = VideoWriter(strcat(p.save,'Figures_Movies_Processing\4_Modelling\Fiber_Fitting_',num2str(ii),'_Video.mp4'),'MPEG-4');
        v.FrameRate=p.movie_framerate;
        v.Quality=100;
        open(v)
    end

    for kk=find(~cellfun('isempty',AllFibers.Centroid(ii,:)))
        [AllFibers.Centroid_Refined{ii,kk},...
            AllFibers.Curvature{ii,kk},...
            AllFibers.Length{ii,kk},...
            AllFibers.Length_effective{ii,kk},...
            AllFibers.EigenVectors_Tensor{ii,kk},...
            AllFibers.Orientation_Tensor{ii,kk},...
            AllFibers.px_Tensor{ii,kk},...
            AllFibers.py_Tensor{ii,kk},...
            AllFibers.pz_Tensor{ii,kk},...
            AllFibers.PrincipalAxisLength{ii,kk}]=Vlad_fiber_modelling_v3(full(AllFibers.Object{ii,kk}),AllFibers.BoundingBox{ii,kk},v,p);

    end
%  plot3(AllFibers.px_Tensor{1,1},AllFibers.py_Tensor{1,1},AllFibers.pz_Tensor{1,1})
%  xlabel('x'); ylabel('y'); zlabel('z')
    disp(strcat('Refined track: ',num2str(ii),' from total of: ',num2str(nr_fibers)))


    if p.make_movie==1
        close(v)
    end
%ii
end

%%% save data
if ~exist(strcat(p.save,'3_Refined_fibers\'),'dir')
    mkdir(strcat(p.save,'3_Refined_fibers\'))
end
if p.save_data_modelling==1
    % save as one big file
    save(strcat(p.save,'3_Refined_fibers\',"AllFibers.mat"),"AllFibers","p")
end

end

%% QUANTITIES COMPUTATION
if p.quantities==1
    if p.modelling==0
        load(strcat(p.save,'3_Refined_fibers\',"AllFibers.mat"),"AllFibers")
    end

clc; clearvars -except p AllFibers

% compute
for ij=1:size(AllFibers.Centroid(:,1),1)
    index_time=find(~cellfun('isempty',AllFibers.Centroid(ij,:)))';

    if p.correct_orientations==1
        for it=index_time'
            switch p.use_which_vectors_for_orientation
                case 'poly'
                    AllFibers.red_tensor{ij,it}=[AllFibers.EigenVectors_Tensor{ij,it}(1,1),AllFibers.EigenVectors_Tensor{ij,it}(1,2),AllFibers.EigenVectors_Tensor{ij,it}(1,3)];
                    AllFibers.green_tensor{ij,it}=[AllFibers.EigenVectors_Tensor{ij,it}(2,1),AllFibers.EigenVectors_Tensor{ij,it}(2,2),AllFibers.EigenVectors_Tensor{ij,it}(2,3)];
                    AllFibers.blue_tensor{ij,it}=[AllFibers.EigenVectors_Tensor{ij,it}(3,1),AllFibers.EigenVectors_Tensor{ij,it}(3,2),AllFibers.EigenVectors_Tensor{ij,it}(3,3)];
                case 'region'
                    AllFibers.red_tensor{ij,it}=[AllFibers.EigenVectors{ij,it}{1,1}(1,1),AllFibers.EigenVectors{ij,it}{1,1}(2,1),AllFibers.EigenVectors{ij,it}{1,1}(3,1)];
                    AllFibers.green_tensor{ij,it}=[AllFibers.EigenVectors{ij,it}{1,1}(1,2),AllFibers.EigenVectors{ij,it}{1,1}(2,2),AllFibers.EigenVectors{ij,it}{1,1}(3,2)];
                    AllFibers.blue_tensor{ij,it}=[AllFibers.EigenVectors{ij,it}{1,1}(1,3),AllFibers.EigenVectors{ij,it}{1,1}(2,3),AllFibers.EigenVectors{ij,it}{1,1}(3,3)];
                case 'svd'
                    PrincAxis = Beppe_PricipalAxis(full(AllFibers.Object{ij,it}));
                    AllFibers.red_tensor{ij,it}=PrincAxis(:,1);
                    AllFibers.green_tensor{ij,it}=PrincAxis(:,2);
                    AllFibers.blue_tensor{ij,it}=PrincAxis(:,3);
            end
        end

        %%% set the red vector in the x direction for the first time-step
        if abs(Vlad_angle_between_two_vectors(AllFibers.red_tensor{ij,index_time(1)},[1 0 0]))>90
            AllFibers.red_tensor{ij,index_time(1)}=-AllFibers.red_tensor{ij,index_time(1)};
        end

        for it=index_time(1:end-1)'

            AllFibers.blue_tensor{ij,it}=cross(AllFibers.red_tensor{ij,it},AllFibers.green_tensor{ij,it});


            if Vlad_angle_between_two_vectors(AllFibers.red_tensor{ij,it},AllFibers.red_tensor{ij,it+1})>p.angle
                AllFibers.red_tensor{ij,it+1}=-AllFibers.red_tensor{ij,it+1};
                if Vlad_angle_between_two_vectors(AllFibers.green_tensor{ij,it},AllFibers.green_tensor{ij,it+1})>p.angle
                    AllFibers.green_tensor{ij,it+1}=-AllFibers.green_tensor{ij,it+1};
                    AllFibers.blue_tensor{ij,it+1}=cross(AllFibers.red_tensor{ij,it+1},AllFibers.green_tensor{ij,it+1});
                end
            elseif Vlad_angle_between_two_vectors(AllFibers.green_tensor{ij,it},AllFibers.green_tensor{ij,it+1})>p.angle
                AllFibers.green_tensor{ij,it+1}=-AllFibers.green_tensor{ij,it+1};
                if Vlad_angle_between_two_vectors(AllFibers.red_tensor{ij,it},AllFibers.red_tensor{ij,it+1})>p.angle
                    AllFibers.red_tensor{ij,it+1}=-AllFibers.red_tensor{ij,it+1};
                    AllFibers.blue_tensor{ij,it+1}=cross(AllFibers.red_tensor{ij,it+1},AllFibers.green_tensor{ij,it+1});
                end
            elseif Vlad_angle_between_two_vectors(AllFibers.blue_tensor{ij,it},AllFibers.blue_tensor{ij,it+1})>p.angle
                AllFibers.blue_tensor{ij,it+1}=-AllFibers.blue_tensor{ij,it+1};
                if Vlad_angle_between_two_vectors(AllFibers.green_tensor{ij,it},AllFibers.green_tensor{ij,it+1})>p.angle
                    AllFibers.green_tensor{ij,it+1}=-AllFibers.green_tensor{ij,it+1};
                    AllFibers.red_tensor{ij,it+1}=-cross(AllFibers.blue_tensor{ij,it+1},AllFibers.green_tensor{ij,it+1});
                end

            end

            AllFibers.Orientation_Tensor{ij,it+1}=Vlad_find_euler_angles(AllFibers.red_tensor{ij,it+1},AllFibers.green_tensor{ij,it+1},AllFibers.blue_tensor{ij,it+1});
        end
        AllFibers.blue_tensor{ij,index_time(end)}=cross(AllFibers.red_tensor{ij,index_time(end)},AllFibers.green_tensor{ij,index_time(end)});

    end

    positions=cell2mat(AllFibers.Centroid_Refined(ij,index_time)') * p.dx; %%% stored as y,x,z
    orientations=cell2mat(AllFibers.Orientation_Tensor(ij,index_time)'); %%% stored as y,x,z
    timesteps=index_time-index_time(1)+1;
    time=(timesteps-1)*p.dt;

    x=positions(:,1);
    y=positions(:,2);
    z=positions(:,3);

    theta=orientations(:,1);
    phi=orientations(:,2);
    psi=orientations(:,3);

    %%%%%%% ----positions and translational velocity---- %%%%%%%%
    [x_fitted,~]=fit_data(time,x,p.fitting_type,p.kernel_trajectory_positions);
    [y_fitted,~]=fit_data(time,y,p.fitting_type,p.kernel_trajectory_positions);
    [z_fitted,~]=fit_data(time,z,p.fitting_type,p.kernel_trajectory_positions);

    x_dot=Vlad_compute_derivative(time,x_fitted,p.derivative_stencil,p.disable_edge_points);
    y_dot=Vlad_compute_derivative(time,y_fitted,p.derivative_stencil,p.disable_edge_points);
    z_dot=Vlad_compute_derivative(time,z_fitted,p.derivative_stencil,p.disable_edge_points);

    x_dot_dot=Vlad_compute_derivative(time,x_dot,p.derivative_stencil,p.disable_edge_points);
    y_dot_dot=Vlad_compute_derivative(time,y_dot,p.derivative_stencil,p.disable_edge_points);
    z_dot_dot=Vlad_compute_derivative(time,z_dot,p.derivative_stencil,p.disable_edge_points);


    %%%%%%% ---- angles and rotational rates---- %%%%%%%%
    if p.use_acos_cos==1
        theta=acosd(abs(cosd(theta)));
        phi=acosd(abs(cosd(phi)));
        psi=acosd(abs(cosd(psi)));
        [theta_fitted,~]=fit_data(time,theta,p.fitting_type,p.kernel_trajectory_angles);
        [phi_fitted,~]=fit_data(time,phi,p.fitting_type,p.kernel_trajectory_angles);
        [psi_fitted,~]=fit_data(time,psi,p.fitting_type,p.kernel_trajectory_angles);

        theta_dot=Vlad_compute_derivative(time,theta_fitted,p.derivative_stencil,p.disable_edge_points);
        phi_dot=Vlad_compute_derivative(time,phi_fitted,p.derivative_stencil,p.disable_edge_points);
        psi_dot=Vlad_compute_derivative(time,psi_fitted,p.derivative_stencil,p.disable_edge_points);
    else

        if p.remove_outliers==1
            theta_filled=filloutliers(theta,"linear","movmedian",p.windows_removal,"SamplePoints",timesteps);
            psi_filled=filloutliers(psi,"linear","movmedian",p.windows_removal,"SamplePoints",timesteps);
            phi_filled=filloutliers(phi,"linear","movmedian",p.windows_removal,"SamplePoints",timesteps);

            [theta_fitted,~]=fit_data(time,theta_filled,p.fitting_type,p.kernel_trajectory_angles);
            [phi_fitted,~]=fit_data(time,phi_filled,p.fitting_type,p.kernel_trajectory_angles);
            [psi_fitted,~]=fit_data(time,psi_filled,p.fitting_type,p.kernel_trajectory_angles);
        else
            [theta_fitted,~]=fit_data(time,theta,p.fitting_type,p.kernel_trajectory_angles);
            [phi_fitted,~]=fit_data(time,phi,p.fitting_type,p.kernel_trajectory_angles);
            [psi_fitted,~]=fit_data(time,psi,p.fitting_type,p.kernel_trajectory_angles);
        end


        theta_dot=Vlad_compute_derivative(time,theta_fitted,p.derivative_stencil,p.disable_edge_points);
        phi_dot=Vlad_compute_derivative(time,phi_fitted,p.derivative_stencil,p.disable_edge_points);
        psi_dot=Vlad_compute_derivative(time,psi_fitted,p.derivative_stencil,p.disable_edge_points);
    end


    %%%%%% ---Tumbling from orientation of principal axis of the fiber - red vector ----- %%%%% 
    red=[]; green=[]; blue=[]; R=[];
    for it=index_time'
        if it==1
            tt=1;
        else
            tt=it-index_time(1)+1;
        end
        % extract the vectors of the fibre
        red(tt,:)=AllFibers.red_tensor{ij,it};
        green(tt,:)=AllFibers.green_tensor{ij,it};
        blue(tt,:)=AllFibers.blue_tensor{ij,it};

    end


    if p.use_fitted_vectors_for_rotation_matrix == 1
        % fit the red vector
        [red_fitted_1,~]=fit_data(time,red(:,1),p.fitting_type,p.kernel_trajectory_vectors);
        [red_fitted_2,~]=fit_data(time,red(:,2),p.fitting_type,p.kernel_trajectory_vectors);
        [red_fitted_3,~]=fit_data(time,red(:,3),p.fitting_type,p.kernel_trajectory_vectors);
        % fit the green vector
        [green_fitted_1,~]=fit_data(time,green(:,1),p.fitting_type,p.kernel_trajectory_vectors);
        [green_fitted_2,~]=fit_data(time,green(:,2),p.fitting_type,p.kernel_trajectory_vectors);
        [green_fitted_3,~]=fit_data(time,green(:,3),p.fitting_type,p.kernel_trajectory_vectors);

        for it=index_time'
                if it==1
                    tt=1;
                else
                    tt=it-index_time(1)+1;
                end
        
                % compute a temporary vector perpendicular to the red and green fitted
                % vectors
                green_temp=cross([red_fitted_1(tt); red_fitted_2(tt); red_fitted_3(tt)],[green_fitted_1(tt); green_fitted_2(tt); green_fitted_3(tt)]);
                green_temp=green_temp/norm(green_temp);
        
                % compute a vector perpendicular to the red and the temporary
                % vector (this is in the same plane as the fitted red and green vectors, but it's perpendicular to the fitted red one)
                green_fitted=cross(green_temp,[red_fitted_1(tt); red_fitted_2(tt); red_fitted_3(tt)]); green_fitted=green_fitted/norm(green_fitted);
                green_fitted_1(tt)=green_fitted(1);
                green_fitted_2(tt)=green_fitted(2);
                green_fitted_3(tt)=green_fitted(3);
            
                blue_fitted=cross([red_fitted_1(tt); red_fitted_2(tt); red_fitted_3(tt)],green_fitted); blue_fitted=blue_fitted/norm(blue_fitted);
                blue_fitted_1(tt)=blue_fitted(1);
                blue_fitted_2(tt)=blue_fitted(2);
                blue_fitted_3(tt)=blue_fitted(3);

                % generate the rotation matrix - each column corresponds to one vector of the fibre
                R(:,1,tt)=[red_fitted_1(tt); red_fitted_2(tt); red_fitted_3(tt)];
                R(:,2,tt)=[green_fitted_1(tt); green_fitted_2(tt); green_fitted_3(tt)];
                R(:,3,tt)=[blue_fitted_1(tt); blue_fitted_2(tt); blue_fitted_3(tt)];
        end
    else             
        for it=index_time'
            if it==1
                tt=1;
            else
                tt=it-index_time(1)+1;
            end
            % generate the rotation matrix - each column corresponds to one vector of the fibre
            R(:,1,tt)=red(tt,:)';
            R(:,2,tt)=green(tt,:)';
            R(:,3,tt)=blue(tt,:)';
        end
    end

    red_dot_1=rad2deg(Vlad_compute_derivative(time,red_fitted_1,p.derivative_stencil,p.disable_edge_points));
    red_dot_2=rad2deg(Vlad_compute_derivative(time,red_fitted_2,p.derivative_stencil,p.disable_edge_points));
    red_dot_3=rad2deg(Vlad_compute_derivative(time,red_fitted_3,p.derivative_stencil,p.disable_edge_points));

     
    %%%%%------%%%%
    % Computation like Mobin
    omega_x=psi_dot + cos(theta_fitted).*phi_dot;
    omega_y=-theta_dot.*sin(psi_fitted) + sin(theta_fitted).*cos(psi_fitted).*phi_dot;
    omega_z=theta_dot.*cos(psi_fitted) + sin(theta_fitted).*sin(psi_fitted).*phi_dot;

    % Computation with quaternions
    eulerAngles = [theta_fitted,phi_fitted,psi_fitted];
    q = quaternion(eulerAngles,'eulerd','ZYX','frame');
    av = rad2deg(angvel(q,p.dt,'frame')); % units in deg/s

    omega_q_x=av(:,1);
    omega_q_y=av(:,3);
    omega_q_z=-av(:,2);

    % Computation with rotation matrix
    % skip some time-steps to lower the effect of noise on the result
    % compute the derivative of the rotation matrix
    R_dot=[]; R_q = []; av2=[];     OMEGA_s_R_dot=[]; omega_s_R_dot=[]; OMEGA_b_R_dot=[]; omega_b_R_dot=[];
    for tt=1:size(R,3)
        av2(tt,:)=[NaN, NaN, NaN];
        omega_b_R_dot(:,tt)=[NaN; NaN; NaN];
        omega_s_R_dot(:,tt)=[NaN; NaN; NaN];
    end
    for it=1:size(R,3)-p.skip_timesteps_rotation_matrix

        R_dot(:,:,it) = (R(:,:,it+p.skip_timesteps_rotation_matrix) - R(:,:,it))/p.skip_timesteps_rotation_matrix;

        R_q(:,:,1) = R(:,:,it);
        R_q(:,:,2) = R(:,:,it+p.skip_timesteps_rotation_matrix);
        q3 = quaternion(R_q,'rotmat','frame');
        av2_temp = rad2deg(angvel(q3,p.skip_timesteps_rotation_matrix,'frame')); % units in deg/frame
        av2(it,:) = av2_temp(2,:);  % disregard the first data point, because it is always wrong -> outlier

    end

    for tt=1:size(R,3)-p.skip_timesteps_rotation_matrix
        OMEGA_s_R_dot(:,:,tt)=rad2deg(R_dot(:,:,tt)*R(:,:,tt)');
    
        omega_s_R_dot(1,tt)=OMEGA_s_R_dot(3,2,tt);
        omega_s_R_dot(2,tt)=OMEGA_s_R_dot(1,3,tt);
        omega_s_R_dot(3,tt)=OMEGA_s_R_dot(2,1,tt);
   
    
        OMEGA_b_R_dot(:,:,tt)=rad2deg(R(:,:,tt)'*R_dot(:,:,tt));
    
        omega_b_R_dot(1,tt)=OMEGA_b_R_dot(3,2,tt);
        omega_b_R_dot(2,tt)=OMEGA_b_R_dot(1,3,tt);
        omega_b_R_dot(3,tt)=OMEGA_b_R_dot(2,1,tt);
        
    end

    omega_q2_x=av2(:,1);
    omega_q2_y=av2(:,3);
    omega_q2_z=-av2(:,2);

    %% write to structure
    for it=index_time'
        if it==1
            tt=1;
        else
            tt=it-index_time(1)+1;
        end

        %%%%%% write positions and euler angles
        AllFibers.x{ij,it}=x(tt);
        AllFibers.y{ij,it}=y(tt);
        AllFibers.z{ij,it}=z(tt);
    
        AllFibers.x_fitted{ij,it}=x_fitted(tt);
        AllFibers.y_fitted{ij,it}=y_fitted(tt);
        AllFibers.z_fitted{ij,it}=z_fitted(tt);
    
        AllFibers.x_dot{ij,it}=x_dot(tt);
        AllFibers.y_dot{ij,it}=y_dot(tt);
        AllFibers.z_dot{ij,it}=z_dot(tt);

        AllFibers.x_dot_dot{ij,it}=x_dot_dot(tt);
        AllFibers.y_dot_dot{ij,it}=y_dot_dot(tt);
        AllFibers.z_dot_dot{ij,it}=z_dot_dot(tt);
    
        AllFibers.theta{ij,it}=theta(tt);
        AllFibers.phi{ij,it}=phi(tt);
        AllFibers.psi{ij,it}=psi(tt);

        AllFibers.theta_dot{ij,it}=theta_dot(tt);
        AllFibers.phi_dot{ij,it}=phi_dot(tt);
        AllFibers.psi_dot{ij,it}=psi_dot(tt);
    
        AllFibers.theta_fitted{ij,it}=theta_fitted(tt);
        AllFibers.phi_fitted{ij,it}=phi_fitted(tt);
        AllFibers.psi_fitted{ij,it}=psi_fitted(tt);
    
        %%%%%% write rotation rates MOBIN
        % components of the angular velocity in body frame
        % coordinates (index b) computed like Mobin
        AllFibers.Omega_b_x_Mobin{ij,it}=omega_x(tt);
        AllFibers.Omega_b_y_Mobin{ij,it}=omega_y(tt);
        AllFibers.Omega_b_z_Mobin{ij,it}=omega_z(tt);

        %%%%%% write rotation rates QUATERNION
        %  components of the angular velocity in body frame
        % coordinates (index b) computed with quaternions
        AllFibers.Omega_b_x_Quaternion{ij,it}=omega_q_x(tt);
        AllFibers.Omega_b_y_Quaternion{ij,it}=omega_q_y(tt);
        AllFibers.Omega_b_z_Quaternion{ij,it}=omega_q_z(tt);
        AllFibers.Spinning_rate_Quaternion{ij,it}=abs(omega_q_x(tt));
        AllFibers.Tumbling_rate_Quaternion{ij,it}=sqrt(omega_q_y(tt)^2 + omega_q_z(tt)^2);
        %AllFibers.Spinning_rate_Quaternion{ij,it}=abs(dot([omega_q_x(tt), omega_q_y(tt), omega_q_z(tt)],AllFibers.red_tensor{ij,it}));
        %AllFibers.Tumbling_rate_Quaternion{ij,it}=norm(cross([omega_q_x(tt), omega_q_y(tt), omega_q_z(tt)],AllFibers.red_tensor{ij,it}));

        %%%%%% write tumbling rate DERIVATIVE of PRINCIPAL AXIS
        AllFibers.Tumbling_rate_pdot{ij,it}=norm([red_dot_1(tt),red_dot_2(tt),red_dot_3(tt)]);


        %%%%% write FITTED ORIENTATION VECTORS
        AllFibers.red_tensor_fitted{ij,it}=[red_fitted_1(tt); red_fitted_2(tt); red_fitted_3(tt)];
        AllFibers.green_tensor_fitted{ij,it}=[green_fitted_1(tt); green_fitted_2(tt); green_fitted_3(tt)];
        AllFibers.blue_tensor_fitted{ij,it}=[blue_fitted_1(tt); blue_fitted_2(tt); blue_fitted_3(tt)];


        %%%%% write rotation rates DERIVATIVE OF ROTATION MATRIX
        AllFibers.Omega_s_x_RotationMatrix{ij,it}=omega_s_R_dot(1,tt);
        AllFibers.Omega_s_y_RotationMatrix{ij,it}=omega_s_R_dot(2,tt);
        AllFibers.Omega_s_z_RotationMatrix{ij,it}=omega_s_R_dot(3,tt);

        AllFibers.Omega_b_x_RotationMatrix{ij,it}=omega_b_R_dot(1,tt);
        AllFibers.Omega_b_y_RotationMatrix{ij,it}=omega_b_R_dot(2,tt);
        AllFibers.Omega_b_z_RotationMatrix{ij,it}=omega_b_R_dot(3,tt);

        AllFibers.Spinning_rate_RotationMatrix{ij,it}=omega_b_R_dot(1,tt);
        AllFibers.Tumbling_rate_RotationMatrix{ij,it}=sqrt(omega_b_R_dot(2,tt)^2 + omega_b_R_dot(3,tt)^2);

        AllFibers.Omega_b_x_Quaternion2{ij,it}=omega_q2_x(tt);
        AllFibers.Omega_b_y_Quaternion2{ij,it}=omega_q2_y(tt);
        AllFibers.Omega_b_z_Quaternion2{ij,it}=omega_q2_z(tt);
        AllFibers.Spinning_rate_Quaternion2{ij,it}=abs(omega_q2_x(tt));
        AllFibers.Tumbling_rate_Quaternion2{ij,it}=sqrt(omega_q2_y(tt)^2 + omega_q2_z(tt)^2);



    end

    disp(strcat('Computed quantities track: ',num2str(ij),' from total of: ',num2str(size(AllFibers.Centroid(:,1),1))))
end

% save data
if ~exist(strcat(p.save,'4_Quantities_Refined_fibers\'),'dir')
    mkdir(strcat(p.save,'4_Quantities_Refined_fibers\'))
end
% % save one big file
% save(strcat(p.save,'4_Quantities_Refined_fibers\',"AllFibers.mat"),"AllFibers","p")

AllFibers.Object=[];

% save one big file of only the data, without the objects of the fibres
save(strcat(p.save,'4_Quantities_Refined_fibers\',"AllFibers_Only_data.mat"),"AllFibers","p")

end

%%% plot figures of tracks
if p.plot==1
    plot_each_track_and_quantities(p,AllFibers)
end



end