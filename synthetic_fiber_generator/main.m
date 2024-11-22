%% Generates fiber moving through space to test tracking code
% Program to generate a fibre and move it through space.
% The fibre can rotate and translate. 
% - INPUT:
%       Fibre geometry and motion parameters: Constant translational and 
%       rotational velocities.
% - OUTPUT:
%       Results are stored as light intensities to mimick the ouput of
%       MART (one I matrix for every time-step).

clc ; clear ; close all
%addpath('/Users/vlad/Owncloud/Research/_codes/___common')
%addpath('C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad')
addpath('H:\Fiber23May2023b\___CODE_2023_10_19_13_39\_common\')
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%%   ----------------------     General parameters   ----------------------
save_folder = 'rot_z_50\';
plotting=1;         % if 1, plots are enabled
skip=20;             % if 1, plots all fibres; if x, plots every x fibre
method_matrix_generation='Squares';

if ~exist(save_folder,'dir')
    mkdir(save_folder)
end

%%   --------------------    INPUT: Fibre geometry   ----------------------
% Fibres geometry
geo.fl = 49 ;      % fiber length [px]
geo.ft = 1 ;       % fiber half thickness [px]
geo.nst = 80 ;     % number of segments used to approximate the fibre
geo.st = linspace(0,1,geo.nst);   % curv. coordinate
%  ---- Examples of fibres shape ---
% straight fibre (exact cylinder)
% geo.x=0*geo.st ; geo.z=0*geo.st ; geo.y=1*geo.st ;
% nearly straight fibre (good for correct identification of orientation)
% geo.x=0*geo.st ; geo.z=0*geo.st+0.01*geo.st.^2 ; geo.y=1*geo.st ;
% curved fibre aligned with lab reference frame
% geo.x=geo.st.^2 ; geo.z=2*(geo.st-1/2).^2 ; geo.y=0*geo.st ;
% curved fibre aligned with lab reference frame and symmetric
geo.x=geo.st-1/2 ; geo.z=1+(geo.st-1/2).^2 ; geo.y=0*geo.st ;
% curved fibre aligned with lab reference frame and symmetric
%geo.x=geo.st-1/2 ; geo.y=1+(geo.st-1/2).^2 ; geo.z=0*geo.st ;
% curved fibre without specific orientation
%geo.x=geo.st.^2 ; geo.z=2*(geo.st-1/2).^2 ; geo.y=1-1*geo.st.^2 ;

% curved fibre aligned with lab reference frame and symmetric
%geo.x=geo.st/3-1/2 ; geo.y=1+(geo.st-1/2).^2 ; geo.z=0*geo.st ;

%%   --------------------    INPUT: Fibre motion     ----------------------
translation_speed_1 = 0; % px per second in x direction
translation_speed_2 = 0; % px per second in y direction
translation_speed_3 = 0; % px per second in z direction

% Define the rotation speed
rotation_speed_1 = deg2rad(0); % x degrees per second SPINNING 1    (RED)
rotation_speed_2 = deg2rad(50); % x degrees per second TUMBLING 2   (GREEN)
rotation_speed_3 = deg2rad(0); % x degrees per second TUMBLING 3    (BLUE)

% Define the rotation duration
rotation_duration = 1; % rotate for x seconds

% Define the time-steps
dt = 0.01;          % time-step in seconds
num_substeps=1;    % number of sub-time-steps to apply smaller increments of rotation

%%   --------------------    INPUT: Graphics     ----------------------
graph.filename = 'fibre_1' ;
% Graphical options
graph.nm = 10 ;         % approx. number of bullets to show on the fibre
graph.sa = geo.fl/4 ;   % set the legth of the ref. frames (axis scale)
graph.lwr = 4 ;         % thickness of reference frame vectors
graph.lws = 3 ;         % thickness of fobres skeleton
graph.ms = 8 ;          % size of markers on the fibres
% fibre color
% graph.mycmap = flip(cool(geo.nst)) ;
% graph.mycmap = flip(gray(geo.nst)) ;
% graph.mycmap(:,2) = graph.mycmap(:,2)*0.5 ;
graph.mycmap = summer(geo.nst) ;
graph.mygreen = [0 0.8 0] ; % my dark green
graph.plot =  1 ;       % plot one every graph.plot points

%%   -----------    GENERATE: Fibre orientation over time   ---------------
%%%%%% find the initial orientation of the fibre ---- TO BE DONE instead of
% length of fibre segments
geo.d = sum(sqrt(...
    (geo.x(1:end-1)-geo.x(2:end)).^2+...
    (geo.y(1:end-1)-geo.y(2:end)).^2+...
    (geo.z(1:end-1)-geo.z(2:end)).^2));
% normalized by desired length
geo.x = geo.x/geo.d*geo.fl ; geo.y = geo.y/geo.d*geo.fl ; geo.z = geo.z/geo.d*geo.fl ;
% translate away from origin
geo.x = geo.x + geo.fl + 5 ; geo.y = geo.y + geo.fl + 5 ; geo.z = geo.z + geo.fl + 5 ;

% Calculate the number of time steps
num_steps = round(rotation_duration / dt);

% Initial fibre position and orientation
% Find the center of mass
dyn.nt= num_steps ;
X = zeros(geo.nst,dyn.nt) ; Y = X ; Z = X ;
X(:,1) = geo.x(:) ; Y(:,1) = geo.y(:) ; Z(:,1) = geo.z(:) ;
dyn.CM = mean([X(:,1) Y(:,1) Z(:,1)]) ; % center of mass
[dyn.V1,dyn.V2,dyn.V3] = find_fibre_frame(X(:,1),Y(:,1),Z(:,1),dyn.CM) ;

% Define the initial orientation of the body-fixed coordinate system
% rows correspond to the vectors
initial_orientation = [dyn.V1'; dyn.V2'; dyn.V3'];
% initial_orientation = [1 0 0; 0 1 0; 0 0 1];
% initial_orientation = axang2rotm([0 0 1 pi/4])*initial_orientation;

% Initial axes
axis_1=initial_orientation(1,:);  % generates spinning (e.g. rotation around red vector)
axis_2=initial_orientation(2,:);  % generates tumbling 2 (e.g. rotation around green vector)
axis_3=initial_orientation(3,:);  % generates tumbling 3 (e.g. rotation around blue vector)

% Normalize the initial axes
axis_1=axis_1/norm(axis_1); axis_2=axis_2/norm(axis_2); axis_3=axis_3/norm(axis_3);

% Initialize the body fixed vectors array
body_fixed_vectors = zeros(3, 3, num_steps);

body_fixed_vectors(:,1,1)=initial_orientation(1,:)';
body_fixed_vectors(:,2,1)=initial_orientation(2,:)';
body_fixed_vectors(:,3,1)=initial_orientation(3,:)';

% Loop through each time step and generate the rotation
for i = 2:num_steps
    % Calculate the rotation angle for this time step
    rotation_angle_1 = rotation_speed_1 * i * dt;
    rotation_angle_2 = rotation_speed_2 * i * dt;
    rotation_angle_3 = rotation_speed_3 * i * dt;

    % Increment the rotation by substeps
    for jj=1:num_substeps*i
        % Compute incremental rotation matrix for rotations around each
        % axis
        rotation_matrix_1=axang2rotm([axis_1 rotation_angle_1/(num_substeps*i)]);
        rotation_matrix_2=axang2rotm([axis_2 rotation_angle_2/(num_substeps*i)]);
        rotation_matrix_3=axang2rotm([axis_3 rotation_angle_3/(num_substeps*i)]);

        % Compute the complete rotation matrix as a product of the
        % rotations around axis 1, 2, and 3. The order of this product
        % matters only if the applied rotation is too large if
        % rotations around more than one axis is desired
        if jj==1
            rotation_matrix=rotation_matrix_1 * rotation_matrix_2 * rotation_matrix_3;
        else
            rotation_matrix=rotation_matrix_1 * rotation_matrix_2 * rotation_matrix_3 * rotation_matrix;
        end

        % Calculate the body-fixed coordinate system for this time step
        body_fixed_vectors(:, :, i) = rotation_matrix * initial_orientation';

        % Update the new axis of rotation 
        axis_1=body_fixed_vectors(:,1,i)'; axis_1=axis_1/norm(axis_1);
        axis_2=body_fixed_vectors(:,2,i)'; axis_2=axis_2/norm(axis_2);
        axis_3=body_fixed_vectors(:,3,i)'; axis_3=axis_3/norm(axis_3);

    end
end

% Compute the orientation matrix in each time-step (equivalent to the
% rotation matrix with respect to the lab ref. frame)
for i = 1:num_steps
    R(:,:,i)=[body_fixed_vectors(:,1,i) body_fixed_vectors(:,2,i) body_fixed_vectors(:,3,i)];
end

% Loop through each time step and generate the translation
for it=2:num_steps
    % define translations
    X(:,it) = X(:,1) + translation_speed_1*(it-1)*dt ; % translation in x
    Y(:,it) = Y(:,1) + translation_speed_2*(it-1)*dt ; % translation in y
    Z(:,it) = Z(:,1) + translation_speed_3*(it-1)*dt ; % translation in z

    % Determine coordinates of the center of mass (CM) and perform rotation
    dyn.P0 = [X(:,it) Y(:,it) Z(:,it)] ; dyn.CM = mean(dyn.P0) ;
    dyn.P1 = R(:,:,it)*(dyn.P0(:,[1 3 2])-dyn.CM(:,[1 3 2]))' ; dyn.P1 = dyn.CM + dyn.P1';
    X(:,it) = dyn.P1(:,1) ; Y(:,it) = dyn.P1(:,2) ; Z(:,it) = dyn.P1(:,3) ;

%     % Determine coordinates of the center of mass (CM) and perform rotation
%     dyn.P0 = [X(:,it) Y(:,it) Z(:,it)] ; dyn.CM = mean(dyn.P0) ;
%     dyn.P1 = R(:,:,it)*(dyn.P0-dyn.CM)' ; dyn.P1 = dyn.CM + dyn.P1';
%     X(:,it) = dyn.P1(:,1) ; Y(:,it) = dyn.P1(:,2) ; Z(:,it) = dyn.P1(:,3) ;
end

%%   -----------    COMPUTE: Fibre angular velocity   ---------------
% Find the change of the orientation
R_dot=diff(R,1,3)/dt; 

for i = 3:num_steps
    OMEGA_s(:,:,i)=rad2deg(R_dot(:,:,i-1)*R(:,:,i)');

    omega_s(1,i)=OMEGA_s(3,2,i);
    omega_s(2,i)=OMEGA_s(1,3,i);
    omega_s(3,i)=OMEGA_s(2,1,i);


    % OMEGA_s is build like this: [0          -omega_s_3        omega_s_2 ; 
    %                           omega_s_3         0           -omega_s_1 ; 
    %                           -omega_s_2        omega_s_1         0]
    % where omega_s_1 means the rotation around the axis 1 of the lab
    % reference frame

    OMEGA_b(:,:,i)=rad2deg(R(:,:,i)'*R_dot(:,:,i-1));

    omega_b(1,i)=OMEGA_b(3,2,i);
    omega_b(2,i)=OMEGA_b(1,3,i);
    omega_b(3,i)=OMEGA_b(2,1,i);
    
    % OMEGA_b is build like this: [0          -omega_b_3        omega_b_2 ; 
    %                           omega_b_3         0           -omega_b_1 ; 
    %                           -omega_b_2        omega_b_1         0]
    % where omega_b_1 means the rotation around the axis 1 of the body
    % reference frame

end

%%   ----    COMPUTE: Fibre angular velocity and translation with KABSCH algorithm   -------
% for it = 1:num_steps-1
%     P = [X(:,it) Y(:,it) Z(:,it)]';
%     Q = [X(:,it+1) Y(:,it+1) Z(:,it+1)]';
% 
%     [U,r,lrms] = Kabsch(P,Q);
% 
%     U_kabsch(:,:,it) = U; 
%     R_kabsch(:,:,it+1) = U_kabsch(:,:,it) * R(:,:,it); 
% end


%%   -----------    GENERATE: Light intensity matrix   ---------------
pix.X = ceil(X)  ; pix.Y = ceil(Y)  ; pix.Z = ceil(Z)  ;
pix.minX=min(min(pix.X)) ; pix.maxX=max(max(pix.X)) ;
pix.minY=min(min(pix.Y)) ; pix.maxY=max(max(pix.Y)) ;
pix.minZ=min(min(pix.Z)) ; pix.maxZ=max(max(pix.Z)) ;

for it = 1:num_steps
    % initialize the light intensity matrix
    I=zeros(pix.maxY+5*geo.ft,pix.maxX+5*geo.ft,pix.maxZ+5*geo.ft);

    switch method_matrix_generation
        case 'Squares'
            for i=1:size(pix.X,1)
                I(pix.Y(i,it)-geo.ft:pix.Y(i,it)+geo.ft,...
                pix.X(i,it)-geo.ft:pix.X(i,it)+geo.ft,...
                pix.Z(i,it)-geo.ft:pix.Z(i,it)+geo.ft) = 1;
            end

        case 'Spheres'
            %I(:,:,:,it)=generate_fiber_from_polynomial(1:size(I,2),1:size(I,1),1:size(I,3),pix.X(:,it),pix.Y(:,it),pix.Z(:,it),geo.ft);
            I(:,:,:)=generate_fiber_from_polynomial(1:size(I,2),1:size(I,1),1:size(I,3),X(:,it),Y(:,it),Z(:,it),geo.ft);

    end


    % Gaussian filtering
    %I(:,:,:,it)=imgaussfilt3(I(:,:,:,it),1);

    % Compress the matrix to save space
    I=ndSparse(I);

    % Save each time-step
    save(strcat(save_folder,'I_',num2str(it),'.mat'),'I')
    clear I

    disp(strcat('Timesteps completed: ',num2str(it)))
end

% save parameters
time = (1:num_steps)*dt;
input_spinning_1 = ones(1,num_steps) * rad2deg(rotation_speed_1);
input_tumbling_2 = ones(1,num_steps) * rad2deg(rotation_speed_2);
input_tumbling_3 = ones(1,num_steps) * rad2deg(rotation_speed_3);
geo.curvature = Vlad_compute_curvature_length(geo.x',geo.y',geo.z');% average curvature [1/vox]
geo.curvature_ref = pi/geo.fl;                                      % reference curvature i.e. if the fibre is placed in a circular shape [1/vox]
geo.length_prime = abs(max(geo.z)-min(geo.z));                      % fibre extension perpendicular to the main axis [vox]
geo.major_axis = abs(max(geo.x)-min(geo.x));                        % fibre extension in the direction of the main axis [vox]

save(strcat(save_folder,'_parameters.mat'),"time","input_tumbling_3","input_tumbling_2","input_spinning_1","geo")

%%   ----    COMPUTE: Fibre angular velocity and translation with ICP algorithm   -------
% jump_steps=10;
% for it = 1:jump_steps:num_steps-jump_steps
% 
%     %%%%% generate the point clouds
%     %%% ---- fixed ---- %%%
%     % find the positions of light intensities
%     [b,a,c]=ind2sub(size(I(:,:,:,it)),find(I(:,:,:,it)==1));
%     
%     % asign the light intensities
%     for ii=1:numel(b)
%         intensity(ii)=full(I(b(ii),a(ii),c(ii),it));
%     end
%     
%     % generate the fixed point cloud of blue color
%     %fixed = pointCloud([a,b,c],'Color',[zeros(numel(intensity),1) zeros(numel(intensity),1) 255*intensity'],"Intensity",intensity');
%     fixed = pointCloud([a,b,c]);
% 
%     clear b a c intensity
% 
%     %%% ---- moving ---- %%%
%     % find the positions of light intensities
%     [b,a,c]=ind2sub(size(I(:,:,:,it+jump_steps)),find(I(:,:,:,it+jump_steps)==1));
%     
%     % asign the light intensities
%     for ii=1:numel(b)
%         intensity(ii)=full(I(b(ii),a(ii),c(ii),it+jump_steps));
%     end
%     
%     % generate the fixed point cloud of blue color
%     %moving = pointCloud([a,b,c],'Color',[255*intensity' zeros(numel(intensity),1) zeros(numel(intensity),1)],"Intensity",intensity');  
%     moving = pointCloud([a,b,c]);
% 
%     % show the pair of the fixed and moving point clouds
%     %pcshowpair(moving,fixed,VerticalAxis='Y',VerticalAxisDir='Down')
% 
%     %%%%% downsample the point clouds
%     grid_step=2;    % increase these to downsample more
%     fixed_downsampled = pcdownsample(fixed,'gridAverage',grid_step);
%     moving_downsampled = pcdownsample(moving,'gridAverage',grid_step);
% 
%     %initial_transform = rigidtform3d([0 0 0],[translation_speed_1*dt*jump_steps translation_speed_2*dt*jump_steps translation_speed_3*dt*jump_steps]);
% 
% %     [tform, movingReg,rmse] = pcregistericp(moving_downsampled,fixed_downsampled,...
% %         "Metric","planeToPlane",'MaxIterations',1000,'Tolerance',[1e-10 1e-10],'Verbose',true,...
% %         'InitialTransform',initial_transform);
% 
%     [tform, movingReg,rmse] = pcregistericp(moving_downsampled,fixed_downsampled);
%     % show the pair of the fixed and moving registered point clouds
%     %pcshowpair(movingReg,fixed,VerticalAxis='Y',VerticalAxisDir='Down')
% 
%     clear b a c intensity
% end

%%   -----------    PLOT   ---------------

if plotting==1
    % Plot the rotated body fixed reference frames
%     figure(1);hold on; box on; grid on; daspect([1 1 1]);
%     xlabel('x');ylabel('y');zlabel('z');
%     quiver3(0,0,0,initial_orientation(1,1),initial_orientation(1,2),initial_orientation(1,3),'Color',[0 0 0],'Linewidth',2)
%     quiver3(0,0,0,initial_orientation(2,1),initial_orientation(2,2),initial_orientation(2,3),'Color',[0.4 0.4 0.4],'Linewidth',2)
%      quiver3(0,0,0,initial_orientation(3,1),initial_orientation(3,2),initial_orientation(3,3),'Color',[0.8 0.8 0.8],'Linewidth',2)
%     for i = 1:num_steps
%         quiver3(i, 0, 0, body_fixed_vectors(1, 1, i), body_fixed_vectors(2, 1, i), body_fixed_vectors(3, 1, i), 'r','LineWidth',2);
%         quiver3(i, 0, 0, body_fixed_vectors(1, 2, i), body_fixed_vectors(2, 2, i), body_fixed_vectors(3, 2, i), 'g','LineWidth',2);
%         quiver3(i, 0, 0, body_fixed_vectors(1, 3, i), body_fixed_vectors(2, 3, i), body_fixed_vectors(3, 3, i), 'b','LineWidth',2);
%     end
%     view(3)

    % Plot the results of spinning and tumbling
    figure(2); hold on; box on; grid minor; %daspect([1 1 1]);

    plot((1:num_steps)*dt,ones(1,num_steps) * rad2deg(rotation_speed_1),'r-','Linewidth',2) % SPINNING 1    (RED)
    plot((1:num_steps)*dt,ones(1,num_steps) * rad2deg(rotation_speed_2),'g-','Linewidth',2) % TUMBLING 2    (GREEN)
    plot((1:num_steps)*dt,ones(1,num_steps) * rad2deg(rotation_speed_3),'b-','Linewidth',2) % TUMBLING 3    (BLUE)

    plot((1:num_steps)*dt,omega_b(1,:),'ro','MarkerSize',10) % SPINNING 1    (RED)
    plot((1:num_steps)*dt,omega_b(2,:),'go','MarkerSize',10) % TUMBLING 2    (GREEN)
    plot((1:num_steps)*dt,omega_b(3,:),'bo','MarkerSize',10) % TUMBLING 3    (BLUE)


    legend('INPUT: Spinning','INPUT: Tumbling 1', 'INPUT: Tumbling 2',...
            'OUTPUT: Spinning','OUTPUT: Tumbling 1', 'OUTPUT: Tumbling 2')

    xlabel('Time [s]')
    ylabel('Rotation rate [deg/s]')
    set(gca,'FontSize',20)

    % generate and plot fibers
%     figure(401) ; clf ; hold on ; daspect([1 1 1]);
%     plot_fibre_frame(X(:,1),Y(:,1),Z(:,1),dyn.V1,dyn.V2,dyn.V3,graph,geo.nst)

    % display fibres
    figure(402) ; clf ; hold on ; daspect([1 1 1]); view(3);  grid off; box on;
    set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'TickDir','none','FontSize',20,'Color','none');
    %set(gca,'XColor','none','YColor','none','ZColor','none')
    xlabel('$x$'); ylabel('$y$'); zlabel('$z$')

    %view([-27.35 17.18])


    % initial orientation vectors
    %quiver3(mean(X(:,1)), mean(Y(:,1)), mean(Z(:,1)),initial_orientation(1,1),initial_orientation(1,2),initial_orientation(1,3),15,'Color',[0 0 0],'Linewidth',2)
    %quiver3(mean(X(:,1)), mean(Y(:,1)), mean(Z(:,1)),initial_orientation(2,1),initial_orientation(2,2),initial_orientation(2,3),15,'Color',[0.4 0.4 0.4],'Linewidth',2)
    %quiver3(mean(X(:,1)), mean(Y(:,1)), mean(Z(:,1)),initial_orientation(3,1),initial_orientation(3,2),initial_orientation(3,3),15,'Color',[0.8 0.8 0.8],'Linewidth',2)

    IPlot = 0;
    for it = 1:skip:num_steps

        % orientation vectors
        quiver3(mean(X(:,it)), mean(Y(:,it)), mean(Z(:,it)), body_fixed_vectors(1, 1, it), body_fixed_vectors(2, 1, it), body_fixed_vectors(3, 1, it),10, 'r','LineWidth',2,'MaxHeadSize',100);
        quiver3(mean(X(:,it)), mean(Y(:,it)), mean(Z(:,it)), body_fixed_vectors(1, 2, it), body_fixed_vectors(2, 2, it), body_fixed_vectors(3, 2, it),10, 'g','LineWidth',2,'MaxHeadSize',100);
        quiver3(mean(X(:,it)), mean(Y(:,it)), mean(Z(:,it)), body_fixed_vectors(1, 3, it), body_fixed_vectors(2, 3, it), body_fixed_vectors(3, 3, it),10, 'b','LineWidth',2,'MaxHeadSize',100);
        
        % light intensity
        load(strcat(save_folder,'I_',num2str(it),'.mat'),'I')
        IPlot = IPlot + I;
        clear I

        % polynomials
        plot3(X(:,it),Y(:,it),Z(:,it),'-k','Linewidth',2)
    end

    % light intensity
    level=0.2;
    pat=patch(isosurface(IPlot,level));
    pat.FaceVertexCData=level;
    pat.FaceColor='flat'; pat.EdgeColor='none'; pat.FaceAlpha=0.3;

end

%%   ------------------     FUNCTIONS    ------------------
function Fiber_poly=generate_fiber_from_polynomial(X,Y,Z,px,py,pz,sphere_radius)
    %%% input
    % px, py, pz are the centers of each section of the fiber
    % sphere_radius is the radius of the generated sphere at the location of
    % each section of the fibre
    
    %X=linspace(floor(min(px))-sphere_radius,ceil(max(px))+sphere_radius,size(px,1));
    %Y=linspace(floor(min(py))-sphere_radius,ceil(max(py))+sphere_radius,size(py,1));
    %Z=linspace(floor(min(pz))-sphere_radius,ceil(max(pz))+sphere_radius,size(pz,1));
    
    [Xs,Ys,Zs]=meshgrid(X,Y,Z);
    
    %  Salami generator
    % produce first sphere
    Fiber_poly = (Xs-px(1,1)).^2 + (Ys-py(1,1)).^2 + (Zs-pz(1,1)).^2 <=    sphere_radius.^2;
    
    for j=2:size(px,1)          % for each center produce one sphere
        % add them together
        Fiber_poly = or(Fiber_poly,((Xs-px(j,1)).^2 + (Ys-py(j,1)).^2 + (Zs-pz(j,1)).^2 <=    sphere_radius.^2));
    end
end

function plot_fibre_frame(X,Y,Z,V1,V2,V3,graph,nst)

% Function to plot the fibre with its own reference frame.
% INPUT:
%       - X,Y,Z : coordinates of fibre points
%       - V1,V2,V3 : fibres ref. frame vectors
%       - graph : graphical options
% OUTPUT:
%       - Plots fibre and reference frame of the fibre in the figure open

plot3(X(:),Y(:),Z(:),'-k','Linewidth',0.001) % graph.lws
for kk=1:round(nst/graph.nm):nst
    mycol = graph.mycmap(kk,:) ;
    graph.hplot = plot3(X(kk),Y(kk),Z(kk),'o','Color','None',...
        'Markersize',graph.ms,'Markerfacecolor',mycol);
end
% plot reference frame of the fibre
quiver3(mean(X),mean(Y),mean(Z),V1(1),V1(2),V1(3),graph.sa,'-r',...
    'linewidth',graph.lwr) % X
quiver3(mean(X),mean(Y),mean(Z),V2(1),V2(2),V2(3),graph.sa,'-', ...
    'linewidth',graph.lwr,'Color',graph.mygreen) % Y
quiver3(mean(X),mean(Y),mean(Z),V3(1),V3(2),V3(3),graph.sa,'-b', ...
    'linewidth',graph.lwr) % Z

end

function [V1,V2,V3] = find_fibre_frame(X,Y,Z,CM)

% Determine fibre reference frame using find_refframe. This function simply
% prepares the data to be used by find_refframe and detrmines in a more
% clever way the orientation of the eigenvectors found.
% INPUT:
%       X,Y,Z fibres coordinates
% OUTPUT:
%       vectors of the fibre reference frame
% Adapts data to be read by find_refframe

main.xs = X' ; main.ys = Y' ; main.zs = Z' ;
main.debug = 0 ; fp.A = CM ; fp.B = CM ;
ref = find_refframe(main,fp) ;
V1 = ref.V1 ; V2 = ref.V2 ;
% V3 = ref.V3 ; automatically determined from previous two
ns2 = round(numel(main.xs)/2) ; % approximate midpoint
% Extrema points (PE) and midpoint to extrema-midpoint (PM) differences
PE = [...
    main.xs(end)-main.xs(1) , ...
    main.ys(end)-main.ys(1) , ...
    main.zs(end)-main.zs(1)] ;
PM = [main.xs(ns2) , main.ys(ns2) , main.zs(ns2)] - 1/2*[...
    main.xs(end)+main.xs(1) , ...
    main.ys(end)+main.ys(1) , ...
    main.zs(end)+main.zs(1)] ;
% Use dot product to check relative orientation of vector and define sign
% Vector 1/x always from x=0 to x=l of the fibre
if sum(V1'.*PE)<0 ; V1 = - V1 ; end
% Vector 2/y always from inner to outer part of the fibre (same plane of 
% the fibre)
if sum(V2'.*PM)>0 ; V2 = - V2 ; end
% Vector 3/z determined as results of previous two (out of the plane of the
% fibre and perpendicular to previous two, given by the cross product to
% obtain a right-handed frame)
V3 = cross(V1,V2) ;

end


