function[nr_fibers,Fibers]=Vlad_fiber_discrimination_3D_v4(I,p,kk)
%%% INPUT %%%
% 3D MART object imported in MATLAB thorugh read_im7 function: e.g. I=I(1:nx,1:ny,1:nz)

%%% OUTPUT %%%
% 3D objects of each fiber and its properties such as volume, centroid position etc.

%% Binarization
if p.imbin_thres>0
    BW = imbinarize(I,p.imbin_thres);
else
    BW = imbinarize(I);
end

% identify connected regions and compute their geometrical properties
s = regionprops3(BW,'PrincipalAxisLength','VoxelList','Centroid','Volume','Extent','Image','BoundingBox','Orientation','EigenVectors');
Length=s.PrincipalAxisLength;       % principal axis length in vox of each object
%volume=s.Volume;                   % volume in vox of each object

%% Fiber discretisation
vn2=find(Length(:,1)>p.l_min);      
% consider only objects which respect the conditions such as minimal length or minimal volume or aspect ratio
% e.g. vn2=find(Length(:,1)>p.l_min && volume>1000 && Length(:,1)/Length(:,3)>5); 

vv2=s.VoxelList(vn2,:);             % store all voxels of each object (fiber) which has passed the above discretisation
nr_fibers=numel(vv2);               % total number of fibers in this timestep

%% Asignment
% create a matrix containing only the fiber intensities
II=zeros(size(I,1),size(I,2),size(I,3));
for il=1:numel(vv2)                 
    try1=vv2{il};
    for ii=1:numel(try1(:,2))
        II(try1(ii,2),try1(ii,1),try1(ii,3))=I(try1(ii,2),try1(ii,1),try1(ii,3));
    end
end

% store each fiber property to a unique row in a structure called "Fibers"
i=1;
for ii=vn2'
    
    Vox=s.VoxelList{ii};                        % list of voxels
    min_z=min(Vox(:,3)); max_z=max(Vox(:,3));   % region in which the fiber is contained
    min_x=min(Vox(:,1)); max_x=max(Vox(:,1)); 
    min_y=min(Vox(:,2)); max_y=max(Vox(:,2));

    % store to structure
    Fibers.Object{i,1}=ndSparse(II(min_y:max_y,min_x:max_x,min_z:max_z).*s.Image{ii});       % intensities of the fiber
    Fibers.Centroid{i,1}=s.Centroid(ii,:);                                % location of the centroid of the fiber
    Fibers.BoundingBox{i,1}=[min_x max_x; min_y max_y; min_z max_z];            % bounding reagion
    Fibers.Orientation{i,1}=s.Orientation(ii,:);                                % orientatation
    Fibers.EigenVectors{i,1}=s.EigenVectors(ii,:);                              % Eigenvectors of the equivalent ellipsoid

    %%% rearrange eigenvectors
    Fibers.EigenVectors{i,1}{1,1}=Fibers.EigenVectors{i,1}{1,1}([2 1 3],:); % first row is the x component
                                                                            % second row is the y component
                                                                            % thrid row is the z component

                                                                            % different columns are different vectors

    i=i+1;
end


%% look at the results
if p.plot==1
% %% show the full object and the distretized object
%     %%% show the full object
%     fig1=figure(1); box on; grid on; hold on;
%     fig1.WindowState='maximized';
%     title('Full Object','Interpreter','latex')
%     for i=1:length(p.levellist)
%         level=p.levellist(i);
%         pat=patch(isosurface(I,level));
%         pat.FaceVertexCData=level;
%         pat.FaceColor='flat';
%         pat.EdgeColor='none';
%         pat.FaceAlpha=p.facealphalist(i);
%     end
%     colormap(hsv(numel(p.levellist)-1))
%     colorbar
%     daspect([1 1 1])
%     view(3)
% 
%     xlabel('$x$ (stream-wise)')
%     ylabel('$y$ (span-wise)')
%     zlabel('$z$ (wall-normal)')
% 
%     quiver3(0,0,0,1,0,0,100,'filled','LineWidth',2,'Color',[255 68 59]/255,'Marker','.','MarkerSize',20,'MarkerFaceColor','k','MarkerEdgeColor',[152 152 157]/255,'MaxHeadSize',150);
%     quiver3(0,0,0,0,1,0,100,'filled','LineWidth',2,'Color',[50 215 75]/255,'MaxHeadSize',150)
%     quiver3(0,0,0,0,0,1,100,'filled','LineWidth',2,'Color',[10 132 255]/255,'MaxHeadSize',150)
% 
%     hold off
% 
%     if p.print==1
%         % make folder where to save the figures of the processing
%         if ~exist(strcat(p.save,'Figures_Movies_Processing\1_Binarization\'),'dir')
%             mkdir(strcat(p.save,'Figures_Movies_Processing\1_Binarization\'))
%         end
%         savefig(strcat(p.save,'Figures_Movies_Processing\1_Binarization\','Timestep_',num2str(kk),'.fig'))
%         print(strcat(p.save,'Figures_Movies_Processing\1_Binarization\','Timestep_',num2str(kk),'.tif'),'-dtiffn')
%     end
% 
%     %%% show the binarized object
%     % fig2=figure(2); fig2.Position=[303 47 591 413];
%     % volshow(BW);
% 
%     %%% show the binarized and discretized object
%     %fig3=figure(3); fig3.Position=[895 40 548 419];
%     %volshow(imbinarize(II,0));
% 
%     %%% show the discretized object
%     fig4=figure(4); box on; grid on; hold on;
%     fig4.WindowState='maximized';
%     title('Discretized Object','Interpreter','latex')
%     for i=1:length(p.levellist)
%         level=p.levellist(i);
%         pat=patch(isosurface(II,level));
%         pat.FaceVertexCData=level;
%         pat.FaceColor='flat';
%         pat.EdgeColor='none';
%         pat.FaceAlpha=p.facealphalist(i);
%     end
%     colormap(hsv(numel(p.levellist)-1))
%     colorbar
%     daspect([1 1 1])
%     view(3)
% 
%     xlabel('$x$ (stream-wise)')
%     ylabel('$y$ (span-wise)')
%     zlabel('$z$ (wall-normal)')
% 
%     quiver3(0,0,0,1,0,0,100,'filled','LineWidth',2,'Color',[255 68 59]/255,'Marker','.','MarkerSize',20,'MarkerFaceColor','k','MarkerEdgeColor',[152 152 157]/255,'MaxHeadSize',150);
%     quiver3(0,0,0,0,1,0,100,'filled','LineWidth',2,'Color',[50 215 75]/255,'MaxHeadSize',150)
%     quiver3(0,0,0,0,0,1,100,'filled','LineWidth',2,'Color',[10 132 255]/255,'MaxHeadSize',150)
% 
%     hold off
% 
%     if p.print==1
%         % make folder where to save the figures of the processing
%         if ~exist(strcat(p.save,'Figures_Movies_Processing\2_Discretization\'),'dir')
%             mkdir(strcat(p.save,'Figures_Movies_Processing\2_Discretization\'))
%         end        
%         savefig(strcat(p.save,'Figures_Movies_Processing\2_Discretization\','Timestep_',num2str(kk),'.fig'))
%         print(strcat(p.save,'Figures_Movies_Processing\2_Discretization\','Timestep_',num2str(kk),'.tif'),'-dtiffn')
%     end
% 


%% look at each fiber
    j=1;
    for ii=vn2'
        
        Vox=s.VoxelList{ii};
        min_z=min(Vox(:,3)); max_z=max(Vox(:,3)); min_x=min(Vox(:,2)); max_x=max(Vox(:,2)); min_y=min(Vox(:,1)); max_y=max(Vox(:,1));
    
        fiber=II(min_x:max_x,min_y:max_y,min_z:max_z);
        
        %%% isosurface plotting
        fig6=figure(6); box on; grid on; hold on;
        fig6.WindowState='maximized';
        set(gca,'SortMethod','childorder')
        title(strcat('Fiber, $Vol=',num2str(s.Volume(ii)), '\ L_1=',num2str(s.PrincipalAxisLength(ii,1)),'\ L_2=',num2str(s.PrincipalAxisLength(ii,2)),...
            '\ L_3=',num2str(s.PrincipalAxisLength(ii,3)),'\ Extent=',num2str(s.Extent(ii)),'$'),'Interpreter','latex')
        for i=1:length(p.levellist)
            level=p.levellist(i);
            pat=patch(isosurface(fiber,level));
            %pat=patch(isosurface(s.Image{ii},0));
            pat.FaceVertexCData=level;
            pat.FaceColor='flat';
            pat.EdgeColor='none';
            pat.FaceAlpha=p.facealphalist(i);
        end
        disp(strcat('Volume Fiber =  ',num2str(s.Volume(ii)),' voxel'))
        disp(strcat('Longest Principal Axis Fiber =  ',num2str(s.PrincipalAxisLength(ii,1)),' voxel'))
        disp(strcat('Shortest Principal Axis Fiber =  ',num2str(s.PrincipalAxisLength(ii,3)),' voxel'))
        disp(strcat('Extent Fiber =  ',num2str(s.Extent(ii))))
        daspect([ 1 1 1 ])
        colormap(hsv(numel(p.levellist)-1))
        colorbar
        view(3)

        x=Fibers.Centroid{j,1}(1)-Fibers.BoundingBox{j,1}(1,1);
        y=Fibers.Centroid{j,1}(2)-Fibers.BoundingBox{j,1}(2,1);
        z=Fibers.Centroid{j,1}(3)-Fibers.BoundingBox{j,1}(3,1);

        red=Fibers.EigenVectors{j,1}{1,1}(:,1);
        green=Fibers.EigenVectors{j,1}{1,1}(:,2);
        blue=Fibers.EigenVectors{j,1}{1,1}(:,3);

        quiver3(x,y,z,red(1),red(2),red(3),10,'LineWidth',2,'Color','r','Marker','.','MarkerSize',30,'MarkerFaceColor','k','MarkerEdgeColor','k','MaxHeadSize',20)
        quiver3(x,y,z,green(1),green(2),green(3),10,'LineWidth',2,'Color','g','MaxHeadSize',20)
        quiver3(x,y,z,blue(1),blue(2),blue(3),10,'LineWidth',2,'Color','b','MaxHeadSize',20)

        xlabel('$x$ (stream-wise)')
        ylabel('$y$ (span-wise)')
        zlabel('$z$ (wall-normal)')
    
        quiver3(0,0,0,1,0,0,10,'filled','LineWidth',2,'Color',[255 68 59]/255,'Marker','.','MarkerSize',20,'MarkerFaceColor','k','MarkerEdgeColor',[152 152 157]/255,'MaxHeadSize',30);
        quiver3(0,0,0,0,1,0,10,'filled','LineWidth',2,'Color',[50 215 75]/255,'MaxHeadSize',30)
        quiver3(0,0,0,0,0,1,10,'filled','LineWidth',2,'Color',[10 132 255]/255,'MaxHeadSize',30)
    
        hold off

        if p.print==1
            if ~exist(strcat(p.save,'Figures_Movies_Processing\2_Discretization\Fibers_Timestep_',num2str(kk),'\'),'dir')
                mkdir(strcat(p.save,'Figures_Movies_Processing\2_Discretization\Fibers_Timestep_',num2str(kk),'\'))
            end
            savefig(strcat(p.save,'Figures_Movies_Processing\2_Discretization\Fibers_Timestep_',num2str(kk),'\','Fiber_',num2str(ii),'.fig'))
            print(strcat(p.save,'Figures_Movies_Processing\2_Discretization\Fibers_Timestep_',num2str(kk),'\','Fiber_',num2str(ii),'.tif'),'-dtiffn')
        end
    
    
        if p.pause_enabled==1
            pause
        end
    
        gcf; clf; clc;

        j=j+1;
    
        %%% voxels plotting
    %     fig5=figure(5); fig5.Position=[1 1 600 600];
    %     intensities=sort(nonzeros(fiber));
    %     int=linspace(Input.imbin_thres,max(intensities),numel(intensities))';
    %     
    %     for ij=1:numel(Vox(:,1))
    %         plotcube([1 1 1],[Vox(ij,2)-0.5 Vox(ij,1)-0.5 Vox(ij,3)-0.5],int(ij)/max(intensities),[0 0 0]);
    %     end
    
    end
end
end