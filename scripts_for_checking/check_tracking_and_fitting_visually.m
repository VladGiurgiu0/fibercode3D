clc; clear;

Input.f2='C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad_2022_11_07_19_48\_processed_data\Spanwise_Fiber_2022_10_17_permuted\';
load(strcat(Input.f2,"Quantities_Refined_fibers.mat"))

%% show trajectories of all fibers in a period of time
starting=1; % starting timestep
ending=80; % ending timestep

%%% plot tracks to check tracking
fig3=figure(3);
fig3.Position=[201 1 600 600];
hold on
for ij=1:size(AllFibers.Centroid(:,1),1)
    color=rand(1,3);
    for ii=starting:ending-1
        if ~isempty(AllFibers.Centroid{ij,ii})
            scatter3(AllFibers.Centroid{ij,ii}(1),AllFibers.Centroid{ij,ii}(2),AllFibers.Centroid{ij,ii}(3),'filled','MarkerFaceColor',color)
        end
    end
end
box on
grid on

%% show (and make movie) of one fiber trajectory in the laboratory reference system
which_fiber=42;

levellist=[0.1 5];      % which levels of intensity to show
alphalist=[0.2 0.5];

for ij=1:size(AllFibers.Centroid(:,1),1)
    timesteps=find(~cellfun('isempty',AllFibers.Centroid(ij,:)')==1);
    time=timesteps*dt;

    figure();
    set(gcf,'Color','white')
    set(gca,'Color','none', 'Visible','off')
    fig=gcf; fig.Position=[1 41 1440 790];
    xlabel('x')
    ylabel('y')
    zlabel('z')
    colormap(cool(numel(levellist)))

    if ij==which_fiber

        v = VideoWriter(strcat('movies\fiber_track_',num2str(ij),'.mp4'),'MPEG-4');
        v.FrameRate=5;
        v.Quality=100;
        open(v)

        max_x=max(cell2mat(AllFibers.x(ij,:))/dx);
        max_y=max(cell2mat(AllFibers.y(ij,:))/dx);
        max_z=max(cell2mat(AllFibers.z(ij,:))/dx);

        

        pos=1;
        for it=timesteps(1:1:end)'
                   

            II=zeros(ceil(max_x)+100,ceil(max_y)+100,ceil(max_z)+100);
            dim=size(full(AllFibers.Object{ij,it}));
            bounds=AllFibers.BoundingBox{ij,it};
            II(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2))=II(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2))+full(AllFibers.Object{ij,it});
            pos=pos+1;

            r=-AllFibers.red_tensor{ij,it};
            g=-AllFibers.green_tensor{ij,it};
            b=AllFibers.blue_tensor{ij,it};

            hold on

            quiver3(AllFibers.Centroid_Refined{ij,it}(1),AllFibers.Centroid_Refined{ij,it}(2),AllFibers.Centroid_Refined{ij,it}(3),...
                r(1),r(2),r(3),25,'LineWidth',2,'Color','r','MaxHeadSize',150)
            
            quiver3(AllFibers.Centroid_Refined{ij,it}(1),AllFibers.Centroid_Refined{ij,it}(2),AllFibers.Centroid_Refined{ij,it}(3),...
                g(1),g(2),g(3),25,'LineWidth',2,'Color','g','MaxHeadSize',150)
            
            quiver3(AllFibers.Centroid_Refined{ij,it}(1),AllFibers.Centroid_Refined{ij,it}(2),AllFibers.Centroid_Refined{ij,it}(3),...
                b(1),b(2),b(3),25,'LineWidth',2,'Color','b','MaxHeadSize',150)

            plot3(AllFibers.px_Tensor{ij,it},AllFibers.py_Tensor{ij,it},AllFibers.pz_Tensor{ij,it},'LineWidth',2,'Color','k','LineStyle','-')
        
            for i=1:length(levellist)
                level=levellist(i);
                p=patch(isosurface(permute(II,[2 1 3]),level));
                p.FaceVertexCData=level;
                p.FaceColor='flat';
                p.EdgeColor='none';
                p.FaceAlpha=(alphalist(i));
            end

            %II=[];

            daspect([1 1 1])
            set(gca,'Color','none', 'Visible','off')
            view([2 1 2])
            
            ylim([0 570])
            xlim([620 780])
            zlim([500 640])


            frame=getframe(gcf);
            writeVideo(v,frame);

            %pause(1)
            it
            hold off
            clf
        end

    end
    close(v)
    close all
end


%% show each fiber (and make movie) in a reference frame moving with the fiber's center of mass
for ij=1:size(AllFibers.Centroid(:,1),1)
    index_time=find(~cellfun('isempty',AllFibers.Centroid(ij,:)))';
    intensities=reshape(full(AllFibers.Object{ij,1}),[],1);
    %intensities=intensities(intensities>threshold);
    %[~,edges]=find_bin_edges_for_same_counts(intensities,nr_bins_intensity);
    %centers=(edges(1:end-1)+edges(2:end))/2;
    centers=[0.1 6 12];

    levellist=centers;                              % which levels of intensity to show
    alphalist=linspace(0.01,0.6,numel(levellist));  % alpha level correpsonding to the intensity values in levellist


    v = VideoWriter(strcat('movies\fiber_',num2str(ij),'.mp4'),'MPEG-4');
    v.FrameRate=5;
    v.Quality=100;
    open(v)
    for it=index_time'

    Fiber=full(AllFibers.Object{ij,it});

    x=AllFibers.x{ij,it}/p.dx;
    y=AllFibers.y{ij,it}/p.dx;
    z=AllFibers.z{ij,it}/p.dx;

    px=AllFibers.px_Tensor{ij,it};
    py=AllFibers.py_Tensor{ij,it};
    pz=AllFibers.pz_Tensor{ij,it};

    red=AllFibers.red_tensor{ij,it}; 
    green=AllFibers.green_tensor{ij,it}; 
    blue=AllFibers.blue_tensor{ij,it}; 

%     a=a+Bounding_Box(1,1);
%     b=b+Bounding_Box(2,1);
%     c=c+Bounding_Box(3,1);

    [A, B, C]=meshgrid(1:size(Fiber,1),1:size(Fiber,2),1:size(Fiber,3));
    A=A+AllFibers.BoundingBox{ij,it}(1,1); %Bounding_Box(1,1);
    B=B+AllFibers.BoundingBox{ij,it}(2,1); %Bounding_Box(2,1);
    C=C+AllFibers.BoundingBox{ij,it}(3,1); %Bounding_Box(3,1);


    figure(1); colormap((cool(numel(levellist))))
    fig=gcf; clf ;fig.Position=[1 601 1800 400];
    hold all
    
    for i=1:length(levellist)
        level=levellist(i);
        p=patch(isosurface(A,B,C,permute(Fiber,[2 1 3]),level));
        p.FaceVertexCData=level;
        p.FaceColor='flat';
        p.EdgeColor='none';
        p.FaceAlpha=(alphalist(i));
    
    end

    p=patch(isosurface(x+80*ones(size(A)),B,C,permute(Fiber,[2 1 3]),min(centers))); p.FaceVertexCData=level; p.FaceColor='flat'; p.EdgeColor='none'; p.FaceColor=[0.9 0.9 0.9]; p.FaceAlpha=(1);
    p=patch(isosurface(A,y+30*ones(size(B)),C,permute(Fiber,[2 1 3]),min(centers))); p.FaceVertexCData=level; p.FaceColor='flat'; p.EdgeColor='none'; p.FaceColor=[0.9 0.9 0.9]; p.FaceAlpha=(1);
    p=patch(isosurface(A,B,z-51*ones(size(C)),permute(Fiber,[2 1 3]),min(centers))); p.FaceVertexCData=level; p.FaceColor='flat'; p.EdgeColor='none'; p.FaceColor=[0.9 0.9 0.9]; p.FaceAlpha=(1);
    
    p=patch(isosurface(x+79.5*ones(size(A)),B,C,permute(Fiber,[2 1 3]),max(centers))); p.FaceVertexCData=level; p.FaceColor='flat'; p.EdgeColor='none'; p.FaceColor=[0.7 0.7 0.7]; p.FaceAlpha=(1);
    p=patch(isosurface(A,y+29.5*ones(size(B)),C,permute(Fiber,[2 1 3]),max(centers))); p.FaceVertexCData=level; p.FaceColor='flat'; p.EdgeColor='none'; p.FaceColor=[0.7 0.7 0.7]; p.FaceAlpha=(1);
    p=patch(isosurface(A,B,z-50.5*ones(size(C)),permute(Fiber,[2 1 3]),max(centers))); p.FaceVertexCData=level; p.FaceColor='flat'; p.EdgeColor='none'; p.FaceColor=[0.7 0.7 0.7]; p.FaceAlpha=(1);
    
    plot3(px,py,z-50*ones(size(pz)),'-','LineWidth',2,'Color',[255 68 59]/255)
    plot3(px,y+29*ones(size(py)),pz,'-','LineWidth',2,'Color',[50 215 75]/255)
    plot3(x+79*ones(size(px)),py,pz,'-','LineWidth',2,'Color',[10 132 255]/255)
    
    daspect([1 1 1]); box on; %grid on; grid minor;
    
    % scatter3(a,b,c,1,'MarkerFaceColor','r','MarkerEdgeColor','none',...
    %     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    plot3(px,py,pz,'k-','LineWidth',4,'Color',[100 100 100]/255)%[255 55 95]/255)
    %scatter3(px_2,py_2,pz_2,100,'g.')
    %scatter3(b,a,surface,100,'g.')
    %[Xs, Ys, Zs]=meshgrid(1:size(Fiber,1),1:size(Fiber,2),1:size(Fiber,3));
    %salami = isosurface(Xs,Ys,Zs,Fiber_salami,0);
    %patch(salami,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0 0 1],'FaceAlpha',0.5);
    
    quiver3(x,y,z,red(1),red(2),red(3),20,'LineWidth',2,'Color','r','Marker','.','MarkerSize',30,'MarkerFaceColor','k','MarkerEdgeColor','k','MaxHeadSize',150)
    quiver3(x,y,z,green(1),green(2),green(3),20,'LineWidth',2,'Color','g','MaxHeadSize',150)
    quiver3(x,y,z,blue(1),blue(2),blue(3),20,'LineWidth',2,'Color','b','MaxHeadSize',150)
    
    quiver3(x-50,y-50,z-50,1,0,0,15,'filled','LineWidth',2,'Color',[255 68 59]/255,'Marker','.','MarkerSize',20,'MarkerFaceColor','k','MarkerEdgeColor',[152 152 157]/255,'MaxHeadSize',150);
    quiver3(x-50,y-50,z-50,0,1,0,15,'filled','LineWidth',2,'Color',[50 215 75]/255,'MaxHeadSize',150)
    quiver3(x-50,y-50,z-50,0,0,1,15,'filled','LineWidth',2,'Color',[10 132 255]/255,'MaxHeadSize',150)
    view(3)
    set(gcf,"Position",[1 1 600 600])
    xlim([x-51 x+80])
    ylim([y-51 y+30])
    zlim([z-51 z+50])
    
    set(gca,'TickLabelInterpreter','latex','FontSize',15,...
        'LineWidth',2,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'TickLength',[0 0],'FontSize',20)
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('$z$','Interpreter','latex')
    
    clb=colorbar();
    clb.FontSize=15;
    clb.LineWidth=2;
    clb.TickLabelInterpreter="latex";
    clb.Label.String='Counts';
    clb.Label.Interpreter="latex"';
    clb.Location="eastoutside";
    clb.Position=[0.93 0.048 0.036 0.22];
    clb.Limits=[0 12];
    
    %legend('2 counts','7 counts', '13 counts','18 counts','a','b','c','f')
    
    frame=getframe(gcf);
    writeVideo(v,frame);
    
    
    %pause

    end
    close(v)
end