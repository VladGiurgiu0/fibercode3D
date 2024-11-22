function plot_each_track_and_quantities(p,AllFibers)
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
%%% plots each track and the quantities like position and velocity


plot_isosurface=0;

for ij=1:size(AllFibers.Centroid(:,1),1)
    % initalize
    timesteps=find(~cellfun('isempty',AllFibers.Centroid(ij,:)')==1);
    time=timesteps*p.dt;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------- figure position, velocity, orientation ------- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1); set(gcf,"Position",[400 50 1000 800]); clf
    color=rand(1,3);

%%% positions, velocities
    subplot(2,3,1);hold on; box on; grid on; 
    plot(time,cell2mat(AllFibers.x(ij,:)),'o','Color',color,'MarkerSize',5,'LineWidth',2)
    plot(time,cell2mat(AllFibers.x_fitted(ij,:)),'k-')
    xlabel('Time [s]','Interpreter','latex'); ylabel('$x$ [mm]','Interpreter','latex');
    legend('raw','filtered')

    subplot(2,3,2);hold on; box on; grid on;
    plot(time,cell2mat(AllFibers.y(ij,:)),'o','Color',color,'MarkerSize',5,'LineWidth',2)
    plot(time,cell2mat(AllFibers.y_fitted(ij,:)),'k-')
    xlabel('Time [s]','Interpreter','latex'); ylabel('$y$ [mm]','Interpreter','latex');
    legend('raw','filtered')

    subplot(2,3,3);hold on; box on; grid on;
    plot(time,cell2mat(AllFibers.z(ij,:)),'o','Color',color,'MarkerSize',5,'LineWidth',2)
    plot(time,cell2mat(AllFibers.z_fitted(ij,:)),'k-')
    xlabel('Time [s]','Interpreter','latex'); ylabel('$z$ [mm]','Interpreter','latex');
    legend('raw','filtered')
%
    subplot(2,3,4); hold on; box on; grid on;
    plot(time,cell2mat(AllFibers.x_dot(ij,:)),'k-')
    xlabel('Time in s','Interpreter','latex'); ylabel('$\dot{x}$ [mm/s]','Interpreter','latex');

    subplot(2,3,5); hold on; box on; grid on;
    plot(time,cell2mat(AllFibers.y_dot(ij,:)),'k-')
    xlabel('Time in s','Interpreter','latex'); ylabel('$\dot{y}$ [mm/s]','Interpreter','latex');

    subplot(2,3,6);hold on; box on; grid on;
    plot(time,cell2mat(AllFibers.z_dot(ij,:)),'k-')
    xlabel('Time in s','Interpreter','latex'); ylabel('$\dot{z}$ [mm/s]','Interpreter','latex');


    if p.print==1
        if ~exist(strcat(p.save,'Figures_Movies_Processing\5_Quantities\'),'dir')
            mkdir(strcat(p.save,'Figures_Movies_Processing\5_Quantities\'))
        end     
        savefig(strcat(p.save,'Figures_Movies_Processing\5_Quantities\','Track_',num2str(ij),'_A','.fig'))
        print(strcat(p.save,'Figures_Movies_Processing\5_Quantities\','Track_',num2str(ij),'_A','.tif'),'-dtiffn')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------- figure orientation, rotation rates ------- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    it=find(~cellfun('isempty',AllFibers.Centroid(ij,:)))';

    for tt=it'
        red(ij,tt,:)=cell2mat(AllFibers.red_tensor(ij,tt));         % red vector - lowest moment of inertia
        green(ij,tt,:)=cell2mat(AllFibers.green_tensor(ij,tt));     % green vector - medium moment of inertia
        blue(ij,tt,:)=cell2mat(AllFibers.blue_tensor(ij,tt));       % blue vector - highest moment of inertia
    
        red_fitted(ij,tt,:)=cell2mat(AllFibers.red_tensor_fitted(ij,tt));         % red vector - lowest moment of inertia
        green_fitted(ij,tt,:)=cell2mat(AllFibers.green_tensor_fitted(ij,tt));     % green vector - medium moment of inertia
        blue_fitted(ij,tt,:)=cell2mat(AllFibers.blue_tensor_fitted(ij,tt));       % blue vector - highest moment of inertia
    
    end

    figure(2); set(gcf,"Position",[800 50 1000 800]); clf
    subplot(2,3,1); hold on; grid on; box on;

    plot(red(ij,it,1),'r.')
    plot(red(ij,it,2),'g.')
    plot(red(ij,it,3),'b.')
    plot(red_fitted(ij,it,1),'r-')
    plot(red_fitted(ij,it,2),'g-')
    plot(red_fitted(ij,it,3),'b-')
    ylabel('$e_{11},\ e_{12},\ e_{13}$')
    xlabel('time-step')
    p1=plot([],[],'k.');
    p2=plot([],[],'k-');
    legend([p1,p2],{'raw','fitted'},'location','best')

    subplot(2,3,2); hold on; grid on; box on;
    plot(green(ij,it,1),'r.')
    plot(green(ij,it,2),'g.')
    plot(green(ij,it,3),'b.')
    plot(green_fitted(ij,it,1),'r-')
    plot(green_fitted(ij,it,2),'g-')
    plot(green_fitted(ij,it,3),'b-')
    ylabel('$e_{21},\ e_{22},\ e_{23}$')
    xlabel('time-step')
    p1=plot([],[],'k.');
    p2=plot([],[],'k-');
    legend([p1,p2],{'raw','fitted'},'location','best')

    subplot(2,3,3); hold on; grid on; box on;
    plot(blue(ij,it,1),'r.')
    plot(blue(ij,it,2),'g.')
    plot(blue(ij,it,3),'b.')
    plot(blue_fitted(ij,it,1),'r-')
    plot(blue_fitted(ij,it,2),'g-')
    plot(blue_fitted(ij,it,3),'b-')
    ylabel('$e_{31},\ e_{32},\ e_{33}$')
    xlabel('time-step')
    p1=plot([],[],'k.');
    p2=plot([],[],'k-');
    legend([p1,p2],{'raw','fitted'},'location','best')

    subplot(2,3,4)
    hold on; box on; grid on;
    pr=plot(time,abs(cell2mat(AllFibers.Omega_b_x_RotationMatrix(ij,:))),'r.-');
    pq=plot(time(2:end),abs(cell2mat(AllFibers.Omega_b_x_Quaternion2(ij,2:end))),'k.-');
    xlabel('Time [s]','Interpreter','latex'); ylabel('$\omega_1$ in deg/s','Interpreter','latex');

    subplot(2,3,5);hold on; box on; grid on;
    plot(time,abs(cell2mat(AllFibers.Omega_b_y_RotationMatrix(ij,:))),'r.-')
    plot(time(2:end),abs(cell2mat(AllFibers.Omega_b_z_Quaternion2(ij,2:end))),'k.-')
    xlabel('Time [s]','Interpreter','latex'); ylabel('$\omega_2$ in deg/s','Interpreter','latex');

    subplot(2,3,6);hold on; box on; grid on;
    plot(time,abs(cell2mat(AllFibers.Omega_b_z_RotationMatrix(ij,:))),'r.-')
    plot(time(2:end),abs(cell2mat(AllFibers.Omega_b_y_Quaternion2(ij,2:end))),'k.-')
    xlabel('Time [s]','Interpreter','latex'); ylabel('$\omega_3$ in deg/s','Interpreter','latex');

    legend([pr,pq],{'Rotation matrix','Quaternion'},'location','best')


    if p.print==1
        if ~exist(strcat(p.save,'Figures_Movies_Processing\5_Quantities\'),'dir')
            mkdir(strcat(p.save,'Figures_Movies_Processing\5_Quantities\'))
        end  
        savefig(strcat(p.save,'Figures_Movies_Processing\5_Quantities\','Track_',num2str(ij),'_B','.fig'))
        print(strcat(p.save,'Figures_Movies_Processing\5_Quantities\','Track_',num2str(ij),'_B','.tif'),'-dtiffn')
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------- figure fiber track, eigenvectors ------- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(3); hold on; box on; grid on;
    fig=gcf; fig.Position=[1 601 1800 400];

    if plot_isosurface==1
    
    max_x=max(cell2mat(AllFibers.x(ij,:))/p.dx);
    max_y=max(cell2mat(AllFibers.y(ij,:))/p.dx);
    max_z=max(cell2mat(AllFibers.z(ij,:))/p.dx);

    II=zeros(ceil(max_x)+100,ceil(max_y)+100,ceil(max_z)+100);
    pos=1;
    for it=timesteps(1:p.skip:end)'
        bounds=AllFibers.BoundingBox{ij,it};
        II(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2))=II(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2))+permute(full(AllFibers.Object{ij,it}),[2 1 3]);
        %II(bounds(2,1):bounds(2,2),bounds(1,1):bounds(1,2),bounds(3,1):bounds(3,2))=II(bounds(2,1):bounds(2,2),bounds(1,1):bounds(1,2),bounds(3,1):bounds(3,2))+full(AllFibers.Object{ij,it});
        
        pos=pos+1;
    end

    for i=1:length(p.levellist)
        level=p.levellist(i);
        pat=patch(isosurface(permute(II,[2 1 3]),level));
        pat.FaceVertexCData=level;
        pat.FaceColor='flat';
        pat.EdgeColor='none';
        pat.FaceAlpha=p.facealphalist(i);
    end
    colormap(hsv(numel(p.levellist)))         
    hold on

    end

    for it=timesteps(1:p.skip:end)'

        r=AllFibers.red_tensor{ij,it};
        g=AllFibers.green_tensor{ij,it};
        b=AllFibers.blue_tensor{ij,it};

        quiver3(AllFibers.Centroid_Refined{ij,it}(1),AllFibers.Centroid_Refined{ij,it}(2),AllFibers.Centroid_Refined{ij,it}(3),...
            r(1),r(2),r(3),30,'LineWidth',2,'Color','r')
        
        quiver3(AllFibers.Centroid_Refined{ij,it}(1),AllFibers.Centroid_Refined{ij,it}(2),AllFibers.Centroid_Refined{ij,it}(3),...
            g(1),g(2),g(3),30,'LineWidth',2,'Color','g')
        
        quiver3(AllFibers.Centroid_Refined{ij,it}(1),AllFibers.Centroid_Refined{ij,it}(2),AllFibers.Centroid_Refined{ij,it}(3),...
            b(1),b(2),b(3),30,'LineWidth',2,'Color','b')

        plot3(AllFibers.px_Tensor{ij,it},AllFibers.py_Tensor{ij,it},AllFibers.pz_Tensor{ij,it},'LineWidth',2,'Color','k','LineStyle','-')
    hold on
    end

    xlabel('$x$'); ylabel('$y$'); zlabel('$z$'); 

    view(3)
    daspect([1 1 1])  

    if p.print==1
        if ~exist(strcat(p.save,'Figures_Movies_Processing\5_Quantities\'),'dir')
            mkdir(strcat(p.save,'Figures_Movies_Processing\5_Quantities\'))
        end        
        savefig(strcat(p.save,'Figures_Movies_Processing\5_Quantities\','Track_',num2str(ij),'_C','.fig'))
        print(strcat(p.save,'Figures_Movies_Processing\5_Quantities\','Track_',num2str(ij),'_C','.tif'),'-dtiffn')
    end


    if p.pause_enabled==1
        disp(' Press any key to see next fiber ')
        pause
    end


    close all
end
end