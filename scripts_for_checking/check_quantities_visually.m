clc; clear;
%Input.f2='C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad\_processed_data\Spanwise_Fiber_2022_10_17\run\ImgPreproc_01\TomographicPIV_02\';
%Input.f2='C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad\_processed_data\Spanwise_Fiber_2022_10_26\Fiber_Re720\TimeFilAvg_fL=5_skip_3\TomographicPIV\';
%Input.f2='C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad\_processed_data\Spanwise_Fiber_2022_10_26\Fiber_Re720\Segmented_5\TomographicPIV\';
%Input.f2='C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad\_processed_data\Spanwise_Fiber_2022_10_26\Fiber_Re720\Segmented_no_sharpening_5\';
%Input.f2='C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad\_processed_data\Spanwise_Fiber_2022_10_17\run\ImgPreproc\TomographicPIV\';
%Input.f2='_virtual_data\';
%Input.f2='C:\Users\admin\Desktop\Vlad\_fibercode\Vlad\_processed_data\Spanwise_Fiber_2022_10_17\run\ImgPreproc\TomographicPIV\';
%Input.f2='D:\processed_data_fiber\';
%Input.f2='C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad_2022_11_07_19_48\_processed_data\Spanwise_Fiber_2022_10_17\run\ImgPreproc\TomographicPIV\';
%Input.f2='C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad_2022_11_07_19_48\_processed_data\Spanwise_Fiber_2022_10_17_permuted\';
%Input.f2='C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad_2022_11_07_19_48\_processed_data\Spanwise_Fiber_2022_10_26\Re_720\TimeFilAvg_fL=5_skip_3\TomographicPIV\';
%Input.f2='C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad_2022_11_07_19_48\_processed_data\Spanwise_Fiber_2022_10_17_permuted\';
%Input.f2='C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad_2022_11_07_19_48\_Marco_virtual_fiber_code\_virtual_data\';
Input.f2='C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad_2022_11_07_19_48\_processed_data\Spanwise_Fiber_2022_10_26\Re_720\TimeFilAvg_fL=5_skip_3\TomographicPIV\';


load(strcat(Input.f2,"Quantities_Refined_fibers.mat"))
%load(strcat(Input.f2,"Tracked_fibers.mat"))

plot_tracks=1;
plot_fiber_positions_and_rates=1;
plot_fiber=1;

%% show tracks between timestep
starting=1;
ending=80;

levellist=[0.1 5];      % which levels of intensity to show
alphalist=[0.2 0.5];


%%% plot tracks to check tracking
if plot_tracks==1
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
end

%%% plot track positions to check fitting of trajectories
if plot_fiber_positions_and_rates==1
    for ij=1:size(AllFibers.Centroid(:,1),1)
          %timesteps=find(~cellfun('isempty',AllFibers.Centroid(ij,:)')==1);
        figure(); box on; grid on;
        fig=gcf; fig.Position=[201 1 600 600];
        color=rand(1,3);

        timesteps=find(~cellfun('isempty',AllFibers.Centroid(ij,:)')==1);
        
        time=timesteps*dt;


        %%%%%%% ------- plotting ------- %%%%%%%
        subplot(3,4,1)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.x(ij,:)),'o','Color',color,'MarkerSize',5,'LineWidth',2)
        plot(time,cell2mat(AllFibers.x_fitted(ij,:)),'k-')
        subplot(3,4,5)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.y(ij,:)),'o','Color',color,'MarkerSize',5,'LineWidth',2)
        plot(time,cell2mat(AllFibers.y_fitted(ij,:)),'k-')
        subplot(3,4,9)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.z(ij,:)),'o','Color',color,'MarkerSize',5,'LineWidth',2)
        plot(time,cell2mat(AllFibers.z_fitted(ij,:)),'k-')

        subplot(3,4,1); xlabel('Time in s','Interpreter','latex'); ylabel('$x$ in mm','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        subplot(3,4,5); xlabel('Time in s','Interpreter','latex'); ylabel('$y$ in mm','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        subplot(3,4,9); xlabel('Time in s','Interpreter','latex'); ylabel('$z$ in mm','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')

        subplot(3,4,2)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.x_dot(ij,:)),'k-')
        subplot(3,4,6)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.y_dot(ij,:)),'k-')
        subplot(3,4,10)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.z_dot(ij,:)),'k-')

        subplot(3,4,2); xlabel('Time in s','Interpreter','latex'); ylabel('$u_x$ in mm/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        subplot(3,4,6); xlabel('Time in s','Interpreter','latex'); ylabel('$u_y$ in mm/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        subplot(3,4,10); xlabel('Time in s','Interpreter','latex'); ylabel('$u_z$ in mm/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')


        subplot(3,4,3)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.theta(ij,:)),'o','Color',color,'MarkerSize',5,'LineWidth',2)
        plot(time,cell2mat(AllFibers.theta_fitted(ij,:)),'k-')
        subplot(3,4,7)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.phi(ij,:)),'o','Color',color,'MarkerSize',5,'LineWidth',2)
        plot(time,cell2mat(AllFibers.phi_fitted(ij,:)),'k-')
        subplot(3,4,11)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.psi(ij,:)),'o','Color',color,'MarkerSize',5,'LineWidth',2)
        plot(time,cell2mat(AllFibers.psi_fitted(ij,:)),'k-')

        subplot(3,4,3); xlabel('Time in s','Interpreter','latex'); ylabel('$\theta$ in deg','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        ylim([-180 180])
        subplot(3,4,7); xlabel('Time in s','Interpreter','latex'); ylabel('$\phi$ in deg','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        ylim([-180 180])
        subplot(3,4,11); xlabel('Time in s','Interpreter','latex'); ylabel('$\psi$ in deg','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        ylim([-180 180])

        subplot(3,4,4)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.theta_dot(ij,:)),'k-')
        subplot(3,4,8)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.phi_dot(ij,:)),'k-')
        subplot(3,4,12)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.psi_dot(ij,:)),'k-')

        subplot(3,4,4); xlabel('Time in s','Interpreter','latex'); ylabel('$\dot{\theta}$ in deg/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        %ylim([-180 180])
        subplot(3,4,8); xlabel('Time in s','Interpreter','latex'); ylabel('$\dot{\phi}$ in deg/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        %ylim([-180 180])
        subplot(3,4,12); xlabel('Time in s','Interpreter','latex'); ylabel('$\dot{\psi}$ in deg/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        %ylim([-180 180])


        figure(); box on; grid on;
        fig=gcf; fig.Position=[801 1 600 600];

        subplot(3,4,1)
        title('Euler angles')
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.theta(ij,:)),'o','Color',color,'MarkerSize',5,'LineWidth',2)
        plot(time,cell2mat(AllFibers.theta_fitted(ij,:)),'k-')
        subplot(3,4,5)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.phi(ij,:)),'o','Color',color,'MarkerSize',5,'LineWidth',2)
        plot(time,cell2mat(AllFibers.phi_fitted(ij,:)),'k-')
        subplot(3,4,9)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.psi(ij,:)),'o','Color',color,'MarkerSize',5,'LineWidth',2)
        plot(time,cell2mat(AllFibers.psi_fitted(ij,:)),'k-')

        subplot(3,4,1); xlabel('Time in s','Interpreter','latex'); ylabel('$\theta$ in deg','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        ylim([-180 180])
        subplot(3,4,5); xlabel('Time in s','Interpreter','latex'); ylabel('$\phi$ in deg','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        ylim([-180 180])
        subplot(3,4,9); xlabel('Time in s','Interpreter','latex'); ylabel('$\psi$ in deg','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        ylim([-180 180])

        subplot(3,4,2)
        title('Euler angles rates')
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.theta_dot(ij,:)),'k-')
        subplot(3,4,6)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.phi_dot(ij,:)),'k-')
        subplot(3,4,10)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.psi_dot(ij,:)),'k-')

        subplot(3,4,2); xlabel('Time in s','Interpreter','latex'); ylabel('$\dot{\theta}$ in deg/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        %ylim([-180 180])
        subplot(3,4,6); xlabel('Time in s','Interpreter','latex'); ylabel('$\dot{\phi}$ in deg/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        %ylim([-180 180])
        subplot(3,4,10); xlabel('Time in s','Interpreter','latex'); ylabel('$\dot{\psi}$ in deg/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        %ylim([-180 180])


        subplot(3,4,3)
        title('Mobin - fiber ref.')
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.Mobin_omega_x(ij,:)),'k-')
        subplot(3,4,7)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.Mobin_omega_y(ij,:)),'k-')
        subplot(3,4,11)
        hold on; box on; grid on;
        plot(time,cell2mat(AllFibers.Mobin_omega_z(ij,:)),'k-')

        subplot(3,4,3); xlabel('Time in s','Interpreter','latex'); ylabel('$\omega_x$ in deg/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        %ylim([-180 180])
        subplot(3,4,7); xlabel('Time in s','Interpreter','latex'); ylabel('$\omega_y$ in deg/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        %ylim([-180 180])
        subplot(3,4,11); xlabel('Time in s','Interpreter','latex'); ylabel('$\omega_z$ in deg/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        %ylim([-180 180])



        subplot(3,4,4)
        plot(time,cell2mat(AllFibers.Quaternion_omega_x(ij,:)),'k-')
        title('Quaternion - lab ref.')
        subplot(3,4,8)
        plot(time,cell2mat(AllFibers.Quaternion_omega_y(ij,:)),'k-')
        subplot(3,4,12)
        plot(time,cell2mat(AllFibers.Quaternion_omega_z(ij,:)),'k-')

        subplot(3,4,4); xlabel('Time in s','Interpreter','latex'); ylabel('$\omega_qx$ in deg/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        %ylim([-180 180])
        subplot(3,4,8); xlabel('Time in s','Interpreter','latex'); ylabel('$\omega_qy$ in deg/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        %ylim([-180 180])
        subplot(3,4,12); xlabel('Time in s','Interpreter','latex'); ylabel('$\omega_qz$ in deg/s','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
        %ylim([-180 180])

        %timesteps=timesteps(1:floor(numel(timesteps)/2));
        if plot_fiber==1
        %%% check eigenvectors
            %II=zeros(2100,1600,700);
            II=[];
            max_x=max(cell2mat(AllFibers.x(ij,:))/dx);
            max_y=max(cell2mat(AllFibers.y(ij,:))/dx);
            max_z=max(cell2mat(AllFibers.z(ij,:))/dx);

            %BB=AllFibers.BoundingBox{ij,timesteps(end)};
            II=zeros(ceil(max_x)+100,ceil(max_y)+100,ceil(max_z)+100);
            pos=1;
            for it=timesteps(1:1:end)'
%                 loc=[AllFibers.Centroid_Refined{ij,it}(1),AllFibers.Centroid_Refined{ij,it}(2),AllFibers.Centroid_Refined{ij,it}(3)];
                dim=size(full(AllFibers.Object{ij,it}));
                bounds=AllFibers.BoundingBox{ij,it};
                %II(floor(loc(2))-floor(dim(2)/2):floor(loc(2))+floor(dim(2)/2),50+floor(loc(1))-floor(dim(1)/2):50+floor(loc(1))+floor(dim(1)/2)-1,50+floor(loc(3))-floor(dim(3)/2):50+floor(loc(3))+floor(dim(3)/2)-1)=permute(full(AllFibers.Object{ij,it}),[2 1 3]);
                II(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2))=II(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2))+full(AllFibers.Object{ij,it});
                pos=pos+1;
            end


            
            figure(); hold on;
            fig=gcf; fig.Position=[1 601 1800 400];
            colormap(cool(numel(levellist)))
            for i=1:length(levellist)
                level=levellist(i);
                p=patch(isosurface(permute(II,[2 1 3]),level));
                p.FaceVertexCData=level;
                p.FaceColor='flat';
                p.EdgeColor='none';
                p.FaceAlpha=(alphalist(i));
            end
            %ax=get(gca,'Children');
            daspect([1 1 1])
            hold on
%             for it=timesteps'
% 
%                 r=[AllFibers.EigenVectors{ij,it}{1,1}(2,1),AllFibers.EigenVectors{ij,it}{1,1}(1,1),AllFibers.EigenVectors{ij,it}{1,1}(3,1)];
%                 g=[AllFibers.EigenVectors{ij,it}{1,1}(2,2),AllFibers.EigenVectors{ij,it}{1,1}(1,2),AllFibers.EigenVectors{ij,it}{1,1}(3,2)];
%                 b=[AllFibers.EigenVectors{ij,it}{1,1}(2,3),AllFibers.EigenVectors{ij,it}{1,1}(1,3),AllFibers.EigenVectors{ij,it}{1,1}(3,3)];
% 
%                 quiver3(AllFibers.Centroid_Refined{ij,it}(1),AllFibers.Centroid_Refined{ij,it}(2),AllFibers.Centroid_Refined{ij,it}(3),...
%                     r(1),r(2),r(3),30,'LineWidth',2,'Color','m')
%                 
%                 quiver3(AllFibers.Centroid_Refined{ij,it}(1),AllFibers.Centroid_Refined{ij,it}(2),AllFibers.Centroid_Refined{ij,it}(3),...
%                     g(1),g(2),g(3),30,'LineWidth',2,'Color','y')
%                 
%                 quiver3(AllFibers.Centroid_Refined{ij,it}(1),AllFibers.Centroid_Refined{ij,it}(2),AllFibers.Centroid_Refined{ij,it}(3),...
%                     b(1),b(2),b(3),30,'LineWidth',2,'Color','c')
%             end
%             daspect([1 1 1])


            for it=timesteps(1:1:end)'

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
            daspect([1 1 1])

            box on
            grid on 
            grid minor
            view(3)
        end

        disp(' Press any key to see next fiber ')
        pause
        close all
    end
end