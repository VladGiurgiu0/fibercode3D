clc; clear; close all;
close all
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

addpath("_common\")

color_1 = colors_vlad('2','CelestialBlueTangerineDreamEnchantedForest');
color_2 = colors_vlad('3','CelestialBlueTangerineDreamEnchantedForest');
color_3 = colors_vlad('1','CelestialBlueTangerineDreamEnchantedForest');

%% load data
%load('G:\__PRL\processed_data_close_wall\region_sgolay_ts_20_fk_20\Loop=2\3_Refined_fibers\AllFibers.mat')
%load('G:\__PRL\processed_data_close_wall\region_sgolay_ts_20_fk_20\Loop=2\4_Quantities_Refined_fibers\AllFibers_Only_data.mat')

load('K:\Matteo_fibers\Results - Copy\all_sets\Loop=4\4_Quantities_Refined_fibers\AllFibers_Only_data.mat')

%% show trajectories of fibers
p.pause_enabled=0;
p.print=0;
plot_tracks(p,AllFibers)
fig1=gcf; %fig1.Renderer='painters';

%% plot fiber by fiber with rotation rates
fig_tracks=figure(); hold all; box on; grid on;
fig_tracks.WindowState='maximized';

skip = 5;

nr_fibers = size(AllFibers.blue_tensor,1);
nr_timesteps = size(AllFibers.blue_tensor,2);

for ij=1:1:nr_fibers
    %clf;
    color=rand(1,3); 
    x=[]; y=[]; z=[];
    omega_s_x =[]; omega_s_y =[]; omega_s_z =[];

    for ii=1:5:nr_timesteps
        if ~isempty(AllFibers.Centroid{ij,ii})
            x(ii)=AllFibers.Centroid{ij,ii}(1);
            y(ii)=AllFibers.Centroid{ij,ii}(2);
            z(ii)=AllFibers.Centroid{ij,ii}(3);

            omega_s_x(ii) = AllFibers.Omega_s_x_RotationMatrix{ij,ii};
            omega_s_y(ii) = AllFibers.Omega_s_y_RotationMatrix{ij,ii};
            omega_s_z(ii) = AllFibers.Omega_s_z_RotationMatrix{ij,ii};

            omega_s_x(omega_s_x==0)=NaN;
            omega_s_y(omega_s_y==0)=NaN;
            omega_s_z(omega_s_z==0)=NaN;

            px = AllFibers.px_Tensor{ij,ii};
            py = AllFibers.py_Tensor{ij,ii};
            pz = AllFibers.pz_Tensor{ij,ii};
            
            subplot(2,2,[1 2]); hold on; grid on; box on;
            % plot polynomials
            plot3(px,py,pz,'LineWidth',3,'Color',color,'LineStyle','-')

            view(3)
            daspect([1 1 1])
            
            xlabel('$x$ (stream-wise)')
            ylabel('$y$ (span-wise)')
            zlabel('$z$ (wall-normal)')

%             % plot shadows of polynomials
%             plot3(px,py,0*ones(size(pz)),'LineWidth',3,'Color',[0.7 0.7 0.7],'LineStyle','-')
%             plot3(px,0*ones(size(py)),pz,'LineWidth',3,'Color',[0.7 0.7 0.7],'LineStyle','-')
%             plot3(0*ones(size(px)),py,pz,'LineWidth',3,'Color',[0.7 0.7 0.7],'LineStyle','-')
        end
    end

%     quiver3(0,0,0,1,0,0,100,'filled','LineWidth',2,'Color',[255 68 59]/255,'Marker','.','MarkerSize',20,'MarkerFaceColor','k','MarkerEdgeColor',[152 152 157]/255,'MaxHeadSize',150);
%     quiver3(0,0,0,0,1,0,100,'filled','LineWidth',2,'Color',[50 215 75]/255,'MaxHeadSize',150)
%     quiver3(0,0,0,0,0,1,100,'filled','LineWidth',2,'Color',[10 132 255]/255,'MaxHeadSize',150)

    subplot(2,2,3); hold on; grid on; box on;
    p1 = plot(x,'.','Color',color_1,'MarkerSize',20);
    p2 = plot(y,'.','Color',color_2,'MarkerSize',20);
    p3 = plot(z,'.','Color',color_3,'MarkerSize',20);

    xlabel('$t$ [-]')
    ylabel('$x$, $y$, $z$ [vox]')
    set(gca,'YScale','log','FontSize',20)

    legend([p1 p2 p3],{'$x$','$y$','$z$'},'location','best','FontSize',20)

    subplot(2,2,4); hold on; grid on; box on;
    p4 = plot(abs(omega_s_x),'.','Color',color_1,'MarkerSize',20);
    p5 = plot(abs(omega_s_y),'.','Color',color_2,'MarkerSize',20);
    p6 = plot(abs(omega_s_z),'.','Color',color_3,'MarkerSize',20);

    xlabel('$t$ [-]')
    ylabel('$|\omega_x|$, $|\omega_y|$, $|\omega_z|$ [$^\circ /$t]')
    set(gca,'YScale','lin','FontSize',20)

    legend([p4 p5 p6],{'$|\omega_x|$','$|\omega_y|$','$|\omega_z|$'},'location','best','FontSize',20)

    pause
end


