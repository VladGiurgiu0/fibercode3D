function plot_tracks(p, AllFibers)
%%% plots tracks as points of the same color

    fig_tracks=figure(); hold on; box on; grid on;
    fig_tracks.WindowState='maximized';

    for tt=1:max(1,floor(p.ki/p.max_track_length))
        starting=(tt-1)*p.max_track_length;
        ending=tt*p.max_track_length;

        for ij=1:size(AllFibers.Centroid(:,1),1)
            color=rand(1,3); x=[]; y=[]; z=[];
            for ii=starting+1:ending
                if ~isempty(AllFibers.Centroid{ij,ii})
                    x(ii)=AllFibers.Centroid{ij,ii}(1);
                    y(ii)=AllFibers.Centroid{ij,ii}(2);
                    z(ii)=AllFibers.Centroid{ij,ii}(3);
                end
            end
            scatter3(x,y,z,'filled','MarkerFaceColor',color)
            %drawnow
            %pause
        end
        %drawnow
        view(3)
        daspect([1 1 1])

        xlabel('$x$ (stream-wise)')
        ylabel('$y$ (span-wise)')
        zlabel('$z$ (wall-normal)')
    
        quiver3(0,0,0,1,0,0,100,'filled','LineWidth',2,'Color',[255 68 59]/255,'Marker','.','MarkerSize',20,'MarkerFaceColor','k','MarkerEdgeColor',[152 152 157]/255,'MaxHeadSize',150);
        quiver3(0,0,0,0,1,0,100,'filled','LineWidth',2,'Color',[50 215 75]/255,'MaxHeadSize',150)
        quiver3(0,0,0,0,0,1,100,'filled','LineWidth',2,'Color',[10 132 255]/255,'MaxHeadSize',150)

        if p.print==1
            if ~exist(strcat(p.save,'Figures_Movies_Processing\3_Tracking\'),'dir')
                mkdir(strcat(p.save,'Figures_Movies_Processing\3_Tracking\'))
            end
            savefig(strcat(p.save,'Figures_Movies_Processing\3_Tracking\','From_timestep_',num2str(starting+1),'_to_',num2str(ending),'.fig'))
            print(strcat(p.save,'Figures_Movies_Processing\3_Tracking\','From_timestep_',num2str(starting+1),'_to_',num2str(ending),'.tif'),'-dtiffn')
        end

        if p.pause_enabled==1
            pause
        end

        %close all
    end 
end