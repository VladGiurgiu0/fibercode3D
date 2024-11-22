clc ; clear ; close all
%addpath('/Users/vlad/Owncloud/Research/_codes/___common')
%addpath('C:\Users\Corsair\Desktop\Vlad\fibrecode\Vlad')
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% load data
output_spinning=NaN(100,8,11);
for ii=1:7
    for jj=[1 2 5 7 10 15 20]
        % load the inputs
        load(strcat("Generated_virtual_fibre_",num2str(ii),"\_parameters.mat"),'geo')
        length_prime(ii) = geo.length_prime; clear geo;
        angular_displacement(ii,jj) = jj;
        
        % raw data
        load(strcat('Generated_virtual_fibre_',num2str(ii),'\',num2str(jj),'_deg_per_timestep\4_Quantities_Refined_fibers\AllFibers.mat'))
        
        % fitered data movmean kernel 10
        %load(strcat('FILTERED_movmean_kernel10_Generated_virtual_fibre_',num2str(ii),'\',num2str(jj),'_deg_per_timestep\4_Quantities_Refined_fibers\AllFibers.mat'))
        
        % fitered data sgolay kernel 10
        %load(strcat('FILTERED_movmean_kernel10_Generated_virtual_fibre_',num2str(ii),'\',num2str(jj),'_deg_per_timestep\4_Quantities_Refined_fibers\AllFibers.mat'))
        
        AllFibers.Omega_b_x_RotationMatrix( cellfun(@isempty, AllFibers.Omega_b_x_RotationMatrix) ) = {NaN};
        AllFibers.Omega_b_x_RotationMatrix(numel(cellfun(@isempty, AllFibers.Omega_b_x_RotationMatrix)):100) = {NaN};
        output_spinning (:,jj,ii) = cell2mat(AllFibers.Omega_b_x_RotationMatrix);
    end
end

output_spinning(:,[3 4 6 8 9 11 12 13 14 16 17 18 19],:)=[];
output_spinning(output_spinning==0)=NaN;
angular_displacement(:,[3 4 6 8 9 11 12 13 14 16 17 18 19])=[];
angular_displacement = angular_displacement(1,:);

load("Generated_virtual_fibre_1\_parameters.mat","input_spinning_1","time")

%% compute error
for ii=1:7
    for jj=[1 2 3 4 5 6 7]
        error_spinning (:,jj,ii) = abs( (output_spinning (:,jj,ii) - input_spinning_1')./input_spinning_1' ) * 100; % in percent
    end
end
mean_error_spinning = squeeze(mean(error_spinning(1:end,:,:),1,'omitnan'));

%% plot
fig1=figure(1); box on; 
%levels=[0:4:20];
levels=[0,4,8,16,32,64];
colormap(parula(64))
[Cont, hCont]=contourf(angular_displacement,length_prime,mean_error_spinning',levels);
clabel(Cont,hCont,levels,'FontSize',20)
%set(gca,'ColorScale','log')
clim([0,64])
%contourLegend(hCont)
%colormap(parula(numel(levels)))
%colorbar
ylabel('$l^{\prime}$ [vox]')
xlabel('$\Delta \alpha \ [^\circ/$time-step]')
set(gca,'FontSize',16)
ylim([5 34])
set(gca,'XTick',[1 5 10 15 20],'YTick',[5 10 15 20 25 30 34],'Linewidth',2)

annotation(fig1,'textbox',...
    [0.530662104362703 0.828571428571429 0.367552181351583 0.0822170763567289],...
    'String','$|{\Delta\alpha}_m - \Delta\alpha|/\Delta\alpha \ [\%]$',...
    'Interpreter','latex',...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 1],'FontSize',16);

