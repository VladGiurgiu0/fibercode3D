clc; clear; close all;
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

addpath("_common\")
%% load data
load('G:\__PRL\processed_data_close_wall\region_sgolay_ts_20_fk_20\Loop=2\3_Refined_fibers\AllFibers.mat')

%% compute
for ij=1:size(AllFibers.Centroid(:,1),1)
    index_time=find(~cellfun('isempty',AllFibers.Centroid(ij,:)))';

    for it=index_time(1:end-1)'

        % region
        r_red = [AllFibers.EigenVectors{ij,it}{1,1}(1,1),AllFibers.EigenVectors{ij,it}{1,1}(2,1),AllFibers.EigenVectors{ij,it}{1,1}(3,1)];
        r_green=[AllFibers.EigenVectors{ij,it}{1,1}(1,2),AllFibers.EigenVectors{ij,it}{1,1}(2,2),AllFibers.EigenVectors{ij,it}{1,1}(3,2)];
        r_blue=[AllFibers.EigenVectors{ij,it}{1,1}(1,3),AllFibers.EigenVectors{ij,it}{1,1}(2,3),AllFibers.EigenVectors{ij,it}{1,1}(3,3)];
    
        % svd           
        PrincAxis = Beppe_PricipalAxis(full(AllFibers.Object{ij,it}));
        s_red=PrincAxis(:,1);
        s_green=PrincAxis(:,2);
        s_blue=PrincAxis(:,3);
 
        disp(strcat('is red equal? ',num2str(isequal(r_red,s_red))))
        disp(strcat('is green equal? ',num2str(isequal(r_green,s_green))))
        disp(strcat('is blue equal? ',num2str(isequal(r_blue,s_blue)))) 

        disp(strcat('region red - svd red = ',num2str(abs(r_red) - abs(s_red)')))
        disp(strcat('region green - svd green = ',num2str(abs(r_green) - abs(s_green'))))
        disp(strcat('region blue - svd blue = ',num2str(abs(r_blue) - abs(s_blue'))))

        pause

    end

end

%% -------------- conclusion is that svd and region are equivalent within machine precision -------------- %%



