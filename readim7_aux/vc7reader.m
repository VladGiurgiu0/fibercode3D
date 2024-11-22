clc
% close all
clear all

numsnap=500;
% fld='C:\Users\admin\Documents\LVExport\TomographicPIV_16x16x16_75%ov';
fld='C:\Users\admin\Documents\LVExport\TomographicPIV_48x48x48_75%ov';

for sn=1:numsnap
fname=strcat(fld,'\B0',num2str(sn,'%.3d'),'.vc7');

imx=readimx(fname);
for il=8:12 %il=1:numel(imx.Frames{1, 1}.Components{1, 1}.Planes)
I(:,:,il)=flip(imx.Frames{1, 1}.Components{1, 1}.Planes{il,1}');
end
I(I==0|I<0)=NaN; 
U_t(:,sn)=nanmean(nanmean(I,2),3);
disp(strcat('Image numebr:',num2str(sn)))
end
U=mean(U_t,2);
bin=numel(U);
delta=0.116;
dy=(22.89+28.8)/bin;
y=0:dy:(bin-1)*dy;
hold on
plot(y/delta,U)
set(gca, 'XScale', 'lin')