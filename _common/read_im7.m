function[I]=read_im7(f1,kk,ik)

%  This function reads and store the information from 
%  im7 formatted 2D planes of raw data to a 3D matrix. 

% Version 1 - 15/02/2020 - Mobin Alipour

fname=strcat(f1,'S',num2str(kk,'%.5d'),'\B0',num2str(ik,'%.4d'),'.im7');
imx=readimx(fname);
II=imx.Frames{1}.Components{1}.Planes{1};
I=zeros(size(II,1),size(II,2),ik);

%parfor ii=1:ik
for ii=1:ik
        
%     fprintf('Reading %d of %d \n',ii,ik)
    
    fname=strcat(f1,'S',num2str(kk,'%.5d'),'\B0',num2str(ii,'%.4d'),'.im7');
    imx=readimx(fname);
    %II=imx.Frames{1}.Components{1}.Planes{1};
    I(:,:,ii)=imx.Frames{1}.Components{1}.Planes{1};
end
