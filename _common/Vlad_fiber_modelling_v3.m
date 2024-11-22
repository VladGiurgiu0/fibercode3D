function [Centroid_refined,Curvature,Length,Length_effective,EigenVectors_Tensor,Orientation_Tensor,px,py,pz,PrincipalAxisLength]=Vlad_fiber_modelling(Fiber,Bounding_Box,v,p)


%% processing
%Fiber=imgaussfilt3(Fiber,1);
switch p.peak_finding_technique
    case "Dilation"
        I_peaks_1=Vlad_peak_finder_dilation_3D(Fiber,'vertical',p.kernel_dilation);
        I_peaks_2=Vlad_peak_finder_dilation_3D(Fiber,'horizontal',p.kernel_dilation);
        I_peaks_3=Vlad_peak_finder_dilation_3D(Fiber,'spanwise',p.kernel_dilation);
        I_peaks=I_peaks_1.*I_peaks_2.*I_peaks_3;
    case "imregionalmax"
        I_peaks=imregionalmax(Fiber,26);
    case "max"
        % find the maxima in each direction
        [~,idx_2] = max(Fiber,[],1);
        [~,idx_1] = max(Fiber,[],2);
        [~,idx_3] = max(Fiber,[],3);
        % delete the edge points
        idx_1(idx_1==1)=NaN;
        idx_2(idx_2==1)=NaN;
        idx_3(idx_3==1)=NaN;

        [X,Y,Z]=meshgrid(1:size(Fiber,2),1:size(Fiber,1),1:size(Fiber,3));
        max_2=[reshape(X(1,:,:),1,[]);reshape(idx_2,1,[]);reshape(Z(1,:,:),1,[])];
        max_1=[reshape(idx_1,1,[]);reshape(Y(:,1,:),1,[]);reshape(Z(:,1,:),1,[])];
        max_3=[reshape(X(:,:,1),1,[]);reshape(Y(:,:,1),1,[]);reshape(idx_3,1,[])];
        %hold all

        I_peaks=zeros(size(Fiber));
        max_1=max_1(:,~isnan(max_1(1,:)));
        max_2=max_2(:,~isnan(max_2(2,:)));
        max_3=max_3(:,~isnan(max_3(3,:)));

        for jj=1:size(max_1(1,:),2)
            I_peaks(max_1(2,jj),max_1(1,jj),max_1(3,jj))=1;
        end

        for jj=1:size(max_2(1,:),2)
            I_peaks(max_2(2,jj),max_2(1,jj),max_2(3,jj))=1;
        end
        
        for jj=1:size(max_3(1,:),2)
            I_peaks(max_3(2,jj),max_3(1,jj),max_3(3,jj))=1;
        end

        I_peaks=bwmorph3(I_peaks,'majority');

    case "Skeletonize"
        bottom_percentile=prctile(reshape(Fiber,1,[]), p.percentile_majority_skeletonize);
        I_peaks= bwskel(imbinarize(Fiber,bottom_percentile),'MinBranchLength',p.min_branch_length_skeletonize);

    case "Majority"
        bottom_percentile=prctile(reshape(Fiber,1,[]), p.percentile_majority_skeletonize);
        I_peaks=bwmorph3(imbinarize(Fiber,bottom_percentile),'majority');
        if sum(sum(sum(I_peaks,3),2),1)==0 % if the fibre is very thin, then 'majority' does not set any voxel to 1 so this condition has to be checked
            I_peaks=imbinarize(Fiber,bottom_percentile);
        end
        %I_peaks=bwskel(I_peaks);

    case "Majority + Skeletonize"
        bottom_percentile=prctile(reshape(Fiber,1,[]), p.percentile_majority_skeletonize);
        I_peaks=bwmorph3(imbinarize(Fiber,bottom_percentile),'majority');
        I_peaks=bwskel(I_peaks);
end

props = regionprops3(imbinarize(Fiber,p.imbin_thres),'PrincipalAxisLength');
PrincipalAxisLength = props.PrincipalAxisLength;

[b,a,c]=ind2sub(size(I_peaks),find(I_peaks>0));
d=Fiber(I_peaks>0);


switch p.polynomial_finding_technique
    case "Mobin"
        [px,poly_x,~,py,poly_y,~,pz,poly_z,~]=Vlad_fit_fiber_curviliniar(a,b,c,d,100,p.use_weights,p.fitting_type_fibre,Bounding_Box);
        %[px_2,poly_x_2,~,py_2,poly_y_2,~,pz_2,poly_z_2,~]=Vlad_fit_fiber_curviliniar(a,b,c,d,100,1);
end


%% find reference frame of the fitted polynomial

[Curvature,Length]=Vlad_compute_curvature_length(px,py,pz);
Length_effective=norm([px(1) py(1) pz(1)]-[px(end) py(end) pz(end)]);

[red,green,blue,Centroid_refined]=Vlad_find_orientation_vectors_angles(px',py',pz','tensor');
x=Centroid_refined(1);
y=Centroid_refined(2);
z=Centroid_refined(3);

Centroid_refined(1)=Centroid_refined(1)+Bounding_Box(1,1);
Centroid_refined(2)=Centroid_refined(2)+Bounding_Box(2,1);
Centroid_refined(3)=Centroid_refined(3)+Bounding_Box(3,1);
EigenVectors_Tensor(1,:)=red;
EigenVectors_Tensor(2,:)=green;
EigenVectors_Tensor(3,:)=blue;



Orientation_Tensor=Vlad_find_euler_angles(red,green,blue);

%%%%% add warning if multiple Fibers in the original object have been identified or in the region there are multiple objects %%%%%%

    px=px+Bounding_Box(1,1);
    py=py+Bounding_Box(2,1);
    pz=pz+Bounding_Box(3,1);



%% plot to check
if p.plot==1
    x=x+Bounding_Box(1,1);
    y=y+Bounding_Box(2,1);
    z=z+Bounding_Box(3,1);
%     px=px+Bounding_Box(1,1);
%     py=py+Bounding_Box(2,1);
%     pz=pz+Bounding_Box(3,1);

    a=a+Bounding_Box(1,1);
    b=b+Bounding_Box(2,1);
    c=c+Bounding_Box(3,1);

    [A, B, C]=meshgrid(1:size(Fiber,1),1:size(Fiber,2),1:size(Fiber,3));
    A=A+Bounding_Box(2,1);
    B=B+Bounding_Box(1,1);
    C=C+Bounding_Box(3,1);


figure(1); 
fig=gcf; clf ;fig.WindowState='maximized';%fig.Position=[20 100 500 500];
daspect([1 1 1])
hold all

%%% plot fiber intensity with isosurface
for i=1:length(p.levellist)
    level=p.levellist(i);
    pat=patch(isosurface(permute(B,[2 1 3]),permute(A,[2 1 3]),permute(C,[2 1 3]),Fiber,level));
    %pat=patch(isosurface(Fiber,level));
    %pat=patch(isosurface(A,B,C,permute(Fiber,[2 1 3]),level));
    pat.FaceVertexCData=level;
    pat.FaceColor='flat';
    pat.EdgeColor='none';
    pat.FaceAlpha=(p.facealphalist(i));

end

%%% plot fiber intensity with plotcube
% Fiber2=permute(Fiber,[1 2 3]);
% intensities=sort(nonzeros(Fiber2));
% int=linspace(0.01,max(intensities),numel(intensities))';
% s=regionprops3(imbinarize(Fiber2,0.01),'all');
% Vox=s.VoxelList{1,1};
% for ij=1:20:numel(Vox(:,1))
%     plotcube([1 1 1],[Bounding_Box(1,1)+Vox(ij,2)-0.5 Bounding_Box(2,1)+Vox(ij,1)-0.5 Bounding_Box(3,1)+Vox(ij,3)-0.5],(int(ij)/max(intensities).^1.5),[255 45 85]/255);
% end
% level=0.2;

%%% plot shadows
% p=patch(isosurface(x+80*ones(size(A)),B,C,permute(Fiber,[2 1 3]),0.2)); p.FaceVertexCData=level; p.FaceColor='flat'; p.EdgeColor='none'; p.FaceColor=[0.9 0.9 0.9]; p.FaceAlpha=(1);
% p=patch(isosurface(A,y+30*ones(size(B)),C,permute(Fiber,[2 1 3]),0.2)); p.FaceVertexCData=level; p.FaceColor='flat'; p.EdgeColor='none'; p.FaceColor=[0.9 0.9 0.9]; p.FaceAlpha=(1);
% p=patch(isosurface(A,B,z-51*ones(size(C)),permute(Fiber,[2 1 3]),0.2)); p.FaceVertexCData=level; p.FaceColor='flat'; p.EdgeColor='none'; p.FaceColor=[0.9 0.9 0.9]; p.FaceAlpha=(1);
% 
% p=patch(isosurface(x+79.5*ones(size(A)),B,C,permute(Fiber,[2 1 3]),0.9)); p.FaceVertexCData=level; p.FaceColor='flat'; p.EdgeColor='none'; p.FaceColor=[0.7 0.7 0.7]; p.FaceAlpha=(1);
% p=patch(isosurface(A,y+29.5*ones(size(B)),C,permute(Fiber,[2 1 3]),0.9)); p.FaceVertexCData=level; p.FaceColor='flat'; p.EdgeColor='none'; p.FaceColor=[0.7 0.7 0.7]; p.FaceAlpha=(1);
% p=patch(isosurface(A,B,z-50.5*ones(size(C)),permute(Fiber,[2 1 3]),0.9)); p.FaceVertexCData=level; p.FaceColor='flat'; p.EdgeColor='none'; p.FaceColor=[0.7 0.7 0.7]; p.FaceAlpha=(1);
% 
% plot3(px,py,z-50*ones(size(pz)),'-','LineWidth',3,'Color',[255 68 59]/255)
% plot3(px,y+29*ones(size(py)),pz,'-','LineWidth',3,'Color',[50 215 75]/255)
% plot3(x+79*ones(size(px)),py,pz,'-','LineWidth',3,'Color',[10 132 255]/255)

daspect([1 1 1]); box on; %grid on; %grid minor;

%%% plot found intensity points
% scatter3(a,b,c,10,'MarkerFaceColor','r','MarkerEdgeColor','none',...
%     'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3)

%%% plot polynomial of the fiber
%plot3(px,py,pz,'k-','LineWidth',5)


%scatter3(px_2,py_2,pz_2,100,'g.')
%scatter3(b,a,surface,100,'g.')
%[Xs, Ys, Zs]=meshgrid(1:size(Fiber,1),1:size(Fiber,2),1:size(Fiber,3));
%salami = isosurface(Xs,Ys,Zs,Fiber_salami,0);
%patch(salami,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0 0 1],'FaceAlpha',0.5);

%%% plot reference frame
% quiver3(x,y,z,red(1),red(2),red(3),20,'LineWidth',2,'Color','r','Marker','.','MarkerSize',30,'MarkerFaceColor','k','MarkerEdgeColor','k','MaxHeadSize',150)
% quiver3(x,y,z,green(1),green(2),green(3),20,'LineWidth',2,'Color','g','MaxHeadSize',150)
% quiver3(x,y,z,blue(1),blue(2),blue(3),20,'LineWidth',2,'Color','b','MaxHeadSize',150)
% 
% quiver3(x-50,y-50,z-50,1,0,0,15,'filled','LineWidth',2,'Color',[255 68 59]/255,'Marker','.','MarkerSize',20,'MarkerFaceColor','k','MarkerEdgeColor',[152 152 157]/255,'MaxHeadSize',150);
% quiver3(x-50,y-50,z-50,0,1,0,15,'filled','LineWidth',2,'Color',[50 215 75]/255,'MaxHeadSize',150)
% quiver3(x-50,y-50,z-50,0,0,1,15,'filled','LineWidth',2,'Color',[10 132 255]/255,'MaxHeadSize',150)
view(3)
%xlim([x-51 x+80])
%ylim([y-51 y+30])
%zlim([z-51 z+50])

%set(gca,'TickLabelInterpreter','latex','FontSize',15,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[])
% xlabel('$x$','Interpreter','latex')
% ylabel('$y$','Interpreter','latex')
% zlabel('$z$','Interpreter','latex')


%campos(1e3*[1.0865, 0.1309, 0.3180])
%camtarget(1e3*[1.3913    0.4502    0.1221])
%camzoom(3)

xlim([x-p.limits_plot x+p.limits_plot])
ylim([y-p.limits_plot y+p.limits_plot])
zlim([z-p.limits_plot z+p.limits_plot])

xlimit=xlim();
ylimit=ylim();
zlimit=zlim();

scatter3(a,b,(zlimit(1))*ones(size(c)),10,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none',...
    'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)

scatter3(a,(ylimit(2))*ones(size(b)),c,10,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none',...
    'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)

scatter3((xlimit(2))*ones(size(a)),b,c,10,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none',...
    'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)

% scatter3(a,b,(zlimit(1)-5)*ones(size(c)),10,-log(1-(d(:,1)/max(d(:,1)))).*[1,0,0],'filled','MarkerEdgeColor','none',...
%     'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
% 
% scatter3(a,(ylimit(2)+10)*ones(size(b)),c,10,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none',...
%     'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
% 
% scatter3((xlimit(2)+15)*ones(size(a)),b,c,10,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none',...
%     'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)

xlimit=xlim();
ylimit=ylim();
zlimit=zlim();
quiver3(xlimit(1)+2,ylimit(1)+2,zlimit(1),1,0,0,10,'filled','LineWidth',2,'Color',[255 68 59]/255,'Marker','.','MarkerSize',20,'MarkerFaceColor','k','MarkerEdgeColor',[152 152 157]/255,'MaxHeadSize',20);
quiver3(xlimit(1)+2,ylimit(1)+2,zlimit(1),0,1,0,10,'filled','LineWidth',2,'Color',[50 215 75]/255,'MaxHeadSize',20)
quiver3(xlimit(1)+2,ylimit(1)+2,zlimit(1),0,0,1,10,'filled','LineWidth',2,'Color',[10 132 255]/255,'MaxHeadSize',20)

%xlabel('$x$ (stream-wise)')
%ylabel('$y$ (span-wise)')
%zlabel('$z$ (wall-normal)')
set(gca,'XTick',[],'YTick',[],'ZTick',[])

plot3(px,py,zlimit(1)*ones(size(pz)),'-','LineWidth',3,'Color',[255 68 59]/255)
plot3(px,ylimit(2)*ones(size(py)),pz,'-','LineWidth',3,'Color',[50 215 75]/255)
plot3(xlimit(2)*ones(size(px)),py,pz,'-','LineWidth',3,'Color',[10 132 255]/255)

plot3(px,py,pz,'k-','LineWidth',5)

% show the plane used for fitting the fibre
% if p.polynomial_finding_technique=="Plane fitting"
%     surf(x_temp+x,y_temp+y,z_temp+z,'EdgeColor','none','FaceColor','k','FaceAlpha',0.2)
% end

quiver3(x,y,z,red(1),red(2),red(3),10,'LineWidth',2,'Color','r','Marker','.','MarkerSize',30,'MarkerFaceColor','k','MarkerEdgeColor','k','MaxHeadSize',150)
quiver3(x,y,z,green(1),green(2),green(3),10,'LineWidth',2,'Color','g','MaxHeadSize',150)
quiver3(x,y,z,blue(1),blue(2),blue(3),10,'LineWidth',2,'Color','b','MaxHeadSize',150)




if p.make_movie==1

frame=getframe(gcf);
writeVideo(v,frame);
end

if p.pause_enabled==1
    disp('Paused')
    pause
end

end
end