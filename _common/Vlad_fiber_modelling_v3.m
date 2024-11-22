function [Centroid_refined,Curvature,Length,Length_effective,EigenVectors_Tensor,Orientation_Tensor,px,py,pz,PrincipalAxisLength]=Vlad_fiber_modelling(Fiber,Bounding_Box,v,p)


%% processing
%Fiber=imgaussfilt3(Fiber,1);
switch p.peak_finding_technique
    case "Dilation"
        kernel=7;
        I_peaks_1=Vlad_peak_finder_dilation_3D(Fiber,'vertical',kernel_dilation);
        I_peaks_2=Vlad_peak_finder_dilation_3D(Fiber,'horizontal',kernel_dilation);
        I_peaks_3=Vlad_peak_finder_dilation_3D(Fiber,'spanwise',kernel_dilation);
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

%[a,b,c]=ind2sub(size(I_peaks),find(I_peaks>0)); % x and y are switched (DON'T uncomment)
[b,a,c]=ind2sub(size(I_peaks),find(I_peaks>0));
d=Fiber(I_peaks>0);


switch p.polynomial_finding_technique
    case "Mobin"
        [px,poly_x,~,py,poly_y,~,pz,poly_z,~]=Vlad_fit_fiber_curviliniar(a,b,c,d,100,p.use_weights,p.fitting_type_fibre,Bounding_Box);
        %[px_2,poly_x_2,~,py_2,poly_y_2,~,pz_2,poly_z_2,~]=Vlad_fit_fiber_curviliniar(a,b,c,d,100,1);
    case "Plane fitting"
%% 

%%%% first approximation
        % make figure
        if p.plot==1
        figure(6); clf; hold on; grid on; box on; daspect([1 1 1]); view(3)
        set(gca,'TickLabelInterpreter','latex')
        xlabel('$x$','Interpreter','latex')
        ylabel('$y$','Interpreter','latex')
        zlabel('$z$','Interpreter','latex')
        % lab coordinate system
        quiver3(0,0,0,1,0,0,5,'filled','LineWidth',2,'Color','r','MaxHeadSize',10)
        quiver3(0,0,0,0,1,0,5,'filled','LineWidth',2,'Color','g','MaxHeadSize',10)
        quiver3(0,0,0,0,0,1,5,'filled','LineWidth',2,'Color','b','MaxHeadSize',10)

        % plot fiber intensity with isosurface
        [A, B, C]=meshgrid(1:size(Fiber,1),1:size(Fiber,2),1:size(Fiber,3));
        for i=1:length(p.levellist)
            level=p.levellist(i);
            pat=patch(isosurface(permute(B,[2 1 3]),permute(A,[2 1 3]),permute(C,[2 1 3]),Fiber,level));
            pat.FaceVertexCData=level;
            pat.FaceColor='flat';
            pat.EdgeColor='none';
            pat.FaceAlpha=(p.facealphalist(i));
        end
        colormap(hsv(numel(p.levellist)))
        clb=colorbar;
        clb.Limits=[min(p.levellist), max(p.levellist)];
        clb.Location="southoutside";

        % plot fiber peaks
        scatter3(a,b,c,10,'MarkerFaceColor',[255 69 58]/255,'MarkerEdgeColor','none',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
        end

        % find approximate orientation of the fiber (through regionprops3)
        ss=regionprops3(imbinarize(Fiber,0),'EigenVectors','BoundingBox','Centroid');
        eigenvectors_peaks=ss.EigenVectors{1,1};
        Bounding_Box_peaks=ss.BoundingBox;
        centroid=ss.Centroid;
        loc_x=centroid(2);
        loc_y=centroid(1);
        loc_z=centroid(3);

        red_peaks=eigenvectors_peaks(:,1); 
        %red_peaks=red_peaks([2 1 3]);           % switch x and y
        green_peaks=eigenvectors_peaks(:,3);
        %green_peaks=green_peaks([2 1 3]);
        blue_peaks=eigenvectors_peaks(:,2);
        %blue_peaks=-blue_peaks([2 1 3]);
        if p.plot==1
        quiver3(loc_x,loc_y,loc_z,red_peaks(1),red_peaks(2),red_peaks(3),10,'r','LineWidth',3,'MaxHeadSize',100)
        quiver3(loc_x,loc_y,loc_z,green_peaks(1),green_peaks(2),green_peaks(3),10,'g','LineWidth',3,'MaxHeadSize',100)
        quiver3(loc_x,loc_y,loc_z,blue_peaks(1),blue_peaks(2),blue_peaks(3),10,'b','LineWidth',3,'MaxHeadSize',100)
        end
        % align the fiber peaks (locations x,y,z are stored in a,b,c) with the lab coordinate system
        A_B_C=[a';b';c'];
        %Rm_peaks=quat2rotm(quaternion([green_peaks';blue_peaks';red_peaks'],'rotmat','frame'));
        %Rm_peaks=quat2rotm(quaternion([red_peaks';blue_peaks';green_peaks'],'rotmat','frame'));
        %Rm_peaks=quat2rotm(quaternion([blue_peaks';green_peaks';red_peaks'],'rotmat','point'));

%         Rm_peaks2 = [dot(red_peaks, [1;0;0]), dot(green_peaks, [1;0;0]), dot(blue_peaks, [1;0;0]);     
%                     dot(red_peaks, [0;1;0]), dot(green_peaks, [0;1;0]), dot(blue_peaks, [0;1;0]);     
%                     dot(red_peaks, [0;0;1]), dot(green_peaks, [0;0;1]), dot(blue_peaks, [0;0;1])];

        % Create quaternions for the two coordinate systems
        q1 = quaternion(rotm2quat([red_peaks';green_peaks';blue_peaks']));
        q2 = quaternion(rotm2quat([1 0 0; 0 1 0; 0 0 1]));
        % Compute the relative quaternion rotation
        qrel = q2 * conj(q1);
        % Convert the relative quaternion rotation to a rotation matrix
        Rm_peaks = rotmat(qrel,'frame');

            Rm_peaks=[red_peaks blue_peaks green_peaks]';
    
            A_B_C_prime=Rm_peaks*A_B_C;
            a_rot=A_B_C_prime(1,:)';
            b_rot=A_B_C_prime(2,:)';
            c_rot=A_B_C_prime(3,:)';
            if p.plot==1
            % plot rotated fiber peaks
            scatter3(a_rot,b_rot,c_rot,10,'MarkerFaceColor',[255 159 10]/255,'MarkerEdgeColor','none',...
                'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
            end
            % sort the locations of the fiber peaks (a,b,c) based on bounding
            % box sizes
    %         d1=Bounding_Box(1,2)-Bounding_Box(1,1);
    %         d2=Bounding_Box(2,2)-Bounding_Box(2,1);
    %         d3=Bounding_Box(3,2)-Bounding_Box(3,1);
    
            d1=Bounding_Box_peaks(1,5)-Bounding_Box_peaks(1,2);
            d2=Bounding_Box_peaks(1,4)-Bounding_Box_peaks(1,1);
            d3=Bounding_Box_peaks(1,6)-Bounding_Box_peaks(1,3);
    
            if d1>d2 && d1>d3
                [a_rot,idx]=sort(a_rot);
                b_rot=b_rot(idx);
                c_rot=c_rot(idx);
                d=d(idx);
            elseif d2>d1 && d2>d3
                [b_rot,idx]=sort(b_rot);
                a_rot=a_rot(idx);
                c_rot=c_rot(idx);
                d=d(idx);
            elseif d3>d1 && d3>d2
                [c_rot,idx]=sort(c_rot);
                a_rot=a_rot(idx);
                b_rot=b_rot(idx);
                d=d(idx); 
            end
                    
            % fit plane through fiber peaks - LeastSquares method
            % for advanced fitting (e.g. Gaussian fitting) https://de.mathworks.com/help/optim/ug/nonlinear-least-squares-problem-based-basics.html        
            [surffit, ~] = Vlad_fit_surface_to_fiber_plane(a_rot, b_rot, c_rot, d);
            surface_coeff=coeffvalues(surffit); % z = (1) + (2)*x + (3)*y
    
            % generate plane for drawing it
            [a_sorted,idx_a]=sort(a_rot);
            b_sorted=b_rot(idx_a);
    
            x_temp=unique(a_sorted);
            y_temp=unique(b_sorted);
    
            x_center=x_temp(floor(numel(x_temp)/2));         % center of the plane for drawing the vectors of the plane
            y_center=y_temp(floor(numel(y_temp)/2));
            z_center=surface_coeff(1) + surface_coeff(2)*x_center + surface_coeff(3)*y_center;
            
            [x_temp,y_temp]=meshgrid(x_temp,y_temp);
            z_temp=surface_coeff(1) + surface_coeff(2)*x_temp + surface_coeff(3)*y_temp;
            
    %         % find the vectors of the plane (1 normal to the plane and 2
    %         % contained in the plane and normal to the other 2)
    %         [x_normal,y_normal,z_normal]=surfnorm(x_temp,y_temp,z_temp);
    %         blue=[x_normal(2,1),y_normal(2,1),z_normal(2,1)]; blue=blue/norm(blue);     % normal vector 
        
            % choose 3 points on the plane
            point_1=[x_temp(1,1), y_temp(1,1), z_temp(1,1)];
            point_2=[x_temp(1,end), y_temp(1,end), z_temp(1,end)];
            point_3=[x_temp(end,1), y_temp(end,1), z_temp(end,1)];
    
            % find two vectors contained on the plane
            vector_1=point_3-point_1;
            vector_2=point_2-point_1;
    
            % compute the normal vector to these vectors and normalize it
            blue = cross(vector_2,vector_1); blue=blue/norm(blue);
    
            if p.plot==1
            % plot the surface which approximates the fiber peaks
            surf(x_temp(1:50:end,1:50:end),y_temp(1:50:end,1:50:end),z_temp(1:50:end,1:50:end),'EdgeColor','none','FaceColor','k','FaceAlpha',0.4)
            clim([min(p.levellist), max(p.levellist)])
            end
    
            % project points onto containing plane
            for jj=1:numel(a)
                [x_proj(jj),y_proj(jj),z_proj(jj)]=Vlad_project_point_on_plane(a_rot(jj),b_rot(jj),c_rot(jj),[blue(1);blue(2);blue(3)],x_center,y_center,z_center);
            end
    
            % green vector aligned with the longest distance between two points
            % on the plane
            red=[max(x_proj)-min(x_proj),y_proj(x_proj==max(x_proj))-y_proj(x_proj==min(x_proj)),z_proj(x_proj==max(x_proj))-z_proj(x_proj==min(x_proj))];
            red=red/norm(red); 
            % red vector generated from the cross product
            if or(size(red,2)==3,size(red,2)==3)
            else 
                blue=[1,0,0];
                red=[0,1,0];
            end
    
            green=-cross(red,blue);
            if p.plot==1
            % plot these vectors of the plane
            quiver3(x_center,y_center,z_center,blue(1),blue(2),blue(3),10,'b','LineWidth',3,'MaxHeadSize',100)
            quiver3(x_center,y_center,z_center,red(1),red(2),red(3),10,'r','LineWidth',3,'MaxHeadSize',100)
            quiver3(x_center,y_center,z_center,green(1),green(2),green(3),10,'g','LineWidth',3,'MaxHeadSize',100)
    
    
            % plot projected points
            scatter3(x_proj,y_proj,z_proj,10,'MarkerFaceColor',[255 214 10]/255,'MarkerEdgeColor','none',...
                'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
            end
            % rotate projected points to x,y plane
            R=[x_proj;y_proj;z_proj];
            %Rm=quat2rotm(quaternion([green',blue',red'],'rotmat','frame'));
    
            % Create quaternions for the two coordinate systems
            q1 = quaternion(rotm2quat([red;green;blue]));
            q2 = quaternion(rotm2quat([1 0 0; 0 1 0; 0 0 1]));
            % Compute the relative quaternion rotation
            qrel = q2 * conj(q1);
            % Convert the relative quaternion rotation to a rotation matrix
            Rm = rotmat(qrel,'frame');
    
            R_prime=Rm*R;
            
            if p.plot==1
            % plot rotated projected points
            %scatter3(R_prime(1,:),R_prime(2,:),zeros(size(R_prime(3,:))),10,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','none',...
            %    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
            scatter3(R_prime(1,:),R_prime(2,:),R_prime(3,:),10,'MarkerFaceColor',[48 209 88]/255,'MarkerEdgeColor','none',...
                'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
            end
            e=R_prime(1,:);
            f=R_prime(2,:);
            % sort the points as above
            if numel(unique(e)) >= numel(unique(f))
                [e,idx]=sort(e);
                f=f(idx);
            else
                [f,idx]=sort(f);
                e=e(idx);
            end
    %         % fit rotated projected points by polynomial (like Mobin)
    %         % compute the curviliniar coordinate
    %         s=zeros(size(e));
    %         for i = 2:length(e)
    %             s(i) = s(i-1) + sqrt( (e(i)-e(i-1))^2 + (f(i)-f(i-1))^2 );
    %         end
    %         s_resolution=100;
    %         % Set up fittype and options.
    %         ft = fittype( p.fitting_type );
    %         opts = fitoptions( 'Method', 'LinearLeastSquares' );
    % 
    %         % x direction
    %         % Fit: 'untitled fit 1'.
    %         [xData, yData, weights] = prepareCurveData( s, e, d' );
    %         if p.use_weights
    %             opts.Weights = weights;
    %         end
    %         
    %         % Fit model to data.
    %         [fitresult_x, gof_x] = fit( xData, yData, ft, opts );
    %         
    %         px_prime=fitresult_x(linspace(0,s(end),s_resolution));
    % 
    %         % y direction
    %         % Fit: 'untitled fit 1'.
    %         [xData, yData, weights] = prepareCurveData( s, f, d' );
    %         if p.use_weights
    %             opts.Weights = weights;
    %         end
    %         
    %         % Fit model to data.
    %         [fitresult_y, ~] = fit( xData, yData, ft, opts );
    %         
    %         py_prime=fitresult_y(linspace(0,s(end),s_resolution));
    
            % fit rotated projected points by polynomial (normal)
            [xData, yData, weights] = prepareCurveData( e, f, d' );
            ft = fittype( p.fitting_type_fibre );
            opts = fitoptions( 'Method', 'LinearLeastSquares' );
    
            if p.use_weights
                opts.Weights = weights;
            end
            
            % Fit polynomial to ro
            [fitresult_y,~] = fit( xData, yData, ft, opts );
            py_prime=fitresult_y(e);
            px_prime=e';
    
            if p.plot==1
            % plot fitted polynomial to the rotated projected points
            plot3(px_prime,py_prime,ones(size(px_prime))*R_prime(3,1),'k-','LineWidth',2)
            end
            % rotate projected points back to plane
            PP_prime=[px_prime';py_prime';ones(size(px_prime,1),1)'*R_prime(3,1)];
            PP=Rm\PP_prime;
    
            px=PP(1,:)';
            py=PP(2,:)';
            pz=PP(3,:)';
    
            if p.plot==1
            % plot the polynomial on the plane
            plot3(px,py,pz,'k-','LineWidth',2)
            end
            % rotate polynomial back to the original fiber
            PX_PY_PZ_prime=[px';py';pz'];
            PX_PY_PZ=Rm_peaks\PX_PY_PZ_prime;
    
            px=PX_PY_PZ(1,:)';
            py=PX_PY_PZ(2,:)';
            pz=PX_PY_PZ(3,:)';
            if p.plot==1
            % plot the polynomial on the original fiber
            plot3(px,py,pz,'k-','LineWidth',2)
            end
            %close all


%%%%%%%%%%%%%%% second approximation %%%%%%%

% make figure
        if p.plot==1
        figure(7); clf; hold on; grid on; box on; daspect([1 1 1]); view(3)
        set(gca,'TickLabelInterpreter','latex')
        xlabel('$x$','Interpreter','latex')
        ylabel('$y$','Interpreter','latex')
        zlabel('$z$','Interpreter','latex')
        % lab coordinate system
        quiver3(0,0,0,1,0,0,5,'filled','LineWidth',2,'Color','r','MaxHeadSize',10)
        quiver3(0,0,0,0,1,0,5,'filled','LineWidth',2,'Color','g','MaxHeadSize',10)
        quiver3(0,0,0,0,0,1,5,'filled','LineWidth',2,'Color','b','MaxHeadSize',10)

        % plot fiber intensity with isosurface
        [A, B, C]=meshgrid(1:size(Fiber,1),1:size(Fiber,2),1:size(Fiber,3));
        for i=1:length(p.levellist)
            level=p.levellist(i);
            pat=patch(isosurface(permute(B,[2 1 3]),permute(A,[2 1 3]),permute(C,[2 1 3]),Fiber,level));
            pat.FaceVertexCData=level;
            pat.FaceColor='flat';
            pat.EdgeColor='none';
            pat.FaceAlpha=(p.facealphalist(i));
        end
        colormap(hsv(numel(p.levellist)))
        clb=colorbar;
        clb.Limits=[min(p.levellist), max(p.levellist)];
        clb.Location="southoutside";

        % plot fiber peaks
        scatter3(a,b,c,10,'MarkerFaceColor',[255 69 58]/255,'MarkerEdgeColor','none',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
        end


        [red_peaks,green_peaks,blue_peaks,Centroid_refined_peaks]=Vlad_find_orientation_vectors_angles(px',py',pz','tensor');

       
        red_peaks=red_peaks';
        green_peaks=green_peaks';
        blue_peaks=blue_peaks';

        loc_x = Centroid_refined_peaks(1);
        loc_y = Centroid_refined_peaks(2);
        loc_z = Centroid_refined_peaks(3);


        if p.plot==1
        quiver3(loc_x,loc_y,loc_z,red_peaks(1),red_peaks(2),red_peaks(3),10,'r','LineWidth',3,'MaxHeadSize',100)
        quiver3(loc_x,loc_y,loc_z,green_peaks(1),green_peaks(2),green_peaks(3),10,'g','LineWidth',3,'MaxHeadSize',100)
        quiver3(loc_x,loc_y,loc_z,blue_peaks(1),blue_peaks(2),blue_peaks(3),10,'b','LineWidth',3,'MaxHeadSize',100)
        end
        % align the fiber peaks (locations x,y,z are stored in a,b,c) with the lab coordinate system
        A_B_C=[a';b';c'];

        Rm_peaks=[red_peaks green_peaks blue_peaks]';

        A_B_C_prime=Rm_peaks*A_B_C;
        a_rot=A_B_C_prime(1,:)';
        b_rot=A_B_C_prime(2,:)';
        c_rot=A_B_C_prime(3,:)';
        if p.plot==1
        % plot rotated fiber peaks
        scatter3(a_rot,b_rot,c_rot,10,'MarkerFaceColor',[255 159 10]/255,'MarkerEdgeColor','none',...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
        end
        % sort the locations of the fiber peaks (a,b,c) based on bounding
        % box sizes
%         d1=Bounding_Box(1,2)-Bounding_Box(1,1);
%         d2=Bounding_Box(2,2)-Bounding_Box(2,1);
%         d3=Bounding_Box(3,2)-Bounding_Box(3,1);

        d1=Bounding_Box_peaks(1,5)-Bounding_Box_peaks(1,2);
        d2=Bounding_Box_peaks(1,4)-Bounding_Box_peaks(1,1);
        d3=Bounding_Box_peaks(1,6)-Bounding_Box_peaks(1,3);

        if d1>d2 && d1>d3
            [a_rot,idx]=sort(a_rot);
            b_rot=b_rot(idx);
            c_rot=c_rot(idx);
            d=d(idx);
        elseif d2>d1 && d2>d3
            [b_rot,idx]=sort(b_rot);
            a_rot=a_rot(idx);
            c_rot=c_rot(idx);
            d=d(idx);
        elseif d3>d1 && d3>d2
            [c_rot,idx]=sort(c_rot);
            a_rot=a_rot(idx);
            b_rot=b_rot(idx);
            d=d(idx); 
        end
                
        % fit plane through fiber peaks - LeastSquares method
        % for advanced fitting (e.g. Gaussian fitting) https://de.mathworks.com/help/optim/ug/nonlinear-least-squares-problem-based-basics.html        
        [surffit, ~] = Vlad_fit_surface_to_fiber_plane(a_rot, b_rot, c_rot, d);
        surface_coeff=coeffvalues(surffit); % z = (1) + (2)*x + (3)*y

        % generate plane for drawing it
        [a_sorted,idx_a]=sort(a_rot);
        b_sorted=b_rot(idx_a);

        x_temp=unique(a_sorted);
        y_temp=unique(b_sorted);

        x_center=x_temp(floor(numel(x_temp)/2));         % center of the plane for drawing the vectors of the plane
        y_center=y_temp(floor(numel(y_temp)/2));
        z_center=surface_coeff(1) + surface_coeff(2)*x_center + surface_coeff(3)*y_center;
        
        [x_temp,y_temp]=meshgrid(x_temp,y_temp);
        z_temp=surface_coeff(1) + surface_coeff(2)*x_temp + surface_coeff(3)*y_temp;
        
%         % find the vectors of the plane (1 normal to the plane and 2
%         % contained in the plane and normal to the other 2)
%         [x_normal,y_normal,z_normal]=surfnorm(x_temp,y_temp,z_temp);
%         blue=[x_normal(2,1),y_normal(2,1),z_normal(2,1)]; blue=blue/norm(blue);     % normal vector 
    
        % choose 3 points on the plane
        point_1=[x_temp(1,1), y_temp(1,1), z_temp(1,1)];
        point_2=[x_temp(1,end), y_temp(1,end), z_temp(1,end)];
        point_3=[x_temp(end,1), y_temp(end,1), z_temp(end,1)];

        % find two vectors contained on the plane
        vector_1=point_3-point_1;
        vector_2=point_2-point_1;

        % compute the normal vector to these vectors and normalize it
        blue = cross(vector_2,vector_1); blue=blue/norm(blue);

        if p.plot==1
        % plot the surface which approximates the fiber peaks
        surf(x_temp(1:50:end,1:50:end),y_temp(1:50:end,1:50:end),z_temp(1:50:end,1:50:end),'EdgeColor','none','FaceColor','k','FaceAlpha',0.4)
        clim([min(p.levellist), max(p.levellist)])
        end

        % project points onto containing plane
        for jj=1:numel(a)
            [x_proj(jj),y_proj(jj),z_proj(jj)]=Vlad_project_point_on_plane(a_rot(jj),b_rot(jj),c_rot(jj),[blue(1);blue(2);blue(3)],x_center,y_center,z_center);
        end

        % green vector aligned with the longest distance between two points
        % on the plane
        red=[max(x_proj)-min(x_proj),y_proj(x_proj==max(x_proj))-y_proj(x_proj==min(x_proj)),z_proj(x_proj==max(x_proj))-z_proj(x_proj==min(x_proj))];
        red=red/norm(red); 
        % red vector generated from the cross product
        if or(size(red,2)==3,size(red,2)==3)
        else 
            blue=[1,0,0];
            red=[0,1,0];
        end

        green=-cross(red,blue);
        if p.plot==1
        % plot these vectors of the plane
        quiver3(x_center,y_center,z_center,blue(1),blue(2),blue(3),10,'b','LineWidth',3,'MaxHeadSize',100)
        quiver3(x_center,y_center,z_center,red(1),red(2),red(3),10,'r','LineWidth',3,'MaxHeadSize',100)
        quiver3(x_center,y_center,z_center,green(1),green(2),green(3),10,'g','LineWidth',3,'MaxHeadSize',100)


        % plot projected points
        scatter3(x_proj,y_proj,z_proj,10,'MarkerFaceColor',[255 214 10]/255,'MarkerEdgeColor','none',...
            'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
        end
        % rotate projected points to x,y plane
        R=[x_proj;y_proj;z_proj];
        %Rm=quat2rotm(quaternion([green',blue',red'],'rotmat','frame'));

        % Create quaternions for the two coordinate systems
        q1 = quaternion(rotm2quat([red;green;blue]));
        q2 = quaternion(rotm2quat([1 0 0; 0 1 0; 0 0 1]));
        % Compute the relative quaternion rotation
        qrel = q2 * conj(q1);
        % Convert the relative quaternion rotation to a rotation matrix
        Rm = rotmat(qrel,'frame');

        R_prime=Rm*R;
        
        if p.plot==1
        % plot rotated projected points
        %scatter3(R_prime(1,:),R_prime(2,:),zeros(size(R_prime(3,:))),10,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','none',...
        %    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
        scatter3(R_prime(1,:),R_prime(2,:),R_prime(3,:),10,'MarkerFaceColor',[48 209 88]/255,'MarkerEdgeColor','none',...
            'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
        end
        e=R_prime(1,:);
        f=R_prime(2,:);
        % sort the points as above
        if numel(unique(e)) >= numel(unique(f))
            [e,idx]=sort(e);
            f=f(idx);
        else
            [f,idx]=sort(f);
            e=e(idx);
        end
%         % fit rotated projected points by polynomial (like Mobin)
%         % compute the curviliniar coordinate
%         s=zeros(size(e));
%         for i = 2:length(e)
%             s(i) = s(i-1) + sqrt( (e(i)-e(i-1))^2 + (f(i)-f(i-1))^2 );
%         end
%         s_resolution=100;
%         % Set up fittype and options.
%         ft = fittype( p.fitting_type );
%         opts = fitoptions( 'Method', 'LinearLeastSquares' );
% 
%         % x direction
%         % Fit: 'untitled fit 1'.
%         [xData, yData, weights] = prepareCurveData( s, e, d' );
%         if p.use_weights
%             opts.Weights = weights;
%         end
%         
%         % Fit model to data.
%         [fitresult_x, gof_x] = fit( xData, yData, ft, opts );
%         
%         px_prime=fitresult_x(linspace(0,s(end),s_resolution));
% 
%         % y direction
%         % Fit: 'untitled fit 1'.
%         [xData, yData, weights] = prepareCurveData( s, f, d' );
%         if p.use_weights
%             opts.Weights = weights;
%         end
%         
%         % Fit model to data.
%         [fitresult_y, ~] = fit( xData, yData, ft, opts );
%         
%         py_prime=fitresult_y(linspace(0,s(end),s_resolution));

        % fit rotated projected points by polynomial (normal)
        [xData, yData, weights] = prepareCurveData( e, f, d' );
        ft = fittype( p.fitting_type );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );

        if p.use_weights
            opts.Weights = weights;
        end
        
        % Fit polynomial to ro
        [fitresult_y,~] = fit( xData, yData, ft, opts );
        py_prime=fitresult_y(e);
        px_prime=e';

        if p.plot==1
        % plot fitted polynomial to the rotated projected points
        plot3(px_prime,py_prime,ones(size(px_prime))*R_prime(3,1),'k-','LineWidth',2)
        end
        % rotate projected points back to plane
        PP_prime=[px_prime';py_prime';ones(size(px_prime,1),1)'*R_prime(3,1)];
        PP=Rm\PP_prime;

        px=PP(1,:)';
        py=PP(2,:)';
        pz=PP(3,:)';

        if p.plot==1
        % plot the polynomial on the plane
        plot3(px,py,pz,'k-','LineWidth',2)
        end
        % rotate polynomial back to the original fiber
        PX_PY_PZ_prime=[px';py';pz'];
        PX_PY_PZ=Rm_peaks\PX_PY_PZ_prime;

        px=PX_PY_PZ(1,:)';
        py=PX_PY_PZ(2,:)';
        pz=PX_PY_PZ(3,:)';
        if p.plot==1
        % plot the polynomial on the original fiber
        plot3(px,py,pz,'k-','LineWidth',2)
        end
        %close all

        %%
end


%% find reference frame of the fitted polynomial

[Curvature,Length]=Vlad_compute_curvature_length(px,py,pz);
Length_effective=norm([px(1) py(1) pz(1)]-[px(end) py(end) pz(end)]);

[red,green,blue,Centroid_refined]=Vlad_find_orientation_vectors_angles(px',py',pz','tensor');
% if Vlad_angle_between_two_vectors(red,[1 0 0])>90
%     red=-red; %green=-green;
% end
% blue=cross(red,green);
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