function [red,green,blue,Centroid]=Vlad_find_orientation_vectors_angles(px,py,pz,Method)
%%% test case
%    s=linspace(0,1,20);
%    px=s.^2; py=s; pz=s;
    
    debug=0;
    
    switch Method
        case "mid-point"
    
            idx_mid=floor(numel(s)/2);
            [a,b,c,d]=Plane_3Points([px(1);py(1);pz(1)],[px(idx_mid);py(idx_mid);pz(idx_mid)],[px(end);py(end);pz(end)]);
            
            mid_point=[px(idx_mid),py(idx_mid),pz(idx_mid)];
            
            x=linspace(0,1,3);
            y=linspace(0,1,3);
            [x,y]=meshgrid(x,y);
            z=(a*x + b*y + d)/(-c);
            [a,b,c]=surfnorm(x,y,z);
            
            surf(x,y,z);
            blue=[a(1,1),b(1,1),c(1,1)];
            
            
            daspect([1 1 1])
            
            px_prime=gradient(px);
            py_prime=gradient(py);
            pz_prime=gradient(pz);
            
            red=[px_prime(idx_mid),py_prime(idx_mid),pz_prime(idx_mid)];
            red=red/norm(red);
            
            
            green=cross(red,blue);
            green=green/norm(green);
    
            if debug==1
                plot3(px,py,pz,'LineWidth',5); box on; grid on; hold on
                quiver3(mid_point(1),mid_point(2),mid_point(3),red(1),red(2),red(3),'LineWidth',3,'Color','red');
                quiver3(mid_point(1),mid_point(2),mid_point(3),green(1),green(2),green(3),'LineWidth',3,'Color','green');
                quiver3(mid_point(1),mid_point(2),mid_point(3),blue(1),blue(2),blue(3),'LineWidth',3,'Color','blue')
            end
    
            main.xs=px;
            main.ys=py;
            main.zs=pz;
    
            nss=numel(main.xs);
            
            ls=sqrt(...
                (main.xs(2:end)-main.xs(1:end-1)).^2+...
                (main.ys(2:end)-main.ys(1:end-1)).^2+...
                (main.zs(2:end)-main.zs(1:end-1)).^2);
            
            % find mass of each segment
            mi=zeros(1,nss);
            mi(2:nss-1)=(ls(2:end)+ls(1:end-1))/2;
            mi(1)=ls(1)/2; mi(nss)=ls(nss-1)/2;
            
            % find center of mass  of the fiber
            %    mi is the mass of each segment
            
            ref.xg=sum(main.xs.*mi)/sum(mi);
            ref.yg=sum(main.ys.*mi)/sum(mi);
            ref.zg=sum(main.zs.*mi)/sum(mi);
            
        case "tensor"
            
            main.xs=px;
            main.ys=py;
            main.zs=pz;
            
            % Version 1 - 06/04/2020 - Marco De Paoli
            % Version 2 - 12/08/2020 - Marco De Paoli
            % Inertia tensor was calculated with repect to extrema, in this version 
            % it has been changed to be computed with respect to Center of Mass 
            
            % CASE 2 - x aligned with first eigenvector, associated to the smallest
            %          eigenvalue of the inertia tensor matrix
            
            nss=numel(main.xs);
            
            % find length of the segments
            ls=sqrt(...
                (main.xs(2:end)-main.xs(1:end-1)).^2+...
                (main.ys(2:end)-main.ys(1:end-1)).^2+...
                (main.zs(2:end)-main.zs(1:end-1)).^2);
            
            % find mass of each segment
            mi=zeros(1,nss);
            mi(2:nss-1)=(ls(2:end)+ls(1:end-1))/2;
            mi(1)=ls(1)/2; mi(nss)=ls(nss-1)/2;
            
            % find center of mass  of the fiber
            %    mi is the mass of each segment
            
            ref.xg=sum(main.xs.*mi)/sum(mi);
            ref.yg=sum(main.ys.*mi)/sum(mi);
            ref.zg=sum(main.zs.*mi)/sum(mi);
            
            % compute with respect to one Center of Mass (A)
            main.xs1=main.xs-ref.xg;
            main.ys1=main.ys-ref.yg;
            main.zs1=main.zs-ref.zg;
            
            % compute with respect to one extrema (A)
            % main.xs1=main.xs-fp.A(1);
            % main.ys1=main.ys-fp.A(2);
            % main.zs1=main.zs-fp.A(3);
            
            A=zeros(3);
            % diagonal part
            A(1,1)=sum(mi.*(main.ys1.^2+main.zs1.^2));
            A(2,2)=sum(mi.*(main.xs1.^2+main.zs1.^2));
            A(3,3)=sum(mi.*(main.xs1.^2+main.ys1.^2));
            % upper part
            A(1,2)=-sum(mi.*(main.xs1.*main.ys1));
            A(1,3)=-sum(mi.*(main.xs1.*main.zs1));
            A(2,3)=-sum(mi.*(main.ys1.*main.zs1));
            A(2,1)=A(1,2); A(3,1)=A(1,3); A(3,2)=A(2,3);
            
            % %%% find the eigenvalues (principal inertia moments, stored in D) and the
            %     eigenvectors (directions of inertia axis, stored in columns of V)
            [V,~] = eig(A); % already normalised
            
            % [V,D] = eig(A); % use this to compare D with Dtheo
            ref.V1=V(:,1); ref.V2=V(:,2); ref.V3=V(:,3);
    
            if debug==1
                quiver3(ref.xg,ref.yg,ref.zg,ref.V1(1),ref.V1(2),ref.V1(3),'Linewidth',3,'Color','magenta')
                quiver3(ref.xg,ref.yg,ref.zg,ref.V2(1),ref.V2(2),ref.V2(3),'Linewidth',3,'Color','yellow')
                quiver3(ref.xg,ref.yg,ref.zg,ref.V3(1),ref.V3(2),ref.V3(3),'Linewidth',3,'Color','cyan')
            end
            
            red=[ref.V1(1),ref.V1(2),ref.V1(3)];
            green=[ref.V2(1),ref.V2(2),ref.V2(3)];
            blue=[ref.V3(1),ref.V3(2),ref.V3(3)];
    
    end
    
    Centroid=[ref.xg,ref.yg,ref.zg];
    if debug==1
        scatter3(Centroid(1),Centroid(2),Centroid(3),5000,'red','.','LineWidth',5)
    end
end
