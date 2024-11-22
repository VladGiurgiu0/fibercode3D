function ref=find_refframe(main,fp)

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

% check: since V is orthogonal, it must be that : inv(V) == V'
% check: for straight and infinitely thin fiber, the moment inertia tensor
%        is defined as follows
% Mtot=sum(mi); Dtheo=Mtot*[0 0 0 ; 0 1/3*sum(ls)^2 0 ; 0 0 1/3*sum(ls)^2]

if main.debug==1
    
    % check consistency of projections:
    % Find projection of parametric curves (xs,ys,zs) on the plane
    % containing the fiber:
    %       (xx,yy,zz)=(x,y,z)-(((x,y,z)-P)\cdot n) n
    % with n=(alpha, beta, gamma) and P point of the plane
    
    ref.xxs=zeros(size(main.xs));
    ref.yys=ref.xxs;
    ref.zzs=ref.xxs;
    for i=1:numel(main.xs)
        ref.xxs(i)=main.xs(i)-...
            ((main.xs(i)-fp.A(1))*ref.V3(1)+(main.ys(i)-fp.A(2))*...
            ref.V3(2)+(main.zs(i)-fp.A(3))*ref.V3(3))*ref.V3(1);
        ref.yys(i)=main.ys(i)-...
            ((main.xs(i)-fp.A(1))*ref.V3(1)+(main.ys(i)-fp.A(2))*...
            ref.V3(2)+(main.zs(i)-fp.A(3))*ref.V3(3))*ref.V3(2);
        ref.zzs(i)=main.zs(i)-...
            ((main.xs(i)-fp.A(1))*ref.V3(1)+(main.ys(i)-fp.A(2))*...
            ref.V3(2)+(main.zs(i)-fp.A(3))*ref.V3(3))*ref.V3(3);
    end
    
    plot3(ref.xxs(1:4:end),ref.yys(1:4:end),ref.zzs(1:4:end),'ok',...
        'Markerfacecolor','y','Markersize',10)
%     
    quiver3(fp.B(1),fp.B(2),fp.B(3),ref.V1(1),ref.V1(2),ref.V1(3),0.5,...
        'Linewidth',2,'color','k')
    quiver3(fp.B(1),fp.B(2),fp.B(3),ref.V2(1),ref.V2(2),ref.V2(3),0.5,...
        'Linewidth',2,'color','k')
    quiver3(fp.B(1),fp.B(2),fp.B(3),ref.V3(1),ref.V3(2),ref.V3(3),0.5,...
        'Linewidth',2,'color','k')
    
    % plot projections to explain reconstruction
    
    plot3(main.xs,main.ys,ones(size(main.zs))*min(main.zs),'-b')
    plot3(main.xs,ones(size(main.ys))*min(main.ys),main.zs,'-b')
    plot3(ones(size(main.xs))*min(main.xs),main.ys,main.zs,'-b')
    ref.a=1/2;
    xlim([min(main.xs)-ref.a max(main.xs)+ref.a])
    ylim([min(main.ys)-ref.a max(main.ys)+ref.a])
    zlim([min(main.zs)-ref.a max(main.zs)+ref.a])
    
end


end