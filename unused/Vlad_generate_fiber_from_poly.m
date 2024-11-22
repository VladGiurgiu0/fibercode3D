function [Fiber_poly]=Vlad_generate_fiber_from_poly(px,py,pz,sphere_radius,X,Y,Z)
%% input
% px, py, pz are the centers of each section of the fiber

% X=linspace(floor(min(px))-sphere_radius,ceil(max(px))+sphere_radius,size(px,1));
% Y=linspace(floor(min(py))-sphere_radius,ceil(max(py))+sphere_radius,size(py,1));
% Z=linspace(floor(min(pz))-sphere_radius,ceil(max(pz))+sphere_radius,size(pz,1));

[Xs,Ys,Zs]=meshgrid(X,Y,Z);

%%  Salami generator
% produce first sphere
Fiber_poly = (Xs-px(1,1)).^2 + (Ys-py(1,1)).^2 + (Zs-pz(1,1)).^2 <=    sphere_radius.^2;

for j=2:size(px,1)          % for each center produce one sphere
    % add them together
    Fiber_poly = or(Fiber_poly,((Xs-px(j,1)).^2 + (Ys-py(j,1)).^2 + (Zs-pz(j,1)).^2 <=    sphere_radius.^2));
end

% % plot for verification
% figure(1)
% fv = isosurface(Xs,Ys,Zs,Fiber_poly,0);
% patch(fv,'FaceColor',[0 0 .7],'EdgeColor',[0 0 1]); 
% daspect([1 1 1]);
% box on; grid on
