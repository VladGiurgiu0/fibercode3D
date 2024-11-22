function [x_proj,y_proj,z_proj]=Vlad_project_point_on_plane(p_x,p_y,p_z,n,q_x,q_y,q_z)
%%% referenceL: https://math.libretexts.org/Courses/Monroe_Community_College/MTH_212_Calculus_III/Chapter_11%3A_Vectors_and_the_Geometry_of_Space/11.5%3A_Equations_of_Lines_and_Planes_in_Space

%%-----------------------------------------------------------------------%%
%%% input
% coordinates of point in 3D space: p_x, p_y, p_z
%                               e.g.: 5,5,5
% normal of plane on which point is to be projected: n 
%                                           e.g.: [1;1;1]
% coordinates of point on plane through which n goes: q_x,q_y,q_z
%                                           e.g.: 0,0,0
%%-----------------------------------------------------------------------%%
%%% output
% coordinates of point projected onto plane: x_proj,y_proj,z_proj
%                                           e.g.: 3,3,3
%%-----------------------------------------------------------------------%%

%% computation
P=[p_x; p_y; p_z];              % vector of point P
Q=[q_x; q_y; q_z];              % vector of point Q (base of vector n)
QP=P-Q;                         % vector between Q and P
d=dot(QP,n)/norm(n);       % distance between point and plane

RP=n*d;                         % vector between R and P
R=(QP-RP)+Q;                        % projected point R on the plane

x_proj=R(1);
y_proj=R(2);
z_proj=R(3);

end