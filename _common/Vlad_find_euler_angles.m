function [Orientation]=Vlad_find_euler_angles(red,green,blue)
use_quaternion=0;
Tait_Bryan=0;

if size(red,1)==1
    red=red';
    green=green';
    blue=blue';
end

if use_quaternion==1
    Orientation=rad2deg(quat2eul(quaternion([red';green';blue'],'rotmat','frame')));
else

ref.V1=red;
ref.V2=green;
ref.V3=blue;
debug=0;
% Version 1 - 26/04/2020 - Marco De Paoli
% Version 2 - 25/05/2020 - Marco De Paoli
% Version 3 - 14/07/2020 - Marco De Paoli
% Version 4 - 20/02/2021 - Mobin Alipour

% %%%   The rotation matrices and the correspondent angles of rotation are
%       here computed.Rotation is computed with respect to the lab
%       reference frame (e1,e2,e3).
% %%%   Angles computed in a different way and matrix used to
%       compute the angular velocities
% %%%   Angles defined as in Marchioli (2010)
% %%%   Angle gamma added to the code

% % Global reference frame
% rot.e1=[1 ; 0 ; 0];
% rot.e2=[0 ; 1 ; 0];
% rot.e3=[0 ; 0 ; 1];

% if debug==1
%     quiver3(fp.B(1),fp.B(2),fp.B(3),rot.e1(1),rot.e1(2),rot.e1(3),0.5,...
%         'Linewidth',2,'color','b')
%     quiver3(fp.B(1),fp.B(2),fp.B(3),rot.e2(1),rot.e2(2),rot.e2(3),0.5,...
%         'Linewidth',2,'color','b')
%     quiver3(fp.B(1),fp.B(2),fp.B(3),rot.e3(1),rot.e3(2),rot.e3(3),0.5,...
%         'Linewidth',2,'color','b')
% end

% The rotation matrix required to obtain the eigen-vectors of the inertia
% tensor of the fiber starting from the reference frame of the system is
% given by the matrix that has with column the eigen vectors
R_fin=[ref.V1 , ref.V2 , ref.V3];


% Compute Euler angles : Possible rotations
%   "ZYX" (default) ? The order of rotation angles is z-axis, y-axis, x-axis.
%   "ZYZ" ? The order of rotation angles is z-axis, y-axis, z-axis.
%   "XYZ" ? The order of rotation angles is x-axis, y-axis, z-axis.
% Before computing the angles, be sure the rotation matrix is orthogonal
if abs(sum(sum(R_fin*R_fin'))-3)>1E-8
    disp('-------- Rotation matrix not orthogonal!')
    out.eul = rotm2eul(R_fin);
    if Tait_Bryan==1
        out.eul=rot2taitbryan(R_fin);
    end
    if debug==1  ; disp(rotm2eul(R_fin)*180/pi) ; end
else
    out.eul = real(rotm2eul(R_fin));
    if debug==1  ; disp(real(rotm2eul(R_fin)*180/pi)) ; end
end

% if debug==1
%     
%     % Check - from global to local reference frame
%     rot.e11=R_fin*rot.e1;
%     rot.e21=R_fin*rot.e2;
%     rot.e31=R_fin*rot.e3;
%     
%     quiver3(fp.B(1)+0.02,fp.B(2)+0.02,fp.B(3)+0.02,...
%         rot.e11(1),rot.e11(2),rot.e11(3),0.5,...
%         'Linewidth',2,'color','r')
%     quiver3(fp.B(1)+0.02,fp.B(2)+0.02,fp.B(3)+0.02,...
%         rot.e21(1),rot.e21(2),rot.e21(3),0.5,...
%         'Linewidth',2,'color','r')
%     quiver3(fp.B(1)+0.02,fp.B(2)+0.02,fp.B(3)+0.02,...
%         rot.e31(1),rot.e31(2),rot.e31(3),0.5,...
%         'Linewidth',2,'color','r')
%     
% end

% % Find fiber orientation with respect to different convenctions
% % Marchioli, Fantoni and Soldati (2010)
% % theta_x
% [out.theta,R1]=find_rotation_matrix(rot.e1,ref.V1);
% if isnan(R1) ; R1=-eye(3) ; out.theta=pi ; end
% % theta_y
% [out.phi,R2]=find_rotation_matrix(rot.e3,ref.V1);
% if isnan(R2) ; R2=eye(3) ; out.phi=pi  ; end
% % theta_z (thetha_y in our case)
% [out.psi,R3]=find_rotation_matrix(rot.e2,ref.V1);
% if isnan(R3) ; R3=eye(3) ; out.psi=pi  ; end
% 
% % gamma_y (y' with y'')
% [out.gamma,R4]=find_rotation_matrix(rot.e2,ref.V2);
% if isnan(R4) ; R4=eye(3) ; out.gamma=pi  ; end
% % gamma_x (y' with x'')
% [out.gamma_x,R5]=find_rotation_matrix(rot.e1,ref.V2);
% if isnan(R5) ; R5=eye(3) ; out.gamma_x=pi  ; end
% % gamma_z (y' with z'')
% [out.gamma_z,R6]=find_rotation_matrix(rot.e3,ref.V2);
% if isnan(R6) ; R6=eye(3) ; out.gamma_z=pi  ; end

% R_rot=R3*R2*R1;

Orientation=rad2deg([out.eul(1),out.eul(2),out.eul(3)]);
end
end

function [phi, theta, psi] = rot2taitbryan(R)
% Extract Tait-Bryan angles from rotation matrix
phi = atan2(R(2,1), R(1,1));
theta = -asin(R(3,1));
psi = atan2(R(3,2), R(3,3));
end
