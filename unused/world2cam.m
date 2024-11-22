%%-----------------------------------------------------------------------%%
%%--Convert world coordinates (U,V,W) to the camera coordinates (X,Y,Z)--%%
%%-----------------------------------------------------------------------%%

function [X,Y,Z] = world2cam(U,V,W,alfa,beta,gama,T)

%Translation matrix: 
%T contains the camera position in the global coordinates (Tx,Ty,Tz)
Tm = [1    0    0  T(1);
      0    1    0  T(2);
      0    0    1  T(3);
      0    0    0   1];
%Rotation matrix:
%alfa, beta, and gama are the rotation in relation to the x, y, and z-axis
rot = rotx(alfa)*roty(beta)*rotz(gama);
Rm = [rot(1,1) rot(1,2) rot(1,3) 0;
      rot(2,1) rot(2,2) rot(2,3) 0;
      rot(3,1) rot(3,2) rot(3,3) 0;
         0        0        0     1];
%Global coordinates vector:
%(U,V,W) is the position of a point in global coordinates
coor = [U;
        V;
        W;  
        1];
%Conversion from global to camera coordinates:
XYZ = Tm*Rm*coor;
X = XYZ(1); Y = XYZ(2); Z = XYZ(3);
 