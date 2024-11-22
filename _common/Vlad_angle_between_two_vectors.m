function [angle]=Vlad_angle_between_two_vectors(a,b)
        angle=atan2d(norm(cross(a,b)),dot(a,b));
end