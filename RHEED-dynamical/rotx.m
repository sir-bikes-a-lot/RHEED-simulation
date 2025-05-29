function [Rt] = rotx(R, theta)
% rotx
% Rotate the input vector (in Cartesian xyz coordinates) by angle theta
% along the x-axis, according to the right-hand hule.
% Inputs:
% R         Vector of Cartesian coordinates (x,y,z)
% theta     Angle of rotation around x-axis, in radians
% Outputs:
% Rt        Rotated vector of Cartesian coordinates (x,y,z)

A = [1, 0, 0;...
     0, cos(theta), -sin(theta);...
     0, sin(theta), cos(theta)];

Rt = A*R;

end