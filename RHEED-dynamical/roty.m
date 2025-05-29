function [Rt] = roty(R, theta)
% roty
% Rotate the input vector (in Cartesian xyz coordinates) by angle theta
% along the y-axis, according to the right-hand hule.
% Inputs:
% R         Vector of Cartesian coordinates (x,y,z)
% theta     Angle of rotation around y-axis, in radians
% Outputs:
% Rt        Rotated vector of Cartesian coordinates (x,y,z)

A = [cos(theta), 0, sin(theta);...
     0, 1, 0;...
     -sin(theta), 0, cos(theta)];

Rt = A*R;

end