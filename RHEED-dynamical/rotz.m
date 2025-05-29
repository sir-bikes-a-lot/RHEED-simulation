function [Rt] = rotz(R, theta)
% rotz
% Rotate the input vector (in Cartesian xyz coordinates) by angle theta
% along the z-axis, according to the right-hand hule.
% Inputs:
% R         Vector of Cartesian coordinates (x,y,z)
% theta     Angle of rotation around z-axis, in radians
% Outputs:
% Rt        Rotated vector of Cartesian coordinates (x,y,z)

A = [cos(theta), -sin(theta), 0;...
     sin(theta), cos(theta), 0;...
     0, 0, 1];

Rt = A*R;

end