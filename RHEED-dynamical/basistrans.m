function [Rt] = basistrans(R, gamma, psi)
% basistrans
% Transform from lattice basis vectors (a,b,c) to Cartesian basis vectors
% (x,y,z).
% Currently only set up for transformation in the (a,b) plane.
% Inputs:
% R         Vector of (a,b,c) coordinates
% psi       Azimuth angle in (a,b) plane, in radians
% gamma     Angle between basis vectors a and b
% Outputs:
% Rt        Vector of Cartesian coordinates (x,y,z)

A = [cos(psi), cos(gamma - psi), 0;...
     -sin(psi), sin(gamma - psi), 0;...
     0, 0, 1];

Rt = A*R;

end