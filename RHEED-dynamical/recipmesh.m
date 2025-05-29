function [Bm, Sz, Bmn] = recipmesh(Na, Nb, recip, S0, K0, K0x, K0y, K0z)
% recipmesh
% Create a mesh of reciprocal lattice vectors centered at the point S0 in
% reciprocal space.
%
% Inputs:
% Na        Number of reciprocal lattice a* vectors = -Na, Na+1, ... Na
% Nb        Number of reciprocal lattice b* vectors = -Nb, Nb+1, ... Nb
% recip     3x3 array of reciprocal lattice vectors, in Angstroms^-1
% S0        Triplet of scattering vector offsets (x,y,z), in Angstroms^-1
% K0        Incident wavevector magnitude, in Angstroms^-1
% K0x, K0y, K0z     (x,y,z) components of incident wavevector, in 
%                   Angstroms^-1
%
% Outputs:
% Bm        (2Na+1)x(2Nb+1)x2 array of reciprocal lattice vectors, in
%           Cartesian (x,y) coordinates.
% Sz        (2Na+1)x(2Nb+1) array of reciprocal space vectors, in Cartesian
%           z coordinates.
% Bmn       (2Na+1)x(2Nb+1)x2x2 array of reciprocal lattice vectors (a*,b*)

% Initialize arrays.
Bmn = zeros(2*Na+1, 2*Nb+1, 2, 2);
Bm = zeros(2*Na+1, 2*Nb+1, 2);
Sz = zeros(2*Na+1, 2*Nb+1);
% Sz = zeros(2*Na+1, 2*Nb+1, 2);          % Two solution branches for Sz

% Loop on each mesh point.
for m=-Na:Na
    for n=-Nb:Nb
        % Array indices
        i = m+Na+1;
        j = n+Nb+1;

        % Assign the reciprocal lattice vectors in-plane.
        Bmn(i,j,1,:) = m*recip(1,1:2);
        Bmn(i,j,2,:) = n*recip(2,1:2);

        % Convert to Cartesian coordinates.
        Bx = Bmn(i,j,1,1) + Bmn(i,j,2,1) + S0(1);
        By = Bmn(i,j,1,2) + Bmn(i,j,2,2) + S0(2);
        Bm(i,j,:) = [Bx, By];

        % Calculate the z-component of the scattering vector.
        Szsq = K0^2 - (Bx + K0x).^2 - (By + K0y).^2;
        Sz(i,j) = sqrt(Szsq) - K0z;         % Angstroms^-1
        % Sz(i,j,1) = sqrt(Szsq) - K0z;         % Angstroms^-1
        % Sz(i,j,2) = -sqrt(Szsq) - K0z;         % Angstroms^-1
        % Sz(i,j) = sqrt(Szsq) - K0z + S0(3);         % Angstroms^-1
    end
end
end