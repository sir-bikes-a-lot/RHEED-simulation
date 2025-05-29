function [recip] = recipLattice(lattice, omega)
% recipLattice
% Calculate the reciprocal lattice vectors for the unit cell (supercell?)
%
% Inputs:
% lattice       3x3 array of lattice vectors. Rows are a, b, c,
%               columns are Cartesian x, y, z. Angstroms
% omega         Unit cell (supercell) volume, in Angstroms^3
%
% Outputs:
% recip         3x3 array of reciprocal lattice vectors. Rows are a, b, c,
%               columns are Cartesian x, y, z. Angstroms^-1

% Unpack lattice.
a = lattice(1,:);
b = lattice(2,:);
c = lattice(3,:);

% Calculate the reciprocal lattice vectors.
ar = 2*pi/omega*cross(b,c);
br = 2*pi/omega*cross(c,a);
cr = 2*pi/omega*cross(a,b);

% Pack up the reciprocal lattice vectors.
recip = [ar; br; cr];

end