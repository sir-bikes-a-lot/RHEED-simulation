function [ crystal ] = LoadCrystal(filename, hkl)
% LoadCrystal
% Load a VASP POSCAR file containing the crystal information, including the
% atom IDs and coordinates. Uses the VASPLAB package: https://github.com/max-radin/VASPLAB
%
% Inputs:
% filename      String containing the .vasp file name.
% hkl           Triplet with supercell size, integers >= 1.
%
% Outputs:
% crystal       Structure containing the atom names and Cartesian
%               coordinates (x,y,z), and the lattice.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the VASPLAB filepath.
addpath('vasplab');

% Load the .vasp file.
geometry = import_poscar(filename);

% Generate a supercell of size (h,k,l).
geometryExt = supercell(geometry, hkl);

% Unpack the lattice basis vectors, atom names and counts, and fractional
% coordinates.
lattice = geometryExt.lattice;          % Lattice basis vectors in Cartesian x,y,z coordinates, in Angstroms
atomnames = geometryExt.symbols;

% Generate a list of the atom names.
N = length(geometryExt.coords);
atoms = strings(N,1);
ind = 0;

for m=1:length(geometryExt.atomcount)
    for n=1:geometryExt.atomcount(m)
        atoms(n + ind) = cell2mat(atomnames(m));
    end
    ind = geometryExt.atomcount(m) + ind;
end

% Calculate the Cartesian coordinates of each atom in the supercell.
r = geometryExt.coords*geometryExt.lattice;

% Offset all the z-coordinates such that z=0 corresponds to the lowest
% fractional coordinate.
z0 = min(r(:,3));
r(:,3) = r(:,3) - z0;

% Calculate the cell volumes.
UComega = dot(geometry.lattice(1,:), cross(geometry.lattice(2,:), geometry.lattice(3,:)));
omega = dot(lattice(1,:), cross(lattice(2,:), lattice(3,:)));

% Pack up the atom names and coordinates to the output structure.
crystal = struct('atoms', atoms, 'r', r, 'lattice', lattice, 'UClattice', geometry.lattice, 'omega', omega, 'UComega', UComega);

end