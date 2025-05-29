function [Ul, mu0] = calcU(atoms, r, lattice, Bl, z, T, elements, E, delEel, delEpl, lamf, delrsq, k0)
% calcU
% Calculate crystal potential at reciprocal lattice vector Bm and height z.
%
% Inputs:
% atoms         String containing the atom names.
% r             Cartesian coordinates of each atom in the supercell, in Angstroms
% lattice       Supercell basis vectors in Cartesian x,y,z coordinates, in Angstroms
% Bl            Reciprocal lattice vector, in Angstroms^-1
% z             Height in multi-slice slab, in Angstroms
% T             Temperature, in K
% elements      Library of elemental parameters
% E             Incident electron beam energy, in eV
% delEel        Single-electron excitation energy loss, in eV
% delEpl        Plasmon scattering energy loss, in eV
% lamf          Electron wavelength at the Fermi surface, in Angstroms
% delrsq        Mean atomic displacement squared, in Angstroms^2
% k0            Magnitude of incident wavevector, in Angstroms^-1
%
% Outputs:
% Ul            Complex crystal potential, in Angstroms^-2
% mu0           Mean absorption coefficient, in Angstroms^-1. Only
%               calculated for |Bl| = 0.

% Get the number of atoms in the supercell.
N = length(atoms);

% Calculate the magnitude of the reciprocal lattice vector Bm.
Bmag = norm(Bl);            % Angstroms^-1
% DEBUG: Turn off Bl-dependence of crystal potential.
% Bmag = 1;

% Calculate area and volume of supercell, A0 and omega.
A0 = norm(cross(lattice(1,:), lattice(2,:)));                   % Angstroms^2
omega = dot(lattice(1,:), cross(lattice(2,:), lattice(3,:)));   % Angstroms^3

% Get the magnitude of the first reciprocal lattice vector in the plane of
% the sample. Approximate this as the geometric mean of the (a,b) vectors.
a = norm(lattice(1,:)); b = norm(lattice(2,:));
sA = sqrt(a*b);                                 % Angstroms^-1

% Initialize auxiliary variable Uk.
Uk = zeros(1,4);

% Calculate the real part of the crystal potential. Loop on the 4 
% Doyle-Turner parameters.
for k=1:4
    % Initialize auxiliary variable Us.
    Us = zeros(1,N);

    % Loop on all the atoms in the supercell.
    for s=1:N
        % Fetch the Doyle-Turner coefficients from the database of 
        % elements.
        element = elements.(atoms(s));
        DT = element.DoyleTurner;
        u11Params = element.u11Params;
        u33Params = element.u33Params;

        % Calculate the Debye-Waller factor B. Neglect the anisotropy and calculate
        % the isotropic factor.
        B = DebyeWaller(u11Params, u33Params, element.M, T);

        % Unpack the Doyle-Turner parameters.
        aRe = [DT.a1, DT.a2, DT.a3, DT.a4, DT.a5];                 % Angstroms
        bRe = [DT.b1, DT.b2, DT.b3, DT.b4, DT.b5];                 % Angstroms^2
        % aIm = [DT.aTDS1, DT.aTDS2, DT.aTDS3, DT.aTDS4, DT.aTDS5];  % Angstroms
        % bIm = [DT.bTDS1, DT.bTDS2, DT.bTDS3, DT.bTDS4, DT.bTDS5];  % Angstroms^2

        % Calculate the k-th part of the potential due to atom s.
        Us(s) = aRe(k)*sqrt(pi/(bRe(k)+B))*exp(-(bRe(k)+B)*Bmag^2/(16*pi^2))*exp(-1i*dot(Bl, r(s,1:2)))*exp((-4*pi^2*(z - r(s,3))^2)/(bRe(k)+B));
    end

    % Sum over all atoms s.
    Uk(k) = sum(Us);
end

% Sum over all Doyle-Turner parameters k.
Ul = (8*pi/A0)*sum(Uk);         % Angstroms^-2

% Calculate the imaginary part of the crystal potential due to inelastic
% scattering. Loop on all the atoms in the supercell.
% Initialize auxiliary variable Us.
Us = zeros(1,N);
nusum = 0;
s0sq = 1/delrsq;                % Angstroms^-2

for s=1:N
    % Get elemental parameters.
    element = elements.(atoms(s));
    u11Params = element.u11Params;
    u33Params = element.u33Params;

    % Calculate the Debye-Waller factor B. Neglect the anisotropy and calculate
    % the isotropic factor.
    B = DebyeWaller(u11Params, u33Params, element.M, T);

    % Calculate the single-electron inelastic scattering factor.
    nuel = calcnuel(E, delEel, element.Z, omega);

    % Calculate the plasmon scattering factor.
    if Bmag == 0
        nupl = calcnupl(E, delEpl, lamf);
    else
        nupl = 0;
    end

    % Calculate the thermal diffuse scattering factor.
    nuTDS = calcnuTDS(E, B, element.Z, omega);

    % Calculate the imaginary part of the crystal potential.
    term1 = pi*nuel*(sA^2)/(sA^2 + Bmag^2)*exp(-(sA^2 + Bmag^2)*(z - r(s,3))^2);
    term2 = nupl;
    term3 = nuTDS*sqrt(pi*s0sq)*exp(-((Bmag^2)/s0sq + (s0sq*(z - r(s,3))^2)/4));
    Us(s) = exp(-1i*dot(Bl, r(s,1:2)))*(term1 + term2 + term3);

    % Calculate the sum of the scattering factors, for calculation of the
    % mean absorption coefficient.
    nusum = nusum + nuel + nuTDS;
end

% Sum over all atoms s. Add to real part of the crystal potential.
Ul = Ul + 1i*2/(A0*N)*sum(Us);              % Angstroms^-2

% Sum the inelastic scattering factors over all atoms s and compute the
% mean absorption coefficient.
m = 9.10938e-31;                % kg
q = 1.60218e-16;                % C
hbar = 1.05457e-34;             % J-s

if Bmag == 0
    mu0 = 2*m*q/(10^20*k0*hbar^2)*(nusum + nupl)/N;     % Angstroms^-1
else
    mu0 = [];
end
end