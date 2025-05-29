function [nuel] = calcnuel(E, delEel, Z, omega)
% calcnuel
% Calculate single-electron excitation scattering factor.
% Inputs:
% E             Incident electron beam energy, in eV
% delEel        Single-electron excitation energy loss, in eV
% Z             Atomic number
% omega         Unit cell volume, in Angstroms^3
%
% Ouputs:
% nuel          Single-electron excitation scattering factor, in V

nuel = log((8.82*sqrt(E)*Z^(1/3))/delEel) - 1/4;
nuel = nuel*147*(Z^(1/3))/(omega*sqrt(E));          % Ang^-3 eV^-1/2

end