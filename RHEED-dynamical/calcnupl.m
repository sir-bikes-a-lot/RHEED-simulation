function [nupl] = calcnupl(E, delEpl, lamf)
% calcnupl
% Calculate surface plasmon scattering factor.
% Inputs:
% E             Incident electron beam energy, in eV
% delEpl        Surface plasmon energy loss, in eV
% lamf          Electron wavelength at the Fermi surface, in Angstroms
%
% Ouputs:
% nupl          Surface plasmon scattering factor, in eV^1/2

nupl = log(lamf*sqrt(E)/12.24);
nupl = nupl*1.96*delEpl/sqrt(E);        % eV^1/2

end