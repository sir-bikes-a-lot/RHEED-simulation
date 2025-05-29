function [nuTDS] = calcnuTDS(E, B, Z, omega)
% calcnupl
% Calculate thermal diffuse scattering factor.
% Inputs:
% E             Incident electron beam energy, in eV
% B             Debye-Waller factor
% Z             Atomic number
% omega         Unit cell volume, in Angstroms^3
%
% Ouputs:
% nuTDS         Thermal diffuse scattering factor, in Ang^-3 eV^-1/2

% Check if Debye-Waller parameter is provided.
nuTDS = 0;
if (B ~= 0)
    nuTDS = log(15.4/(B*Z^(2/3)) + 1);
    nuTDS = nuTDS*4.74*B*(Z^2)/(omega*sqrt(E));     % Ang^-3 eV^-1/2
end

end