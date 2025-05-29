function [Kout] = calck0(E, theta, psi)
% calck0
% Calculate electron beam wavevector k0.
% Inputs:
% E         Electron beam energy, in eV
% theta     Incident beam angle, in radians
% psi       Incident beam azimuth, in radians
%
% Ouputs:
% Kout      Structure with the magnitude and (x,y,z) components, in Ang^-1

hbar = 1.0546e-34; % J*s
m0 = 9.1094e-31;    % kg
c = 3e8;            % m/s
eV = 1.6022e-19;    % J = kg*m^2/s^2
K0 = 1/hbar*sqrt(2*m0*E*eV + ((E*eV).^2)/(c^2))*1e-10;       % Ang^-1

% Calculate the (x,y,z) components.
K0x = K0*cos(theta)*cos(psi);       % Angstroms^-1
K0y = K0*cos(theta)*sin(psi);       % Angstroms^-1
K0z = -K0*sin(theta);               % Angstroms^-1

% Return as a structure.
Kout = struct('mag', K0, 'x', K0x, 'y', K0y, 'z', K0z);

end