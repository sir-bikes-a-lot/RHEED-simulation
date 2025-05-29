function [ B ] = DebyeWaller( u11Params, u33Params, M, T )
% DebyeWaller
% Calculate Debye-Waller factor B from the fitting parameters for wurtzite
% compound semiconductors. Based on M. Schowalter, A. Rosenauer,
% J. T. Titantah, and D. Lamoen, "Temperature-dependent Debye–Waller 
% factors for semiconductors with the wurtzite-type structure", Acta Cryst.
% (2009). A65, 227–231.
%
% Inputs:
% u11Params         Fitting parameters for atomic displacements in the a=b
%                   direction.
% u33Params         Fitting parameters for atomic displacements in the c
%                   direction.
% M                 Atomic mass, in Daltons
% T                 Temperature, in Kelvin
%
% Outputs:
% B                 Debye-Waller factor, in Angstroms^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if Debye-Waller parameters are provided. if not, then just return
% 0.
B = 0;
if ~isempty(u11Params.A)
    % Define some physical constants.
    hbar = 1.05457e-34;             % J-s = kg*m^2*s^-1
    kB = 1.380649e-23;              % J/K = kg*m^2*s^-1*K^-1
    Da = 1.660539e-27;              % kg
    
    % Unpack the fitting parameters.
    sigma11 = u11Params.sigma;
    A11 = u11Params.A;
    B11 = u11Params.B;
    
    sigma33 = u33Params.sigma;
    A33 = u33Params.A;
    B33 = u33Params.B;
    
    % Calculate the static correlation functions from the fitting parameters.
    num = coth(hbar/(2*kB*T)*(A11*exp(-(T^2)/sigma11^2) + B11));    % dimensionless
    den = A11*exp(-(T^2)/sigma11) + B11;                            % Hz
    u11 = hbar/(2*M*Da)*num/den*(1e20);                             % Angstroms^2
    
    num = coth(hbar/(2*kB*T)*(A33*exp(-(T^2)/sigma33^2) + B33));    % dimensionless
    den = A33*exp(-(T^2)/sigma33) + B33;                            % Hz
    u33 = hbar/(2*M*Da)*num/den*(1e20);                             % Angstroms^2
    
    % Calculate the isotropic static correlation function and convert to B.
    usq = 2*u11 + u33;              % Angstroms^2
    B = 8*(pi^2)*usq;               % Angstroms^2
end
end