function [xd, yd, thetaf, psif] = projscreen(d, K0, Bm, Sz, theta, psi)
% projscreen
% Project reciprocal lattice scattering vectors onto the RHEED screen.
%
% Inputs:
% d             Distance from sample to RHEED screen, in cm.
% K0            Scattered wavevector amplitude, assumed equal to the
%               incident wavevector for elastic scattering.
% Bm            Array of reciprocal lattice vectors parallel to sample
%               surface. In Cartesian (x,y) coordinates.
% Sz            Vector of reciprocal lattice vectors perpendicular to 
%               sample surface. In Cartesian (z) coordinates.
% theta         RHEED beam incident angle, in radians.
% psi           RHEED beam incident azimuth, in radians.
%
% Outputs:
% xd            Array of screen positions in x-direction, in cm
% yd            Array of screen positions in y-direction, in cm
% thetaf        Final scattered wavevector takeoff angle, in radians
% psif          Final scattered wavevector azimuth, in radians

% Calculate final scattered wavevector takeoff angle.
thetaf = asin(Sz/K0 - sin(theta));                      % radians

% Get the reciprocal lattice vector x and y components.
Bx = squeeze(Bm(:,:,1));
By = squeeze(Bm(:,:,2));

% Calculate final scattered wavevector azimuth. Offset by the incident beam
% azimuth.
psif = atan((By/K0 + sin(psi))./(Bx/K0 + cos(psi)));    % radians
psif = psif - psi;                                      % radians

% Calculate the RHEED screen positions from the final scattered wavevector
% angles.
xd = d*tan(psif);           % cm
yd = d*tan(thetaf);         % cm

end