function [xd, yd, S, x0ind, y0ind] = calcSmesh(r, d, thetai, phii, k0, Nx, Ny)
% calcSmesh
% Setup a mesh in reciprocal space sitting on the Ewald sphere with radius
% k0. Do this by working backwards from the RHEED screen size r.
%
% Inputs:
% r         RHEED screen radius, in cm
% d         Distance from sample to RHEED screen
% thetai    Incident RHEED beam angle, in radians
% phii      Incident RHEED beam azimuth with respect to x, in radians
% k0        Incident electron beam wavevector, in Angstroms^-1
% Nx        Number of mesh points in x direction
% Ny        Number of mesh points in y direction
%
% Outputs:
% xd, yd    2D arrays of RHEED screen coordinates, in cm       
% S         2D mesh of reciprocal space coordinates (Sx, Sy, Sz), in
%           Angstroms^-1
% x0ind     x-index corresponding to specular reflection (Sx,Sy) = (0,0)
% y0ind     y-index corresponding to specular reflection (Sx,Sy) = (0,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin by setting up the RHEED screen coordinates. Make a symmetric x-mesh
% with a point at (0,d*tan(thetai)). Use an odd number of points.
if(mod(Nx, 2)==0)
    Nx = Nx + 1;
end

if(mod(Ny, 2)==0)
    Ny = Ny + 1;
end

% Square mesh
if Nx==1
    x = 0;
    x0ind = 1;
else
    x = linspace(-r, r, Nx);
    x0ind = (Nx-1)/2 + 1;
end
if Ny==1
    y = d*tan(thetai);
    y0ind = 1;
else
    ystep = 2*r/Ny;
    y1 = linspace(0, d*tan(thetai), d*tan(thetai)/ystep+1);
    y2 =linspace(d*tan(thetai), 2*r, (2*r-d*tan(thetai))/ystep+1);
    y = [y1, y2(2:end)];
    y0ind = length(y1);
    % y = linspace(d*tan(thetai), 2*r, Ny);
end
[xd, yd] = meshgrid(x,y);
xd = transpose(xd); yd = transpose(yd);

% Circular mesh
% rd = linspace(0, r, Nx);
% psi = linspace(-pi, pi, Ny);
% xd = zeros(Nx, Ny); yd = xd;
% for i=1:Nx
%     for j=1:Ny
%         xd(i,j) = rd(i)*cos(psi(j));
%         yd(i,j) = rd(i)*sin(psi(j));
%     end
% end

% Calculate the origin of the Ewald sphere, adopting the convention that
% the tip of the incident wavevector ends at (0,0,0).
Sx0 = -k0*cos(thetai)*cos(phii);
Sy0 = -k0*cos(thetai)*sin(phii);
Sz0 = k0*sin(thetai);

% Loop on each mesh point.
S = zeros(Nx, Ny, 3);

for n=1:Ny
    % Calculate Sz.
    % S(:,n,3) = k0*(sin(thetai) + yd(:,n)./sqrt(yd(:,n).^2 + d^2));      % Angstroms^-1
    S(:,n,3) = k0*(sin(thetai) + yd(:,n)./sqrt(yd(:,n).^2 + d^2));      % Angstroms^-1

    for m=1:Nx
        % Solve quadratic equation for Sx.
        a = 1 + (xd(m,n)/d)^2;
        b = 2*(k0*cos(phii)*(xd(m,n)/d)^2 - Sx0 - (xd(m,n)/d)*(k0*sin(phii) + Sy0));
        c = Sx0^2 + (xd(m,n)/d)^2*k0^2*(cos(phii))^2 - 2*(xd(m,n)/d)*k0*cos(phii)*(k0*sin(phii) + Sy0) + (k0*sin(phii) + Sy0)^2 + (S(m,n,3) - Sz0)^2 - k0^2;

        % Take positive solution branch only.
        S(m,n,1) = (-b + sqrt(b^2 - 4*a*c))/(2*a);                          % Angstroms^-1

        % Solve for Sy.
        S(m,n,2) = (xd(m,n)/d)*(S(m,n,1) + k0*cos(phii)) - k0*sin(phii);    % Angstroms^-1
    end
end
end