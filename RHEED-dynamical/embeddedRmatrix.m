function [Rr] = embeddedRmatrix(Ul, Gamma, Gamma0, z, m0)
% embeddedRmatrix
% Perform the embedded R-matrix calculations as outlined in:
% A. Ichimiya and P. I. Cohen, "Dynamical theory - embedded R-matrix
% method", Chapter 13 in "Reflection High-Energy Electron Diffraction",
% Cambridge University Press (2004).
%
% Inputs:
% Ul            Array of Fourier components of the crystal potential, in
%               Angstroms^-2
% Gamma         Array of scattered wavevector z-components, in Angstroms^-1
% Gamma0        Scattered wavevector z-component of zeroth beam, in Ang^-1
% z             Vector of z-values for the slices, in Angstroms
% m0            Index of zeroth beam
%
% Outputs:
% Rr            Final reflectivity matrix array elements (RHEED intensity)

% Initialize the transfer matrix calculation at the bottom of the slab.
sz = size(Ul);
M = sz(1);
b = max(z);
h = z(2)-z(1);
Nz = length(z);
Qlast = ones(M,M);
Rlast = -1i*diag(Gamma.^2);

% Loop through all the z slices.
for j=1:Nz
    % Construct the transfer matrix.
    A = constructA(Gamma, Ul(:,j));
    
    % Diagonalize the matrix A to find eigenvector matrix Q and eigenvalues
    % Lam.
    [Q, Lam] = eig(A);
    gammasq = diag(Lam);                % Angstroms^-2
    gamma = sqrt(gammasq);              % Angstroms^-1
    
    % Calculate the matrix T.
    T = Q'*Qlast;
    
    % Calculate the matrices R1-R4.
    R1 = diag(-gamma.*cot(gamma*h));        % dimensionless
    R2 = diag(gamma.*csc(gamma*h));         % dimensionless
    R3 = diag(-gamma.*csc(gamma*h));        % dimensionless
    R4 = diag(gamma.*cot(gamma*h));         % dimensionless
    
    % Calculate the R-matrix.
    R = R4 + R3*(T*Rlast*T' - R1)\R2;       % dimensionless

    % Assign the Q- and R-matrix values to the last slice before proceeding
    % to the next slice.
    Qlast = Q;
    Rlast = R;
end

% Construct the matrices Imn and Gmn.
Imn = diag(exp(-1i*Gamma*b));       % dimensionless
Gmn = diag(Gamma);                  % Angstroms^-1

% Calculate the final R-matrix.
Rf = Q*R*Q';

% Calculate the reflectivity matrix M.
Mm = Imn*(1i*Gmn - Rf)\(1i*Gmn + Rf)*Imn;

% Calculate reflectivity ratio.
Rr = Mm.^2;
Rr = abs(Rr(:,m0));
Rr = Gamma./Gamma0.*Rr;

end