function [A] = constructA(Gamma, Um)
% constructA
% Construct the transfer matrix A from the out-of-plane scattered
% wavevector Gamma and the Fourier components of the crystal potential, Um.
%
% Inputs:
% Gamma         Vector of out-of-plane scattered wavevectors, in Angstroms^-1
% Um            Vector of Fourier components of the crystal potential, in
%               Angstroms^-2
%
% Outputs:
% A             Transfer matrix, in Angstroms^-2

% Initialize matrix A.
M = length(Um);
N = (M-1)/2;
A = zeros(M, M);

% Loop through the matrix elements and assign them.
for m=-N:N
    for n=-N:N
        if m==n
            A(m+N+1,n+N+1) = Gamma(m+N+1)^2 + Um(m-n+N+1);
        elseif ((m-n)>=-N)&&((m-n)<=N)
            A(m+N+1,n+N+1) = Um(m-n+N+1);
        else
            A(m+N+1,n+N+1) = 0;
        end
    end
end

end