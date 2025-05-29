% Development script to calculate RHEED patterns.
addpath('Crystal models');
addpath('Libraries');
addpath('vasplab');

% Import parameters for each atom from library.
Ga = elements('Ga', 'GaN');
N = elements('N', 'GaN');
ElementsLib = struct('Ga', Ga, 'N', N);
% Zn = elements('Zn', 'ZnO');
% Ti = elements('Ti', 'ZnO');
% N = elements('N', 'ZnO');
% ElementsLib = struct('Zn', Zn, 'Ti', Ti, 'N', N);

% RHEED beam incident angle and azimuth.
theta = 3.0*pi/180;       % radians
psi = 60*pi/180;        % radians

% Specify temperature.
T = 250+273;            % Kelvin

% RHEED beam vacuum energy.
E0 = 20e3;          % eV

% Make some assumptions about inelastic scattering energy losses.
delEel = 70;        % single-electron excitation energy loss, in eV
delEpl = 10;        % surface plasmon scattering energy loss, in eV
lamf = 10;          % Electron wavelength at Fermi surface, in Angstroms
delrsq = 0.1^2;     % Mean atomic displacement squared, in Angstroms^2      

% Calculate electron beam wavevector.
K0 = calck0(E0, theta, psi);
lambda = 2*pi./K0.mag;      % Angstroms

% Import crystal model.
filename = 'bulk_GaN.vasp';
% filename = 'GaN_000-1_3x3_AOA.vasp';
% filename = 'ZnTiN wurtzite.vasp';
% filename = 'ZnTiN2-225_scan_disorderRS.vasp';

hkl = [3 3 1];
crystal = LoadCrystal(filename, hkl);

% Rotate the crystal to orient the (111) direction along the z-axis.
% A = [sqrt(2)/2, 0, -sqrt(2)/2;...
%      -sqrt(6)/6, sqrt(2)/sqrt(3), -sqrt(6)/6;...
%      sqrt(3)/3, sqrt(3)/3, sqrt(3)/3];

% crystal.r = transpose(A*transpose(crystal.r));

% Loop through all the atoms in the list.
M = length(crystal.atoms);
r = zeros(M, 3);

for m=1:M
    % Rotate crystal by azimuthal angle psi.
    r(m,:) = basistrans(transpose(crystal.r(m,:)), pi/2, psi);  % Angstroms
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the reciprocal lattice vectors.
% recip = recipLattice(A*crystal.lattice, crystal.omega);           % Angstroms^-1
recip = recipLattice(crystal.lattice, crystal.omega);           % Angstroms^-1
% recip = recipLattice(crystal.UClattice, crystal.UComega);           % Angstroms^-1

% Make a mesh of reciprocal lattice vectors in the plane of the sample 
% (a*,b*).
Na = 20;
Nb = 20;
[Bl, Sz, Bmn] = recipmesh(Na, Nb, recip, [0,0,0], K0.mag, K0.x, K0.y, K0.z);
x0ind = Na + 1;
y0ind = Nb + 1;

% Set up RHEED screen coordinate mesh from the reciprocal space mesh.
d = 30;             % cm
[xd, yd, thetaf, psif] = projscreen(d, K0.mag, Bl, Sz, theta, psi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the slices for the embedded R-matrix recursion.
Nz = 150;
z = linspace(min(r(:,3)), max(r(:,3)), Nz);     % Angstroms
h = z(2) - z(1);                                % Angstroms
   
% Initialize the Fourier components of the crystal potential.
U = zeros(2*Na+1, 2*Nb+1, Nz);
mu0 = zeros(1, Nz);

% DEBUG
% [U, mu0] = calcU("Ti", [0 0 0], crystal.lattice, [0 0], h, T, ElementsLib, E0, delEel, delEpl, lamf, delrsq, K0.mag);

% Loop on each reciprocal lattice vector.
for m=1:2*Na+1
    for n=1:2*Nb+1
        % Loop on each slice.
        for k=1:Nz
            % Calculate the Fourier components of the crystal potential for each
            % reciprocal lattice scattering vector and each slice.
            if norm(squeeze(Bl(m,n,:))) == 0
                [U(m,n,k), mu0(k)] = calcU(crystal.atoms, r, crystal.lattice, squeeze(Bl(m,n,:)), z(k), T, ElementsLib, E0, delEel, delEpl, lamf, delrsq, K0.mag);
            else
                [U(m,n,k), ~] = calcU(crystal.atoms, r, crystal.lattice, squeeze(Bl(m,n,:)), z(k), T, ElementsLib, E0, delEel, delEpl, lamf, delrsq, K0.mag);
            end
        end
    end
end

% Put all of the beams of the potential into a 1-D list.
sz = size(U);
Ul = []; Gamma = [];

for m=1:sz(1)
    for n=1:sz(2)
        Ul = [Ul, U(m,n,:)];
        Gamma = [Gamma, Sz(m,n) + K0.z];

        % Get the point corresponding to Bm = (0,0).
        if (m==x0ind)&&(n==y0ind)
            Gamma0 = Sz(m,n) + K0.z;
            m0 = length(Gamma);
        end
    end
end
Ul = squeeze(Ul);
Gamma = transpose(squeeze(Gamma));

% Perform the recursive embedded R-matrix calculations.
Rr = embeddedRmatrix(Ul, Gamma, Gamma0, z, m0);

% Pull the reflectivity ratios out of the 1-D list and put them back in a
% square matrix corresponding to the reciprocal space mesh.
I = zeros(2*Na+1, 2*Nb+1);
k = 1;
for m=1:2*Na+1
    for n=1:2*Nb+1
        if imag(Sz(m,n)) == 0
            I(m,n) = Rr(k);
        end
        k = k + 1;
    end
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate the simulated pattern to a finer (x, y) resolution and then
% apply 2D Gaussian broadening at each point.
% xq = min(min(xd)):0.05:max(max(xd));
% yq = min(min(real(yd))):0.05:max(max(real(yd)));
% [xq, yq] = meshgrid(xq, yq);
% 
% % Put all of the RHEED screen coordinates and intensities into column
% % vectors.
% xc = [];
% yc = [];
% Ic = [];
% 
% for m=1:2*Na+1
%     for n=1:2*Nb+1
%         if imag(yd(m,n))==0
%             xc = [xc, xd(m,n)];
%             yc = [yc, yd(m,n)];
%             Ic = [Ic, I(m,n)];
%         end
%     end
% end
% 
% % Set up interpolation.
% F = scatteredInterpolant(transpose(xc), transpose(yc), transpose(Ic));
% Iq = F(xq, yq);
% 
% % Set up a RHEED screen mesh for 2D Gaussian broadening.
% xlim = 1;                                   % cm
% ylim = 1;                                   % cm
% 
% Nl = 101;
% x = linspace(-xlim, xlim, Nl);              % cm
% y = linspace(-ylim, ylim, Nl);              % cm
% 
% % Specify the broadening of the RHEED screen points.
% sigx = 0.01;                                % cm
% sigy = 0.01;                                % cm
% 
% % Create the 2D Gaussian convolution function.
% H = zeros(Nl, Nl);
% 
% for i=1:Nl
%     for j=1:Nl
%         H(i,j) = exp(-((x(i))^2)/(2*sigx^2) -((y(j))^2)/(2*sigy^2));
%     end
% end
% 
% % Broaden the reciprocal lattice points using convolution (faster).
% L0 = conv2(Iq,H,'same');

% Apply a max threshold to the intensity map.
thresh = 5e-5;
% sz = size(L0);
% L = L0;
sz = size(I);
L = I;

for m=1:sz(1)
    for n=1:sz(2)
        % if (L0(m,n) > thresh)
        if (I(m,n) > thresh)
            L(m,n) = thresh;
        end
    end
end

% Plot RHEED pattern.
radius = 9.3/2;             % cm
figure;
contourf(xd, real(yd), L, 'LineStyle', 'none'); hold on;

% Plot the specular reflection.
scatter(0, d*tan(theta), 40, 'MarkerEdgeColor',[0 0.5 0.5],...
              'MarkerFaceColor',[0 0.7 0.7],...
              'LineWidth',1.5);

hold off;
xlabel('x_{d} (cm)'); ylabel('y_{d} (cm)');
title('RHEED pattern');
colormap gray;
% clim([0, max(max(max(L)), thresh)]);
% axis([-radius, radius, 0, 2*radius]);
axis([-3, 3, 0, 8]);
axis square
% savestr = strcat('GaN_000-1_3x3_AOA_1-100_20keV_3p0deg_hkl_221_Nz_150_Nab_20_5e-5');
% saveas(gcf, savestr, 'fig');
% print(gcf, strcat(savestr, '.png'),'-dpng','-r600'); 

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the magnitude of the potential vs z at a fixed reciprocal space
% vector.
i = x0ind; j = y0ind;

Ang = char(197);
FS = 10;

figure;
plot(z, abs(squeeze(U(i,j,:))), 'k-'); hold on;
% plot(z, real(squeeze(U(i,j,:))), 'b--');
% plot(z, imag(squeeze(U(i,j,:))), 'r:');
xlabel(strcat('z (', Ang, ')'), 'FontSize', FS);
ylabel(strcat('|Um| (', Ang, '^{-2})'), 'FontSize', FS);
grid on;
% savestr = strcat('GaN_000-1_3x3_AOA_20keV_hkl_331_Nz_150');
% saveas(gcf, savestr, 'fig');
% print(gcf, strcat(savestr, '.png'),'-dpng','-r600'); 