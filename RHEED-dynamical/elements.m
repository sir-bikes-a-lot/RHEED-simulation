function [ params ] = elements( element, material )
% ELEMENTS
% Get the parameters for a specified element. Return a structure of
% the following form:
%
% params = struct(  'Z', Z,...                      Atomic mass
%                   'DoyleTurner', DoyleTurner,...  A structure containing
%                                                   the Doyle-Turner
%                                                   coefficients.
%                   'u11Params', u11Params,...      A structure containing
%                                                   the static correlation
%                                                   function fit parameters.
%                   'u33Params', u33Params);        A structure containing
%                                                   the static correlation
%                                                   function fit parameters.
%
% The real part of the Doyle-Turner coefficients is taken from Table 3 in:
% L.M. Peng, "Electron atomic scattering factors and scattering potentials
% of crystals", Micron 30 (1999) 625–648
%
% The imaginary part of the Doyle-Turner coefficients is taken from Tables
% 3-4 in:
% S.L. Dudarev, L.-M. Peng, M.J. Whelan, "On the Doyle-Turner
% representation of the optical potential for RHEED calculations", Surface
% Science 330 (1995) 86-100
%
% The static correlation function fit parameters are taken from Tables 3
% and 4 in:
% M. Schowalter, A. Rosenauer, J. T. Titantah, and D. Lamoen, "Temperature-
% dependent Debye–Waller factors for semiconductors with the wurtzite-type 
% structure", Acta Cryst. (2009). A65, 227–231.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the sheet names.
sheets = sheetnames("elements.xlsx");

% Pull in the elemental scattering data from the spreadsheet.
T1 = readtable("elements.xlsx", 'Sheet', 'scattering');

% Pull in the atomic form factor and incoherent scattering function data
% from the spreadsheet.
if contains(sheets, element)
    T2 = readtable("elements.xlsx", 'Sheet', element);

    Incoherent = struct('S', T2.S,...
                    'FF', T2.FF,...
                    'Inc', T2.Inc);
else
    Incoherent = [];
end

% Correct the NaNs in the Doyle-Turner factor table.
[I,J] = size(T1);
for i=1:I
    for j=2:J
        if(isnan(T1{i,j}))
            T1{i,j} = 0;
        end
    end
end

T3 = readtable("elements.xlsx", 'Sheet', 'u11', 'Range', 'A1:E11');
T4 = readtable("elements.xlsx", 'Sheet', 'u33', 'Range', 'A1:E11');

% Now get just the element specified in the input argument.
ind = strcmp(T1.Element, element);
T1 = T1(ind,:);

ind = strcmp(T3.Material, material);
T3 = T3(ind,:);
ind = strcmp(T3.Element, element);
T3 = T3(ind,:);

ind = strcmp(T4.Material, material);
T4 = T4(ind,:);
ind = strcmp(T4.Element, element);
T4 = T4(ind,:);

% Build the element parameters structure.
Z = T1.Z;
M = T1.M;
DoyleTurner = struct('a1', T1.a1,...
                    'a2', T1.a2,...
                    'a3', T1.a3,...
                    'a4', T1.a4,...
                    'a5', T1.a5,...
                    'b1', T1.b1,...
                    'b2', T1.b2,...
                    'b3', T1.b3,...
                    'b4', T1.b4,...
                    'b5', T1.b5,...
                    'aTDS1', T1.aTDS1,...
                    'aTDS2', T1.aTDS2,...
                    'aTDS3', T1.aTDS3,...
                    'aTDS4', T1.aTDS4,...
                    'aTDS5', T1.aTDS5,...
                    'bTDS1', T1.bTDS1,...
                    'bTDS2', T1.bTDS2,...
                    'bTDS3', T1.bTDS3,...
                    'bTDS4', T1.bTDS4,...
                    'bTDS5', T1.bTDS5);
u11Params =  struct('sigma', T3.sigma,...
                    'A', T3.A,...
                    'B', T3.B);
u33Params =  struct('sigma', T4.sigma,...
                    'A', T4.A,...
                    'B', T4.B);
params = struct('Z', Z, 'M', M,...
                'DoyleTurner', DoyleTurner,...
                'Incoherent', Incoherent,...
                'u11Params', u11Params,...
                'u33Params', u33Params);

end

