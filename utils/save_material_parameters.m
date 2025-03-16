% ------------------------------------------------------------------------ 
% save_material_paramaters.m
%
% Description: Save the parameters for the material as a mat file
%
% Last modified: 02/09/2024, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
clear; 
close all; 
clc

% Define FEM parameters
elementType = 'TET4';
FORMULATION = 'N1t';  % N1/N1t/N0
VOLUME = 1;           % integration over defected (1) or nominal volume (0)
USEJULIA = 1;

% Material properties
E = 260000;   % Young's modulus [Pa]
rho = 1070;   % Density [kg/m^3]
nu = 0.4;     % Poisson's ratio 

% Create the material and element constructor
myMaterial = KirchoffMaterial();
set(myMaterial, 'YOUNGS_MODULUS', E, ...
                'DENSITY', rho, ...
                'POISSONS_RATIO', nu);
myMaterial.PLANE_STRESS = true;	    % set "false" for plane_strain
myElementConstructor = @() Tet4Element(myMaterial);

% Define which proportion of the fish is rigid
propRigid = 0.6;

% Save all variables to a .mat file
save('parameters.mat');
