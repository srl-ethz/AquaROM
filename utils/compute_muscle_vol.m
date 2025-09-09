% compute_muscle_vol
%
% Synthax:
% vol = compute_muscle_vol(Mesh, muscleElements)
%
% Description: Computes the total volume of the muscle(s)
%
% INPUTS: 
% (1) Mesh:             FE mesh object
% (2) muscleElements:   array of size (nElements,1) with 1 if an element is
%                       part of the mucle(s), 0 otherwise 
%
% OUTPUTS:   
% (1) vol:              the total volume of the muscle, in cm^3
%
% Last modified: 17/05/2025, Mathieu Dubied, ETH Zürich
function vol = compute_muscle_vol(Mesh, muscleElements)
    vol = 0;
    nel = size(muscleElements,1);
    vol = 0;
    for el=1:nel
        if muscleElements(el) == 1
            vol = vol + Mesh.Elements(el).Object.vol * 1e6; %cm^3
        end    
    end
end