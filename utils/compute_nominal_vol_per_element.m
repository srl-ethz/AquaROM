% compute_nominal_vol_per_element.m
%
% Synthax:
% vol = compute_nominal_muscle_vol_per_element(Mesh, muscleElements)
%
% Description: Computes the volume of each muscle element
%
% INPUTS: 
% (1) Mesh:             FE mesh object
% (2) muscleElements:   array of size (nElements,1) with 1 if an element is
%                       part of the mucle(s), 0 otherwise 
%
% OUTPUTS:   
% (1) volVector:        a vector of element volumes of the size of the
%                       total number of elements.
%
% Last modified: 05/02/2026, Mathieu Dubied, ETH Zürich
function volVector = compute_nominal_vol_per_element(Mesh,nel)
    volVector = zeros(nel,1);
    for el=1:nel
        volVector(el) = Mesh.Elements(el).Object.vol;
    end
end