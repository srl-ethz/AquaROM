% get_muscle_elements_and_vol
%
% Synthax:
% [muscleElements, vol] = get_muscle_elements_and_vol(Mesh, nodes, ...
%                                elements, muscleBoundaries)
%
% Description: Get the muscle elements and computes the total volume of 
%              the muscle(s)
%
% INPUTS: 
% (1) Mesh:             FE mesh object
% (2) nodes:            nodes and their coordinates
% (3) elements:         elements described by their nodes
%
% OUTPUTS:   
% (1) muscleElements:   array of size (nElements,1) with 1 if an element is
%                       part of the mucle(s), 0 otherwise 
% (2) vol:              the total volume of the muscle, in cm^3
%
% Last modified: 17/05/2025, Mathieu Dubied, ETH Zürich
function [muscleElements, vol] = get_muscle_elements_and_vol(Mesh, nodes, ...
                                elements, muscleBoundaries)
    [Lx,~,~] = mesh_dimensions_3D(nodes);
    nel = size(elements,1);
    muscleElements = zeros(nel,1);
    vol = 0;
    for el=1:nel
        elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
        if  elementCenterX<-Lx*muscleBoundaries(2) && elementCenterX>-Lx*muscleBoundaries(1)
            muscleElements(el) = 1;
            vol = vol + Mesh.Elements(el).Object.vol* 1e6; %cm^3
        end    
    end
end