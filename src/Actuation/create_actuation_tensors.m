% create_actuation_tensors
%
% Synthax:
% Actu = create_actuation_tensors(Assembly, elements, nodes)
%
% Description: This function returns the vectors and matrices related to
% the actuation muscles in the FOM
%
% INPUTS: 
% (1) Assembly:         Assembly from YAFEC               
% (2) elements:         matrix that describes which nodes (columns) in each
%                       element (rows)
% (3) nodes:            position x y of the nodes
%
% OUTPUTS:
% (1) TActu:            struct containing the B1 vectors and B2 matrices
%                       for the top and bottom muscles
%     
%
% Additional notes: -
%
% Last modified: 10/10/2023, Mathieu Dubied, ETH Zürich
function TActu = create_actuation_tensors(Assembly, elements, nodes)
    
    % get size of fish and important features _____________________________
    Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of the nominal fish
    Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of the nominal fish

    nel = size(elements,1);
    actuationDirection = [-1;0;0];               %[1;0]-->[1;0;0] (Voigt notation)

    disp(' FOM ACTUATION TENSORS:')
    fprintf(' Assembling %d elements ...\n', nel)
    
    % top muscle __________________________________________________________
    topMuscle = zeros(nel,1);
    for el=1:nel
        elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
        elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
        if elementCenterY>0.00 &&  elementCenterX <-Lx*0.25
            topMuscle(el) = 1;
        end    
    end

    tic
    B1_top = Assembly.vector('B1', 'weights', topMuscle, actuationDirection);
    fprintf('   B1 top: %.2f s\n',toc)

    tic;
    B2_top = Assembly.matrix_actuation('B2', 'weights', topMuscle, actuationDirection);
    fprintf('   B2 top: %.2f s\n',toc)

    
    % bottom muscle _______________________________________________________
    bottomMuscle = zeros(nel,1);
    for el=1:nel
        elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2))/3;
        elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1))/3;
        if elementCenterY<0.00 &&  elementCenterX <-Lx*0.25
            bottomMuscle(el) = 1;
        end    
    end
    
    tic
    B1_bottom = Assembly.vector('B1', 'weights', topMuscle, actuationDirection);
    fprintf('   B1 bottom: %.2f s\n',toc)

    tic;
    B2_bottom = Assembly.matrix_actuation('B2', 'weights', topMuscle, actuationDirection);
    fprintf('   B2 bottom: %.2f s\n',toc)

    TActu.B1_top = B1_top;
    TActu.B2_top = B2_top;
    TActu.B1_bottom = B1_bottom;
    TActu.B2_bottom = B2_bottom;
end