% ------------------------------------------------------------------------ 
% Performance comparison for different ROB choices.
% 
% Last modified: 18/03/2024, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc

elementType = 'TET4';

FORMULATION = 'N1t'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

USEJULIA = 1;

%% PREPARE MODEL __________________________________________________________                                                   

% DATA ____________________________________________________________________
E       = 260000;      % Young's modulus [Pa]
rho     = 1070;         % density [kg/m^3]
nu      = 0.4;        % Poisson's ratio 

% material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	    % set "false" for plane_strain
myElementConstructor = @()Tet4Element(myMaterial);

% MESH ____________________________________________________________________

filename =  '3d_fish_for_mike';%'3d_rectangle_660el' ;%'3d_fish_for_mike';
[nodes, elements, ~, elset] = mesh_ABAQUSread(filename);
% nodes = 0.01*nodes;
% nodes(:,2)=0.8*nodes(:,2);
nel = size(elements,1);

MeshNominal = Mesh(nodes);
MeshNominal.create_elements_table(elements,myElementConstructor);

Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil
Lz = abs(max(nodes(:,3))-min(nodes(:,3)));  % vertical length of airfoil

% plot nominal mesh
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMeshAxis(nodes, elementPlot, 0);
hold off

% boundary conditions of nominal mesh


% 
% for l=1:length(nset)
%     MeshNominal.set_essential_boundary_condition([nset{l}],3,0)
% end   


% desiredFixedPoint = - 0.45*Lx;
% fixedPoint = find_fixed_point(nodes,desiredFixedPoint);
% nel = size(elements,1);
% nset = {};
% marginFixedPoint = 0.02;
% 
% nset = {};
% for el=1:nel   
%     for n=1:size(elements,2)
%         if  ~any(cat(2, nset{:}) == elements(el,n))
%             nset{end+1} = elements(el,n);
%         end
%     end   
% end
% 
% for l=1:length(nset)
%     MeshNominal.set_essential_boundary_condition([nset{l}],3,0)
% end 
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

% for testing
U = [z_tail,y_head];

%% SIMULATION PARAMETERS __________________________________________________
h = 0.01;
tmax = 2;

%% ROB ____________________________________________________________________

% build ROM
nel = size(elements,1);

fixedPortion = 0.6;
nset1 = {};
fixedElements = zeros(nel,1);
for el=1:nel   
    for n=1:size(elements,2)
        if  nodes(elements(el,n),1)>-Lx*fixedPortion && ~any(cat(2, nset1{:}) == elements(el,n))
            nset1{end+1} = elements(el,n);    
        end
        
        if  nodes(elements(el,n),1)>-Lx*fixedPortion
            fixedElements(el)=1;
        end
    end   
end

% desiredFixedPoint = - 0.45*Lx;
% fixedPoint = find_fixed_point(nodes,desiredFixedPoint);
% nel = size(elements,1);
% nset = {};
% marginFixedPoint = 0.02;
% 
% nset1 = {};
% for el=1:nel   
%     for n=1:size(elements,2)
%         if  nodes(elements(el,n),1)> fixedPoint-marginFixedPoint && ...
%                 nodes(elements(el,n),1)< fixedPoint+marginFixedPoint && ...
%                 ~any(cat(2, nset1{:}) == elements(el,n))
%             nset1{end+1} = elements(el,n);
%         end
%     end   
% end

% for l=1:length(nset)
%     MeshNominal.set_essential_boundary_condition([nset{l}],3,0)
% end 

nset2 = {};
% for el=1:nel   
%     for n=1:size(elements,2)
%         if  ~any(cat(2, nset2{:}) == elements(el,n))
%             nset2{end+1} = elements(el,n);
%         end
%     end   
% end
translation = 1;    %0,1,or 2


figure
set(gcf, 'Position',  [100, 100, 1200, 500])
PlotFieldonDeformedMeshActuation2Muscles(nodes,elements,fixedElements,1,fixedElements,1,zeros(length(nodes),3));

%%
tStartROM = tic;
fprintf('____________________\n')
fprintf('Building ROM ... \n')
[V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
build_ROM_3D_extended(nodes,elements,nset1,nset2,translation,myElementConstructor,USEJULIA);  

% solve EoMs 
tic 
fprintf('____________________\n')
fprintf('Solving EoMs ...\n') 
TI_NL_ROM = solve_EoMs(V,ROM_Assembly,tensors_ROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax); 
toc


fprintf('Time needed to solve the problem using ROM: %.2fsec\n',toc(tStartROM))
timeROM = toc(tStartROM);




%% PLOT ___________________________________________________________________

% displacement of the tail node
uTail = zeros(3,tmax/h);
timePlot = linspace(0,tmax-h,tmax/h);
x0Tail = min(nodes(:,1));
for a=1:tmax/h
    uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_ROM.Solution.q(:,a);
end

f1 = figure;
subplot(2,1,1);
plot(timePlot,x0Tail+uTail(1,:),'DisplayName','k=0')
hold on
grid on
xlabel('Time [s]')
ylabel('x-position tail node')
legend('Location','northwest')
subplot(2,1,2);
plot(timePlot,uTail(2,:),'DisplayName','k=0')
hold on
grid on
xlabel('Time [s]')
ylabel('y-position tail node')
legend('Location','southwest')
sgtitle('Tail node displacement')

documentTitle = ['tail_node_',num2str(fixedPortion,'%.2f '),'_translation_',num2str(translation,'%.0f '),'.svg'];
print(f1,documentTitle,'-dsvg','-r800');

%%
% envelope of the oscillation 
[spineNodes, ~, ~, ~] = find_spine_TET4(elements,nodes);

f2 = figure;
hold on
for t = 1:length(TI_NL_ROM.Solution.q(1,:))
    xPlotUnsorted = nodes(spineNodes,1);
    [xPlot,I] = sort(xPlotUnsorted);
    yPlotUnsorted = V(spineNodes*3-1,:)*TI_NL_ROM.Solution.q(:,t);
    yPlot = yPlotUnsorted(I);
    plot(xPlot,yPlot);
end
grid on
axis equal
xlabel('x-position along fish spine (head at 0)')
ylabel('y-position of oscilation')
title('Kinematics of the fish spine')
xlim([-0.35,0])
ylim([-0.05,0.05])

documentTitle = ['spine_kinematics_',num2str(fixedPortion,'%.2f '),'_translation_',num2str(translation,'%.0f '),'.svg'];
print(f2,documentTitle,'-dsvg','-r800');


%% ANIMATION ______________________________________________________________
elementPlot = elements(:,1:4); 
nel = size(elements,1);

% top muscle
topMuscle = zeros(nel,1);

for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY>0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
        topMuscle(el) = 1;
    end    
end

% bottom muscle
bottomMuscle = zeros(nel,1);
for el=1:nel
    elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
    elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
    if elementCenterY<0.00 &&  elementCenterX < -Lx*0.25 && elementCenterX > -Lx*0.8
        bottomMuscle(el) = 1;
    end    

end

actuationValues = zeros(size(TI_NL_ROM.Solution.u,2),1);
for t=1:size(TI_NL_ROM.Solution.u,2)
    actuationValues(t) = 0;
end

actuationValues2 = zeros(size(TI_NL_ROM.Solution.u,2),1);
for t=1:size(TI_NL_ROM.Solution.u,2)
    actuationValues2(t) = 0;
end
sol = TI_NL_ROM.Solution.u(:,1:2:end);

% documentTitle = ['results_video_',num2str(fixedPortion,'%.2f '),'_translation_',num2str(translation,'%.0f ')];
documentTitle = 'results_video';
AnimateFieldonDeformedMeshActuation2Muscles(nodes, elementPlot,topMuscle,actuationValues,...
    bottomMuscle,actuationValues2,sol, ...
    'factor',1,'index',1:3,'filename',documentTitle,'framerate',1/h*0.5)






