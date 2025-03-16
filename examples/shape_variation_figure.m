% ------------------------------------------------------------------------ 
% Creation of a figure containing all shape variation used
% 
% Last modified: 27/09/2024, Mathieu Dubied, ETH Zurich
%
% ------------------------------------------------------------------------
clear; 
close all; 
clc
set(groot,'defaulttextinterpreter','latex');

%% PREPARE NOMINAL MESH ___________________________________________________                                                   

% load material parameters
load('parameters.mat') 

% specify and create FE mesh
filename ='3d_rectangle_8086el' ;%'3d_rectangle_8086el'; %'3d_rectangle_8086el'
%'3d_rectangle_1272el';%'3d_rectangle_1272el';%'3d_rectangle_660el'; 
kActu = 1.0;    % multiplicative factor for the actuation forces, dependent on the mesh
[MeshNominal, nodes, elements, nsetBC, esetBC] = create_mesh(filename, myElementConstructor, propRigid);
[Lx, Ly, Lz] = mesh_dimensions(nodes);


% plot nominal mesh [Can be removed]
elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
figure('units','normalized','position',[.2 .1 .6 .8])
PlotMeshAxis(nodes, elementPlot, 0);
hold off

%% DEFINE SHAPE VARIATIONS ________________________________________________
[y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, x_concaveTail] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);

U = [z_tail, z_head, y_thinFish, ... %z_tail
    y_linLongTail, y_head,y_ellipseFish, ... % _linLongTail
    z_smallFish, z_notch,z_notch, x_concaveTail];

% U = [z_tail,z_tail,z_tail,z_tail,z_tail,z_tail,z_tail,z_tail,z_tail];


%% PLOT SHAPE VARIATIONS __________________________________________________
f1 = figure('units','centimeters','position',[3 3 19 10.8]);
elementPlot = elements(:,1:4); hold on
L = [Lx,Ly,Lz];
O = [-Lx,-Ly/2,-Lz/2];
textPosX = -0.2;
textPosY = -0.05;
textPosZ = -0.14;
textPos = [textPosX, textPosY, textPosZ];

% Define tiled layout with tighter spacing
t = tiledlayout(3, 6, 'TileSpacing', 'none', 'Padding', 'tight');
n_shapeVariations = 9;
n_subplot = 2 * n_shapeVariations;

% Preallocate an array for the axes handles
axesHandles = gobjects(1, n_subplot); 
widthVec = zeros(n_shapeVariations,1);
heightVec = zeros(n_shapeVariations,1);
xVec = zeros(n_shapeVariations,1);
yVec = zeros(n_shapeVariations,1);
yCenterVec = zeros(n_shapeVariations,1);
xCenterVec = zeros(n_shapeVariations,1);

for i = 1:n_shapeVariations
    % Shape variation i, positive
    axesHandles((i-1)*2 + 1) = nexttile(t, (i-1)*2 + 1);
    create_subplot(U(:,i), 0.5, num2str(i), nodes, elementPlot, L, O, textPos);
    
    % Shape variation i, negative
    axesHandles((i-1)*2 + 2) = nexttile(t, (i-1)*2 + 2);
    create_subplot(U(:,i), -0.5, num2str(i), nodes, elementPlot, L, O, textPos);
    
    % Get the positions of the two subplots
    pos1 = get(axesHandles((i-1)*2 + 1), 'Position');
    pos2 = get(axesHandles((i-1)*2 + 2), 'Position');
    
    % Calculate the rectangle position that encloses both subplots
    x = min(pos1(1), pos2(1)); % Minimum x position
    y = min(pos1(2), pos2(2)); % Minimum y position
    width = max(pos1(1) + pos1(3), pos2(1) + pos2(3)) - x;  % Total width
    height = max(pos1(2) + pos1(4), pos2(2) + pos2(4)) - y; % Total height
    
    % store the rectangle position 
    xVec(i) = x;
    yVec(i) = y;
    yCenterVec(i) = y + height/2;
    xCenterVec(i) = x + width/2;
    widthVec(i) = width;
    heightVec(i) = height;
    
end

% Set axis limits for all the axes
axis(axesHandles, [-0.40 0 -0.061 0.061 -0.16 0.16]);

% Find boxes size and placement
x1All = min(xVec(1:3:end));
x2All = min(xVec(2:3:end));
x3All = min(xVec(3:3:end));
xAll = [x1All; x2All; x3All];
y1All = min(yVec(1:3));
y2All = min(yVec(4:6));
y3All = min(yVec(7:9));
yCenterAll = [yCenterVec(1);yCenterVec(4);yCenterVec(7)];
xCenterAll = [xCenterVec(1);xCenterVec(2);xCenterVec(3)];

y1MaxAll = max(yMaxVec(1:3));
y2MaxAll = max(yMaxVec(4:6));
y3MaxAll = max(yMaxVec(7:9));

yAll = [y1All; y2All; y3All];
yMaxAll = [y1MaxAll,y2MaxAll,y3MaxAll];
widthAll = max(widthVec);
heightAll = max(heightVec);

% Define additional box quantity
extraHorSpace = 0.005 ;%0.02;
extraVertSpace = 0.02;% 0.01;
widthAll = widthAll + extraHorSpace;
heightAll = heightAll + extraVertSpace;


for i = 1:n_shapeVariations
    % Draw a rectangle around the two subplots
    xIdx = mod(i-1,3);
    yIdx = ceil(i/3);
    xBox = xAll(xIdx+1);
    yBox = yCenterAll(yIdx)-heightAll/2;
    annotation('rectangle', [xBox-extraHorSpace/2 yBox-extraVertSpace/2, ...
        widthAll+extraHorSpace,heightAll+extraVertSpace], ...
        'Color', 'k', 'LineWidth', 1.0, 'LineStyle', '--');
    
    varName = strcat('$$\xi_',num2str(i),'$$');
    annotation('textbox',[xBox-extraHorSpace/2 yBox-extraVertSpace/2+0.18, 0.1,0.1], ...
        'String',varName, 'FontSize', 14,...
        'LineStyle','none','FitBoxToText','on','interpreter','latex')
    
end


% % Draw a few lines to help designing the figure
% for i = 1:n_shapeVariations
%     xIdx = mod(i-1,3);
%     yIdx = ceil(i/3);
%     xBox = xAll(xIdx+1);
%     yBox = (yMaxAll(yIdx)-yAll(yIdx));
%     annotation('line', [xVec(i) xVec(i)+0.05], ...
%         [yVec(i),yVec(i)], ...
%         'Color', 'blue');
%     annotation('line', [xVec(i) xVec(i)+0.05], ...
%         [yVec(i)+heightVec(i),yVec(i)++heightVec(i)], ...
%         'Color', 'blue');
%     annotation('line', [xVec(i) xVec(i)+0.25], ...
%         [yCenterAll(yIdx),yCenterAll(yIdx)], ...
%         'Color', 'red');
% end



% annotation('line', [0,1], [0.5+heightAll/2,0.5+heightAll/2],'color', 'blue')
% annotation('line', [0,1], [0.5-heightAll/2,0.5-heightAll/2],'color', 'blue')

% exportgraphics(f1,'shape_variations.pdf','Resolution',1200)
% exportgraphics(f1,'SO3_shapes_V1.jpg','Resolution',600)


%%
function create_subplot(subU,xiPlot,nr_sv,nodes,elementPlot,L,O, textPos)
    v = reshape(subU*xiPlot, 3, []).';
    dm = PlotFieldonDeformedMesh(nodes, elementPlot, v, 'factor', 1);
    plotcube(L,O,.05,[0 0 0]);
end

