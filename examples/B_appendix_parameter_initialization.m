% ------------------------------------------------------------------------ 
% B_appendix_parameter_initialization.m
%
% Description: Analyze the results obtain for SO1-SO3 under different
% initializuation parameters. The optimization itself is performed in
% B_shape_optimisation.m
% 
% Last modified: 08/02/2026, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
clear; 
close all; 
clc
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 


%% OPTIMIZATION PARAMETERS
kActu = 6.0*1e4;   
n_elements = 8086;

%% RETRIEVE RESULTS _______________________________________________________

% main results
filename = sprintf('Results/Data/B_shape_optimization/SO1_%d_el_kActu_%.3f.mat', ...
                    n_elements, kActu);
load(filename)
filename = sprintf('Results/Data/B_shape_optimization/SO2_%d_el_kActu_%.3f.mat', ...
                    n_elements, kActu);
load(filename)
filename = sprintf('Results/Data/B_shape_optimization/SO3_%d_el_kActu_%.3f.mat', ...
                    n_elements, kActu);
load(filename)

outmain_1 = out1;
outmain_2 = out2;
outmain_3 = out3;

% init 1
filename = sprintf('Results/Data/Appendix/Param_init/SO1_pInit1__%d_el_kActu_%.3f.mat', ...
                    n_elements, kActu);
load(filename)
filename = sprintf('Results/Data/Appendix/Param_init/SO2_pInit1__%d_el_kActu_%.3f.mat', ...
                    n_elements, kActu);
load(filename)
filename = sprintf('Results/Data/Appendix/Param_init/SO3_pInit1__%d_el_kActu_%.3f.mat', ...
                    n_elements, kActu);
load(filename)

out1_1 = out1;
out1_2 = out2;
out1_3 = out3;

% init 2
filename = sprintf('Results/Data/Appendix/Param_init/SO1_pInit2__%d_el_kActu_%.3f.mat', ...
                    n_elements, kActu);
load(filename)
filename = sprintf('Results/Data/Appendix/Param_init/SO2_pInit2__%d_el_kActu_%.3f.mat', ...
                    n_elements, kActu);
load(filename)
filename = sprintf('Results/Data/Appendix/Param_init/SO3_pInit2__%d_el_kActu_%.3f.mat', ...
                    n_elements, kActu);
load(filename)

out2_1 = out1;
out2_2 = out2;
out2_3 = out3;

%% PRINT RESULTS AND PERFORMANCE STATS ____________________________________
SOIdx = 3;
initIdx = 2;
switch initIdx
    case 0
        switch SOIdx
            case 1
                outPrint = outmain_1;  
            case 2
                outPrint = outmain_2;
            case 3
                outPrint = outmain_3;
        end      
    case 1
        switch SOIdx
            case 1
                outPrint = out1_1;  
            case 2
                outPrint = out1_2;
            case 3
                outPrint = out1_3;
        end    
    case 2
        switch SOIdx
            case 1
                outPrint = out2_1;  
            case 2
                outPrint = out2_2;
            case 3
                outPrint = out2_3;
        end    
end
xiStar = outPrint.xiStar;
tOpti = outPrint.tOpti;
nIt = outPrint.nIt;
tPerIt = outPrint.tPerIt; 
fprintf('Computation time: %.2f min\n', tOpti);
fprintf('Number of built models and solved EoMs: %d\n', nIt);
fprintf('Computation time per built models/EoMs: %.2f min\n', tPerIt)
disp('Optimal parameters:')
disp(xiStar)
%% PLOT COST FUNCTION FOR THE 3 EXPERIMENTS _______________________________
f_cost = figure('units','centimeters','position',[3 3 9 6]);
hold on
plot(outmain_1.LwoBEvo,'LineWidth',1)
plot(out1_1.LwoBEvo,'--','LineWidth',1)
plot(out2_1.LwoBEvo,'-.','LineWidth',1)
grid on
ylabel('$$L$$')
xlabel('Iterations')
legend('main','1','2')
hold off
fig_title = sprintf('Results/Figures/Appendix/Param_Init/SO1_cost_evolution_%d_el.pdf', n_elements);
exportgraphics(f_cost,fig_title,'Resolution',1200)

f_cost = figure('units','centimeters','position',[3 3 9 6]);
hold on
plot(outmain_2.LwoBEvo,'LineWidth',1)
plot(out1_2.LwoBEvo,'--','LineWidth',1)
plot(out2_2.LwoBEvo,'-.','LineWidth',1)
grid on
ylabel('$$L$$')
xlabel('Iterations')
legend('main','1','2')
hold off
fig_title = sprintf('Results/Figures/Appendix/Param_Init/SO2_cost_evolution_%d_el.pdf', n_elements);
exportgraphics(f_cost,fig_title,'Resolution',1200)

f_cost = figure('units','centimeters','position',[3 3 9 6]);
hold on
plot(outmain_3.LwoBEvo,'LineWidth',1)
plot(out1_3.LwoBEvo,'--','LineWidth',1)
plot(out2_3.LwoBEvo,'-.','LineWidth',1)
grid on
ylabel('$$L$$')
xlabel('Iterations')
legend('main','1','2')
hold off
fig_title = sprintf('Results/Figures/Appendix/Param_Init/SO3_cost_evolution_%d_el.pdf', n_elements);
exportgraphics(f_cost,fig_title,'Resolution',1200)


