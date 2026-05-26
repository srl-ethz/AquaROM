% ------------------------------------------------------------------------ 
% B_appendix_grid_search_estimation.m
%
% Description: Estimates the time needed to perform a grid search
% optimization using the ROM.
% 
% Last modified: 06/02/2026, Mathieu Dubied, ETH Zurich
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


%% GRID PARAMETRIZATION ___________________________________________________
SOIdx =   3;
switch SOIdx
    case 1
        % Constraint for SO1 (upper bound, -lower bound)
        nParam = 3;
        b = [0.5;0.5;
            0.5;0.5;
            0.15;0.15];
    case 2
        % Constraint for SO2 (upper bound, -lower bound)
        nParam = 5;
        b = [0.5;0.5;
            0.5;0.5;
            0.4;0.4;
            0.5;0.5;
            0.5;0.5]; 
    case 3
        % Constraint for SO3 (upper bound, -lower bound)
        nParam = 8;
        b = [0.5;0.5;
            0.5;0.5;
            0.3;0.3;
            0.4;0.4;
            0.3;0.3;
            0.3;0.3 ;
            0.3;0.3;
            0.2;0.01];
end

A = zeros(2 * nParam, nParam);
for i = 1:nParam
    A(2*i-1:2*i,i) =[1;-1];
end

% Define search grid 
granularity = 0.15; % rebuild criterion of the PROM pipeline
[xi_lists, step_used] = build_param_grid_lists(A, b, granularity);


%% GRID STATS _____________________________________________________________
[nComb, nCorners] = grid_combo_stats(xi_lists);
timePerSolve = 31; %in sec
fprintf('\nGrid of SO%d (%d parameters):\n', SOIdx, nParam);
fprintf('   - N combinations: %d\n', nComb);
fprintf('   - N corners of hyperbox: %d\n', nCorners);
fprintf('   - Step granularity: %.2f\n', step_used);

timeAllCombi = nComb*timePerSolve/60;
timeAllCorners = nCorners*timePerSolve/60;

fprintf('   - Time for all combination: %.2f hours\n', timeAllCombi/60);
fprintf('   - Time for all corners: %.2f min\n', timeAllCorners);


%% UTILS __________________________________________________________________
function [xi_lists, step_used] = build_param_grid_lists(A, b, granularity)
%BUILD_PARAM_GRID_LISTS Build per-parameter 1D grids from axis-aligned A*xi<=b bounds.
%
% Returns:
%   xi_lists   : 1×n cell, xi_lists{i} is the vector of grid values for parameter i
%   step_used  : n×1 actual step size used (<= granularity unless range=0)
%
% Extremes are always included. If granularity does not divide the range
% exactly, a finer step is chosen.

    validateattributes(A, {'numeric'}, {'2d','nonempty'});
    validateattributes(b, {'numeric'}, {'column','numel',size(A,1)});
    validateattributes(granularity, {'numeric'}, {'scalar','positive','finite'});

    n  = size(A,2);
    lb = -inf(n,1);
    ub =  inf(n,1);

    % Infer bounds per variable (axis-aligned constraints)
    for r = 1:size(A,1)
        nz = find(A(r,:) ~= 0);
        if numel(nz) ~= 1
            error("Row %d of A is not axis-aligned (expected exactly 1 nonzero).", r);
        end
        j = nz(1);
        a = A(r,j);
        if a > 0
            ub(j) = min(ub(j), b(r)/a);   % a*xi_j <= b
        else
            lb(j) = max(lb(j), b(r)/a);   % a<0 flips inequality
        end
    end

    if any(~isfinite(lb)) || any(~isfinite(ub)) || any(lb > ub)
        error("Invalid or infeasible bounds inferred from A and b.");
    end

    xi_lists  = cell(1,n);
    step_used = zeros(n,1);

    % Build grids
    for j = 1:n
        range = ub(j) - lb(j);
        if range == 0
            xi_lists{j} = lb(j);
            step_used(j) = 0;
        else
            nIntervals   = floor(range / granularity);
            step_used(j) = range / nIntervals;
            xi_lists{j}  = linspace(lb(j), ub(j), nIntervals + 1);
        end
    end
end

function [nComb, nCorners] = grid_combo_stats(xi_lists)
% Number of grid combinations and hyperbox corners.
%
% Inputs:
%   xi_lists : 1×n cell, xi_lists{i} is the vector of values for parameter i
%
% Outputs:
%   nComb    : total number of combinations = prod_i numel(xi_lists{i})
%   nCorners : number of hyperbox corners    = 2^n  (assuming each dimension has >=2 values)
%
% Notes:
% - If some dimension has only 1 value, it is degenerate; corners reduce accordingly:
%   nCorners = 2^(#dims with >=2 values)

    if ~iscell(xi_lists) || isempty(xi_lists)
        error("xi_lists must be a non-empty cell array.");
    end

    n = numel(xi_lists);

    lens = zeros(1,n);
    for i = 1:n
        if ~isnumeric(xi_lists{i}) || isempty(xi_lists{i})
            error("xi_lists{%d} must be a non-empty numeric vector.", i);
        end
        lens(i) = numel(xi_lists{i});
    end

    nComb = prod(lens);

    nActive = sum(lens >= 2);   % non-degenerate dimensions
    nCorners = 2^nActive;
end

