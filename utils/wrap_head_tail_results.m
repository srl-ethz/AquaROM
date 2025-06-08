% ------------------------------------------------------------------------ 
% wrap_head_tail_results.m
%
% Description: Wrap head and tail results from the different models in a
% struct
%
% Last modified: 08/06/2025, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
function [uHead_sol, uTail_sol] = wrap_head_tail_results(...
    uHead_FOM, ...
    uHead_ROM, ...
    uHead_PROM_3, ...
    uHead_PROM_5, ...
    uHead_PROM_8, ...
    uTail_FOM, ...
    uTail_ROM, ...
    uTail_PROM_3, ...
    uTail_PROM_5, ...
    uTail_PROM_8)

    % Store Head results
    uHead_sol.uHead_FOM     = uHead_FOM;
    uHead_sol.uHead_ROM     = uHead_ROM;
    uHead_sol.uHead_PROM_3  = uHead_PROM_3;
    uHead_sol.uHead_PROM_5  = uHead_PROM_5;
    uHead_sol.uHead_PROM_8  = uHead_PROM_8;

    % Store Tail results
    uTail_sol.uTail_FOM     = uTail_FOM;
    uTail_sol.uTail_ROM     = uTail_ROM;
    uTail_sol.uTail_PROM_3  = uTail_PROM_3;
    uTail_sol.uTail_PROM_5  = uTail_PROM_5;
    uTail_sol.uTail_PROM_8  = uTail_PROM_8;
end
