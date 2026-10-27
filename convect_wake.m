function wake = convect_wake(wake,flow,dt)

%CONVECT_WAKE
% Convect all existing wake nodes using freestream velocity only.
%
% First-order explicit wake-convection model:
%
%       x_w^(n+1) = x_w^n + U_inf*dt
%
% INPUTS
%
%   wake
%       Persistent wake structure
%
%       wake.Bp
%           Wake boundary-point coordinates, N x 2
%
%       wake.mu
%           Wake-panel doublet strengths
%
%   flow
%       Flow structure containing flow.Uvec
%
%   dt
%       Timestep
%
% OUTPUT
%
%   wake
%       Wake structure with updated boundary-point locations


%% ========================================================================
%  1) EMPTY WAKE CHECK
% =========================================================================

if isempty(wake.Bp)
    return
end


%% ========================================================================
%  2) FREESTREAM DISPLACEMENT
% =========================================================================

dX = flow.Uvec .* dt;


%% ========================================================================
%  3) CONVECT ALL WAKE NODES
% =========================================================================

wake.Bp(:,1) = wake.Bp(:,1) + dX(1);

wake.Bp(:,2) = wake.Bp(:,2) + dX(2);


end