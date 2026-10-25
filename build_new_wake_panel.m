function newWakePanel = build_new_wake_panel(body, flow, wakeOld, dt)
% BUILD_NEW_WAKE_PANEL
%
% Construct geometry of a single trailing-edge wake panel to enforce kutta
% condition
%
% New panel begins and CURRENT trailing edge
%
% If an old wake exists:
%           - youngest old wake node is treated as TE
%
% If no wake exists yet, first panel initialized using free-stream
% convection over one time step:
%           - current TE = current TE + U_inf*dt

% ===== INPUTS =====
% ---> body = current body panel geometry
% ---> flow = flow struct containing flow.Uvec (free stream velocity)
% ---> wakeOld = persistent wake from previous timestep
% ---> dt = current timestep
%
%
% ===== OUTPUTS =====
% ---> newWakePanel = panel geometry structure for exactly one new wake
% panel

%% ========== 1) LOAD CURRENT TRAILING EDGE ==========

% same convention as steady SDPM solver
TE = body.Bp(end,:);

%% ========== 2) DETERMINE DOWNSTREAM END OF NEW PANEL ==========

if isempty(wakeOld.Bp)

    % First timestep:
    %   - no previous wake node to connect to, so initialize wake panel
    %     with a length corresponding to one timestep of free-stream
    %     convection

    wakePoint = TE + flow.Uvec .* dt;

else

    % Existing wake:
    %   - Connect current TEE edge to youngest existing wake node

    wakePoint = wakeOld.Bp(1,:);

end

%% ========== 3) CHECK NEW PANEL LENGTH ==========

panelLength = norm(wakePoint-TE);

if panelLength <= 1e-12
    error(['New panel has essentialy zero length.', ...
        'The existing wake must be convected away from the ' ...
        'trailing edge before another weake panel can be shed'])
end

%% ========== 4) BUILD PANEL GEOMETRY ==========
BpNew = [TE; wakePoint];

newWakePanel = panel_geometry(BpNew);

end