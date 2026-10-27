function [result, state] = solve_sdpm_unsteady(cfg, bodyState, state, dt)
%SOLVE_SDPM_QUASISTEADY
% Solve one instantaneous quasi-steady moving-body SDPM state.
%
% INPUTS
%
%   cfg
%       Simulation configuration
%
%   bodyState
%       Current rigid-body state:
%
%           bodyState.x
%           bodyState.z
%           bodyState.theta
%
%           bodyState.u
%           bodyState.w
%           bodyState.omega
%% ========================================================================
%  0) PERSISTENT UNSTEADY STATE
% =========================================================================

% These variables contain aerodynamic information carried from the
% previous timestep.

if ~isfield(state,'time')
    state.time = 0;
end

% wake properties from previous time step
if ~isfield(state,'wake')
    state.wake.Bp = []; % wake panel boundary points from previous timestep
    state.wake.mu = []; % wake panel strengths from previous timestep
end

% bound vorticity from previous timestep
if ~isfield(state,'GammaBoundPrev')
    state.GammaBoundPrev = [];
end

% panel potentials (evaluated per collocation method) from previous time step
if ~isfield(state,'phiPrev')
    state.phiPrev = [];
end

%% ========================================================================
%  1) FLOW
% =========================================================================

flow = cfg.flow;

flow.alphaRad = deg2rad(flow.alphaDeg);

flow.Uvec = flow.U .* ...
    [cos(flow.alphaRad), sin(flow.alphaRad)];


%% ========================================================================
%  2) REFERENCE AIRFOIL GEOMETRY
% =========================================================================

% Body-frame / reference coordinates

Bp0 = naca4_airfoil( ...
    cfg.airfoil.code, ...
    cfg.airfoil.chord, ...
    cfg.airfoil.nPanels, ...
    cfg.airfoil.spacing, ...
    cfg.airfoil.closedTE);


%% ========================================================================
%  3) TRANSFORM BODY INTO CURRENT POSITION
% =========================================================================

Bp = transform_body( ...
    Bp0, ...
    bodyState);


%% ========================================================================
%  4) BUILD CURRENT PANEL GEOMETRY
% =========================================================================

body = panel_geometry(Bp);


%% ========================================================================
%  5) CALCULATE LOCAL BODY VELOCITY
% =========================================================================

% Returns N x 2 velocity at each body collocation point.

Vbody = body_velocity( ...
    body, ...
    bodyState);


%% ========================================================================
%  6) LOAD PERSISTENT WAKE
% =========================================================================
wakeOld = state.wake;

% Convect all wake elements shed during previous timesteps
if ~isempty(wakeOld.Bp)

    wakeOld = convect_wake( ...
        wakeOld, ...
        flow, ...
        dt);

end
% wakeCfg = cfg.wake;
% 
% wakeCfg.length = ...
%     wakeCfg.lengthChord .* cfg.airfoil.chord;
% 
% wake = build_steady_wake( ...
%     body, ...
%     flow, ...
%     wakeCfg);


%% ========================================================================
%  7) MOVING-BODY SOURCE STRENGTHS
% =========================================================================

sigma = source_strengths( ...
    body, ...
    flow, ...
    Vbody);


%% ========================================================================
%  8) INFLUENCE COEFFICIENTS
% =========================================================================



inflBB = body_influence(body);

% Influence from OLD persistent wake

if ~isempty(wakeOld.mu)

    wakeOld.panels = panel_geometry(wakeOld.Bp);

    inflBWold = wake_influence(body,wakeOld.panels);

    phiWakeOld = inflBWold.A * wakeOld.mu;

        VtWakeOld = ...
        inflBWold.BtJ * wakeOld.mu;

            VnWakeOld = ...
        inflBWold.BnJ * wakeOld.mu;

else

    phiWakeOld = zeros(body.nPanels,1);
    VtWakeOld  = zeros(body.nPanels,1);
    VnWakeOld  = zeros(body.nPanels,1);

end


% ---------------------------------------------------------
% Construct geometry of ONE candidate new TE wake panel
% ---------------------------------------------------------

newWakePanel = build_new_wake_panel(body,flow, wakeOld, dt);

inflBWnew = wake_influence(body,newWakePanel);



%% ========================================================================
%  9) ASSEMBLE AND SOLVE
% =========================================================================

[A,b] = assemble_dirichlet_system( ...
    inflBB, ...
    inflBWnew, ...
    sigma, ...
    phiWakeOld);

muSolve = A \ b;
n =  body.nPanels;
muBody = muSolve(1:n);
muNew  = muSolve(end);

%% ========================================================================
%  10) COMPLETE CURRENT WAKE CONTRIBUTIONS
% =========================================================================
phiWakeNew = inflBWnew.A * muNew;
VtWakeNew = inflBWnew.BtJ * muNew;
VnWakeNew = inflBWnew.BnJ * muNew;

phiWake = phiWakeOld + phiWakeNew;
VtWake = VtWakeOld + VtWakeNew;
VnWake = VnWakeOld + VnWakeNew;

%% ========================================================================
%  11) CURRENT BODY SURFACE POTENTIAL
% =========================================================================

phiDoublet = inflBB.A * muBody;

phiSource = inflBB.SourceInfluence * sigma;

phiBody = phiDoublet ...
  + phiSource ...
  + phiWake;

%% ========================================================================
%  12) APPEND NEW WAKE PANEL
% =========================================================================
wake = shed_wake(...
    wakeOld, ...
    newWakePanel, ...
    muNew);

%% ========================================================================
%  13) POSTPROCESS
% =========================================================================

aero = postprocess_aero( ...
    body, ...
    muBody, ...
    sigma, ...
    inflBB, ...
    VtWake, ...
    flow, ...
    cfg);


%% ========================================================================
%  14) UPDATE PERSISTENT STATE
% =========================================================================

state.wake = wake;
state.GammaBoundPrev = aero.Gamma;

state.phiPrev = phiBody;

state.time = state.time + dt;
%% ========================================================================
%  15) PACKAGE RESULTS
% =========================================================================

result.cfg = cfg;

result.bodyState = bodyState;

result.flow = flow;

result.Bp0 = Bp0;

result.body = body;

result.Vbody = Vbody;

result.wake = wake;

result.sigma = sigma;

result.mu = muSolve;
result.muBody = muBody;
result.muNew = muNew;

result.system.A = A;

result.system.b = b;

result.aero = aero;

result.Cp = aero.Cp;
result.Vt = aero.Vt;
result.Gamma = aero.Gamma;

result.CL_gamma = aero.CL_gamma;
result.CL_pressure = aero.CL_pressure;

result.CD = aero.CD_pressure;
result.CM = aero.CM_pressure;

result.xCp = ...
    body.cp(:,1) ./ cfg.airfoil.chord;


end