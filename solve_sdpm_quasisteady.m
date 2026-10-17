function result = solve_sdpm_quasisteady(cfg, bodyState)
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
%  6) BUILD CURRENT QUASI-STEADY WAKE
% =========================================================================

wakeCfg = cfg.wake;

wakeCfg.length = ...
    wakeCfg.lengthChord .* cfg.airfoil.chord;

wake = build_steady_wake( ...
    body, ...
    flow, ...
    wakeCfg);


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

inflBW = wake_influence( ...
    body, ...
    wake.panels);

phiWake = inflBW.A   * wake.mu;
VtWake  = inflBW.BtJ * wake.mu;
VnWake  = inflBW.BnJ * wake.mu;
%% ========================================================================
%  9) ASSEMBLE AND SOLVE
% =========================================================================

[A,b] = assemble_dirichlet_system( ...
    inflBB, ...
    inflBW, ...
    sigma);

mu = A \ b;


%% ========================================================================
%  10) POSTPROCESS
% =========================================================================

aero = postprocess_aero( ...
    body, ...
    mu, ...
    sigma, ...
    inflBB, ...
    inflBW, ...
    flow, ...
    cfg);


%% ========================================================================
%  11) PACKAGE RESULTS
% =========================================================================

result.cfg = cfg;

result.bodyState = bodyState;

result.flow = flow;

result.Bp0 = Bp0;

result.body = body;

result.Vbody = Vbody;

result.wake = wake;

result.sigma = sigma;

result.mu = mu;

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