function aero = postprocess_aero( ...
    body, muBody, sigma, inflBB, VtWake, flow, cfg)
%POSTPROCESS_AERO
% Compute surface velocity, pressure coefficient, circulation, and
% integrated aerodynamic loads for the steady SDPM solution.
%
% INPUTS
%   body
%       Body panel geometry from panel_geometry()
%
%   mu
%       Doublet strengths:
%
%           mu(1:N) = body doublet strengths
%           mu(end) = steady wake doublet strength
%
%   sigma
%       Prescribed body-panel source strengths
%
%   inflBB
%       Body-body influence coefficients
%
%   inflBW
%       Wake-body influence coefficients
%
%   flow
%       Flow structure containing:
%
%           flow.U
%           flow.alphaRad
%           flow.rho
%
%   cfg
%       Solver configuration structure
%
%
% OUTPUT
%   aero
%       Structure containing surface velocity, Cp, circulation,
%       dimensional loads, and nondimensional coefficients.


%% ========================================================================
%  1) BASIC QUANTITIES
% =========================================================================

n = body.nPanels;

U     = flow.U;
alpha = flow.alphaRad;
rho   = flow.rho;

c = cfg.airfoil.chord;

% Separate body and wake doublet strengths
n = body.nPanels;

muBody = muBody(:);
VtWake = VtWake(:);


%% ========================================================================
%  2) WAKE-INDUCED TANGENTIAL VELOCITY
% =========================================================================

% All steady wake panels presently share the same doublet strength.
%
% inflBW.BtJ is:
%
%       body receiver panel x wake observer panel
%
% Summing across wake panels gives the total influence coefficient
% corresponding to the single common wake strength.

Vt_wake = VtWake;


%% ========================================================================
%  3) BODY DOUBLET TANGENTIAL VELOCITY
% =========================================================================

Vt_doublet = ...
    inflBB.BtJ * muBody;


%% ========================================================================
%  4) BODY SOURCE TANGENTIAL VELOCITY
% =========================================================================

Vt_source = ...
    inflBB.BtG * sigma;


%% ========================================================================
%  5) FREESTREAM TANGENTIAL VELOCITY
% =========================================================================

Vt_freestream = ...
    U .* cos(body.theta - flow.alphaRad);


%% ========================================================================
%  6) NUMERICAL DOUBLET-STRENGTH DERIVATIVE CORRECTION
% =========================================================================

% IMPORTANT:
% This formulation is intentionally preserved from the original solver
% during the structural refactor.
%
% It should be audited separately for nonuniform/cosine panel spacing,
% because mu is associated with panel locations and the denominator should
% ultimately represent the correct distance between those locations.

dmu_ds = zeros(n,1);


% Interior panels
dmu_ds(2:n-1) = ...
    (muBody(3:end) - muBody(1:end-2)) ./ ...
    (body.length(2:end-1) + body.length(1:end-2));


% Leading endpoint
dmu_ds(1) = ...
    (muBody(2) - muBody(1)) ./ ...
    body.length(1);


% Trailing endpoint
dmu_ds(end) = ...
    (muBody(end) - muBody(end-1)) ./ ...
    body.length(end);


%% ========================================================================
%  7) TOTAL SURFACE TANGENTIAL VELOCITY
% =========================================================================

Vt = ...
      Vt_doublet ...
    + Vt_wake ...
    + Vt_source ...
    + Vt_freestream ...
    + 0.5 .* dmu_ds;


%% ========================================================================
%  8) PRESSURE COEFFICIENT
% =========================================================================

Cp = ...
    1 - (Vt ./ U).^2;


%% ========================================================================
%  9) CIRCULATION-BASED LIFT
% =========================================================================

% Circulation obtained from integration of surface tangential velocity.
Gamma = ...
    sum(Vt .* body.length);


% General nondimensional circulation-based lift coefficient:
%
%             2 Gamma
%       CL = ---------
%              U c
%
CL_gamma = ...
    2 .* Gamma ./ (U .* c);


%% ========================================================================
%  10) DYNAMIC PRESSURE
% =========================================================================

qInf = ...
    0.5 .* rho .* U.^2;


%% ========================================================================
%  11) PANEL PRESSURE FORCES
% =========================================================================

% panel_geometry() already supplies the outward unit normal:
%
%       n = [-sin(theta), cos(theta)]
%
% for the clockwise body-panel convention.

nx = body.normal(:,1);
nz = body.normal(:,2);

ds = body.length;


% Pressure force per unit span on each panel:
%
%       dF = -q_inf Cp n ds

dFx = ...
    -qInf .* Cp .* nx .* ds;

dFz = ...
    -qInf .* Cp .* nz .* ds;


%% ========================================================================
%  12) GLOBAL FORCE COMPONENTS
% =========================================================================

Fx = sum(dFx);

Fz = sum(dFz);


%% ========================================================================
%  13) RESOLVE INTO LIFT / DRAG AXES
% =========================================================================

% Freestream-aligned drag direction:
%
%       e_D = [cos(alpha), sin(alpha)]
%
% Lift direction:
%
%       e_L = [-sin(alpha), cos(alpha)]

D = ...
      Fx .* cos(alpha) ...
    + Fz .* sin(alpha);


L = ...
     -Fx .* sin(alpha) ...
    + Fz .* cos(alpha);


%% ========================================================================
%  14) MOMENT ABOUT REFERENCE POINT
% =========================================================================

% Reference location is stored nondimensionally in cfg and converted here
% to dimensional coordinates.

xRef = ...
    cfg.reference.xMomentChord .* c;

zRef = ...
    cfg.reference.zMomentChord .* c;


% Position from reference point to panel collocation point
rx = ...
    body.cp(:,1) - xRef;

rz = ...
    body.cp(:,2) - zRef;


% 2-D moment per unit span about the y-axis:
%
%       dM = r_x dF_z - r_z dF_x

dM = ...
    rx .* dFz ...
    - rz .* dFx;


M = sum(dM);


%% ========================================================================
%  15) NONDIMENSIONAL LOAD COEFFICIENTS
% =========================================================================

CL_pressure = ...
    L ./ (qInf .* c);


CD_pressure = ...
    D ./ (qInf .* c);


CM_pressure = ...
    M ./ (qInf .* c.^2);


%% ========================================================================
%  16) OPTIONAL CONSTANT-PRESSURE DEBUG CHECK
% =========================================================================

% A spatially uniform pressure acting around a closed contour should
% contribute approximately zero net force.
%
% Keeping this is useful as a geometry / normal-orientation diagnostic.

dFx_const = ...
    -qInf .* nx .* ds;

dFz_const = ...
    -qInf .* nz .* ds;

Fx_const = ...
    sum(dFx_const);

Fz_const = ...
    sum(dFz_const);


%% ========================================================================
%  17) PACKAGE SURFACE QUANTITIES
% =========================================================================

aero.Vt = Vt;

aero.Cp = Cp;

aero.dmu_ds = dmu_ds;


% Individual surface-velocity contributions
aero.Vt_doublet = Vt_doublet;

aero.Vt_wake = Vt_wake;

aero.Vt_source = Vt_source;

aero.Vt_freestream = Vt_freestream;


%% ========================================================================
%  18) PACKAGE CIRCULATION / COEFFICIENTS
% =========================================================================

aero.Gamma = Gamma;

aero.CL_gamma = CL_gamma;

aero.CL_pressure = CL_pressure;

aero.CD_pressure = CD_pressure;

aero.CM_pressure = CM_pressure;


%% ========================================================================
%  19) PACKAGE DIMENSIONAL LOADS
% =========================================================================

% Panel force distributions
aero.dFx = dFx;

aero.dFz = dFz;

aero.dM = dM;


% Integrated global forces
aero.Fx = Fx;

aero.Fz = Fz;


% Lift / drag / moment
aero.L = L;

aero.D = D;

aero.M = M;


%% ========================================================================
%  20) PACKAGE DEBUG / REFERENCE DATA
% =========================================================================

aero.Fx_const = Fx_const;

aero.Fz_const = Fz_const;

aero.qInf = qInf;

aero.cRef = c;

aero.xRef = xRef;

aero.zRef = zRef;


end