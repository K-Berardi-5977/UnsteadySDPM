% ===== QUASI-STEADY / UNSTEADY SDPM ROOT ===== %

clc
clear


%% ========================================================================
%  1) CONFIGURATION SELECTION
% =========================================================================

choice = questdlg( ...
    'How would you like to configure the SDPM solver?', ...
    'SDPM Setup', ...
    'Use Defaults', ...
    'Enter Inputs', ...
    'Cancel', ...
    'Use Defaults');


switch choice

    case 'Use Defaults'

        cfg = default_sdpm_config();

    case 'Enter Inputs'

        cfg = prompt_sdpm_config();

    otherwise

        disp('SDPM execution cancelled.')
        return

end


%% ========================================================================
%  2) TIME DEFINITION
% =========================================================================

% For now these are prescribed here.
% Later these can move into cfg.time.

% dt   = 0.01;
% tEnd = 1;
% 
% t = 0:dt:tEnd;
% 
% Nt = numel(t);


%% ========================================================================
%  3) PRESCRIBED BODY MOTION
% =========================================================================

% Example:
%
% Harmonic pitching motion
%
%       theta(t) = thetaAmp*sin(omegaMotion*t)
%
%       omega(t) = dtheta/dt
%                = thetaAmp*omegaMotion*cos(omegaMotion*t)

thetaAmp = deg2rad(5);

freq = 0.05;                   % Hz
redFreq = pi*freq;
omegaMotion = 2*pi*freq;

T = 1/freq;

nCycles = 2;

dt = T/200;
tEnd = nCycles*T;

t = 0:dt:tEnd;

Nt = numel(t);

% ----- Position histories -----

bodyState.x = zeros(1,Nt);

bodyState.z = zeros(1,Nt);

bodyState.theta = ...
    thetaAmp .* sin(omegaMotion .* t);


% ----- Velocity histories -----

bodyState.u = zeros(1,Nt);

bodyState.w = zeros(1,Nt);

% bodyState.omega = zeros(1, Nt);

bodyState.omega = ...
    thetaAmp .* omegaMotion .* cos(omegaMotion .* t);


%% ========================================================================
%  4) PREALLOCATE TIME-HISTORY STORAGE
% =========================================================================

results = cell(1,Nt);

CL = zeros(1,Nt);

CM = zeros(1,Nt);

Gamma = zeros(1,Nt);


%% ========================================================================
%  5) QUASI-STEADY TIME LOOP
% =========================================================================

for n = 1:Nt

    %% ----- Current rigid-body state -----

    stateNow.x = bodyState.x(n);

    stateNow.z = bodyState.z(n);

    stateNow.theta = bodyState.theta(n);

    stateNow.u = bodyState.u(n);

    stateNow.w = bodyState.w(n);

    stateNow.omega = bodyState.omega(n);


    %% ----- Instantaneous aerodynamic solution -----

    results{n} = ...
        solve_sdpm_unsteady(cfg, stateNow);


    %% ----- Store selected aerodynamic quantities -----

    CL(n) = ...
        results{n}.CL_pressure;

    CM(n) = ...
        results{n}.CM;

    Gamma(n) = ...
        results{n}.Gamma;

end
%% ========================================================================
%  QUASI-STEADY ANIMATION
% =========================================================================

% Reference/body-frame collocation-point locations.
% These do NOT move with pitch/translation and are therefore preferable
% for plotting Cp vs x/c.

body0 = panel_geometry(results{1}.Bp0);

xCp_body = ...
    body0.cp(:,1) ./ cfg.airfoil.chord;


%% ===== FIGURE SETUP =====

fig = figure(10);
clf(fig)

tiledlayout(2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% ------------------------------------------------------------------------
%  1) MOVING AIRFOIL WITH Cp COLORING
% -------------------------------------------------------------------------

ax1 = nexttile;
hold(ax1,'on')
axis(ax1,'equal')
grid(ax1,'on')

airfoilLine = plot( ...
    ax1, ...
    results{1}.body.Bp(:,1), ...
    results{1}.body.Bp(:,2), ...
    '-k', ...
    'LineWidth',1.5);

cpScatter = scatter( ...
    ax1, ...
    results{1}.body.cp(:,1), ...
    results{1}.body.cp(:,2), ...
    35, ...
    results{1}.Cp, ...
    'filled');

cb = colorbar(ax1);
cb.Label.String = 'C_p';

xlabel(ax1,'x')
ylabel(ax1,'z')

title(ax1,'Airfoil Surface Pressure')

% Set limits based on expected motion.
% Adjust these if you later prescribe large translation.
xlim(ax1,[-0.25 1.25] .* cfg.airfoil.chord)
ylim(ax1,[-0.6 0.6] .* cfg.airfoil.chord)


%% ------------------------------------------------------------------------
%  2) INSTANTANEOUS Cp DISTRIBUTION
% -------------------------------------------------------------------------

ax2 = nexttile;
hold(ax2,'on')
grid(ax2,'on')

set(ax2,'YDir','reverse')

cpLine = plot( ...
    ax2, ...
    xCp_body, ...
    results{1}.Cp, ...
    '-o', ...
    'LineWidth',1.2, ...
    'MarkerSize',3);

xlabel(ax2,'x/c')
ylabel(ax2,'C_p')

title(ax2,'Instantaneous Pressure Distribution')

xlim(ax2,[0 1])


%% ------------------------------------------------------------------------
%  3) LIFT HISTORY
% -------------------------------------------------------------------------

ax3 = nexttile;
hold(ax3,'on')
grid(ax3,'on')

plot( ...
    ax3, ...
    t, ...
    CL, ...
    '-', ...
    'LineWidth',1.2);

CLmarker = plot( ...
    ax3, ...
    t(1), ...
    CL(1), ...
    'o', ...
    'MarkerFaceColor','auto');

xlabel(ax3,'Time [s]')
ylabel(ax3,'C_L')

title(ax3,'Lift History')

xlim(ax3,[t(1) t(end)])


%% ------------------------------------------------------------------------
%  4) PRESCRIBED MOTION
% -------------------------------------------------------------------------

ax4 = nexttile;
hold(ax4,'on')
grid(ax4,'on')

thetaDeg = ...
    rad2deg(bodyState.theta);

plot( ...
    ax4, ...
    t, ...
    thetaDeg, ...
    '-', ...
    'LineWidth',1.2);

thetaMarker = plot( ...
    ax4, ...
    t(1), ...
    thetaDeg(1), ...
    'o', ...
    'MarkerFaceColor','auto');

xlabel(ax4,'Time [s]')
ylabel(ax4,'\theta [deg]')

title(ax4,'Prescribed Pitch Motion')

xlim(ax4,[t(1) t(end)])


%% ===== OVERALL TITLE =====

sgtitle( ...
    sprintf( ...
    'NACA %s Quasi-Steady SDPM', ...
    cfg.airfoil.code));


%% ========================================================================
%  ANIMATION LOOP
% =========================================================================

for n = 1:Nt

    current = results{n};


    %% ----- Moving airfoil geometry -----

    airfoilLine.XData = ...
        current.body.Bp(:,1);

    airfoilLine.YData = ...
        current.body.Bp(:,2);


    %% ----- Cp-colored surface points -----

    cpScatter.XData = ...
        current.body.cp(:,1);

    cpScatter.YData = ...
        current.body.cp(:,2);

    cpScatter.CData = ...
        current.Cp;


    %% ----- Instantaneous Cp distribution -----

    cpLine.YData = ...
        current.Cp;


    %% ----- Moving CL marker -----

    CLmarker.XData = ...
        t(n);

    CLmarker.YData = ...
        CL(n);


    %% ----- Moving pitch-angle marker -----

    thetaMarker.XData = ...
        t(n);

    thetaMarker.YData = ...
        thetaDeg(n);


    %% ----- Dynamic title -----

    sgtitle( ...
        sprintf( ...
        ['NACA %s Quasi-Steady SDPM   |   ', ...
         't = %.3f s   |   \\theta = %.2f^\\circ   |   C_L = %.3f'], ...
        cfg.airfoil.code, ...
        t(n), ...
        thetaDeg(n), ...
        CL(n)));


    drawnow

end

%% ========================================================================
%  6) BASIC VALIDATION / TIME-HISTORY PLOTS
% =========================================================================

figure(1)
clf

plot( ...
    t, ...
    rad2deg(bodyState.theta), ...
    'LineWidth', 1.2)

xlabel('Time [s]')
ylabel('\theta [deg]')

title('Prescribed Pitch Motion')

grid on


figure(2)
clf

plot( ...
    t, ...
    CL, ...
    'LineWidth', 1.2)

xlabel('Time [s]')
ylabel('C_L')

title('Quasi-Steady Lift History')

grid on


figure(3)
clf

plot( ...
    rad2deg(bodyState.theta), ...
    CL, ...
    'LineWidth', 1.2)

xlabel('\theta [deg]')
ylabel('C_L')

title('Quasi-Steady Lift vs Pitch Angle')

grid on


alphaEffective = ...
    cfg.flow.alphaDeg - rad2deg(bodyState.theta);

figure
plot(t,alphaEffective)
xlabel('Time [s]')
ylabel('\alpha_{eff} [deg]')
grid on

figure
plot(rad2deg(bodyState.theta), CL, 'LineWidth', 1.2)
xlabel('\theta [deg]')
ylabel('C_L')
title(sprintf('C_l vs \\theta (\\pi f_c / U  = %3d)' , redFreq))
grid on