function wake = shed_wake(wakeOld,newWakePanel, muNew)
%SHED_WAKE 
%
% Appends the newly solved TE wake panel to the persistent wake state
% history
%
% wake ordering convention:
%           - wake.Bp(1,:) = youngest wake node
%
%           - wake.mu(1) = newest wake panel strength
%
% Older wake panels follow downstream
%
%
% ===== INPUTS =====
% ---> wakeOld = wake state that enters the current timestep
% ---> newWakePanel = Geometry of the newly created TE wake panel
% ---> muNew = solved doublet strength of new wake panel
%
%
% ===== OUTPUTS =====
% ---> wake = updated persistent wake which will be input to the next
% timestep


%% ========== 1) CHECK NEW PANEL ==========

if newWakePanel.nPanels ~= 1
    error('newWakePanel must contain exactly one panel.')
end

if ~isscalar(muNew)
    error('muNew must be a scalar.')
end

%% ========== 2) FIRST WAKE PANEL ==========

if isempty(wakeOld.Bp)

    % no previous wake exists - boundary points of newWakePanel become
    % initial wake geometry

    wake.Bp = newWakePanel.Bp;
    wake.mu = muNew;

    return

end

%% ========== 3) VERIFY CONNECTION TO PERSISTENT WAKE ===========

% downstream point of newest wake panel must be coincident with the
% youngest node of the existing wake

newDownstreamPoint = newWakePanel.Bp(end,:);

oldYoungestPoint = wakeOld.Bp(1,:);

connectionError = norm(newDownstreamPoint - oldYoungestPoint);

if connectionError > 1e-10

    error(['New wake panel does not connect to the existing youngest wake node.']);

end

%% ========== 4) PREPEND NEW WAKE NODE ==========

% newWakePanel.Bp is =  [ current TE     
%                           old youngest wake node]

% second row already exists as wakeOld.Bp(1,:). so only NEW current TE node
% needs to be prepended

wake.Bp = [...
    newWakePanel.Bp(1,:);
    wakeOld.Bp];

%% ========== 5) PREPEND NEW WAKE STRENGTH ==========

wake.mu = [...
    muNew;
    wakeOld.mu(:)];

%% ========== 6) CONSISTENCY CHECK ==========

nWakePanels = size(wake.Bp,1) - 1;

if numel(wake.mu) ~= nWakePanels

    error('Wake geometry and strength arrays are inconsistent.');
end

end

