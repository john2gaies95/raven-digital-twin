%% Project RAVEN — 3DOF Point-Mass Boost+Coast Trajectory (with 3D Animation)
% Clean baseline you can put on GitHub and iterate:
% - Vertical-plane 3DOF point-mass (x downrange, z altitude)
% - Thrust during burn, coast after burnout
% - Drag with simple Cd(Mach) and exponential density
% - Variable gravity with altitude
% - Simple pitch program (hold near-vertical then pitch over)
% - Mach computed with ISA speed of sound (non-negative)
% - 3D animation (maps 2D trajectory into 3D with y=0)
%
% Run:
%   >> traj_3dof_pointmass
%
% Outputs:
%   - Console: Max Mach, Max Altitude, Downrange, Burn time, Final time
%   - Figures: Mach vs Time, Altitude vs Downrange
%   - Animation: 3D trajectory playback (optional)

clear; clc; close all;

%% -----------------------------
% 0) Output options
% ------------------------------
MAKE_ANIMATION = true;     % set false if you don't want animation
SAVE_FRAMES    = false;    % if true, saves PNG frames to Results/animation_frames
FRAME_STRIDE   = 3;        % animate every Nth sample (bigger = faster)

% Create output folders relative to this file
repoRoot   = fileparts(mfilename('fullpath'));
resultsDir = fullfile(repoRoot, '..', 'Results');
figDir     = fullfile(resultsDir, 'figures');
animDir    = fullfile(resultsDir, 'animation_frames');

if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
if ~exist(figDir,'dir'),     mkdir(figDir);     end
if SAVE_FRAMES && ~exist(animDir,'dir'), mkdir(animDir); end

%% -----------------------------
% 1) Environment / Atmosphere
% ------------------------------
env.Re   = 6371e3;       % m Earth radius
env.g0   = 9.80665;      % m/s^2
env.rho0 = 1.225;        % kg/m^3 sea-level density
env.H    = 7200;         % m scale height (simple exponential)
env.gamma = 1.4;         % air
env.R     = 287.05;      % J/(kg K)

%% -----------------------------
% 2) Vehicle / Aero / Propulsion
% ------------------------------
veh.m0   = 250;          % kg initial mass
veh.mp   = 120;          % kg propellant mass
veh.tb   = 35;           % s burn time
veh.T    = 18000;        % N thrust (constant during burn)
veh.Aref = 0.06;         % m^2 reference area
veh.Cd0  = 0.25;         % base Cd (subsonic)
veh.Cd1  = 0.40;         % supersonic Cd (rough)
veh.mdot = veh.mp / veh.tb; % kg/s

%% -----------------------------
% 3) Initial conditions + Sim time
% ------------------------------
x0    = 0;          % m downrange
z0    = 0;          % m altitude
v0    = 0;          % m/s
gamma0 = deg2rad(90); % rad (flight-path angle from +x axis; 90 = vertical up)

tFinal = 300;       % s
dt     = 0.05;      % s (smaller = more stable)
N      = floor(tFinal/dt) + 1;

t     = (0:N-1)' * dt;
x     = zeros(N,1);
z     = zeros(N,1);
v     = zeros(N,1);
gam   = zeros(N,1);
m     = zeros(N,1);
Mach  = zeros(N,1);

% init
x(1)=x0; z(1)=z0; v(1)=v0; gam(1)=gamma0; m(1)=veh.m0;

%% -----------------------------
% 4) Pitch program (simple and stable)
% ------------------------------
% Hold near-vertical early, then ramp to a shallower angle by burnout.
pitch.t_hold      = 2.0;            % s hold initial gamma
pitch.gamma_init  = gamma0;         % rad
pitch.gamma_final = deg2rad(25);    % rad by end of burn (tune this)
pitch.t_ramp_end  = veh.tb;         % s

gamma_cmd = @(tk) pitch_profile(tk, pitch);

%% -----------------------------
% 5) Integrate EOM (Euler-Cromer)
% ------------------------------
impactIndex = N;

for k = 1:N-1
    tk = t(k);

    % Stop if vehicle hits ground
    if z(k) < 0
        impactIndex = k;
        break;
    end

    % Mass + thrust schedule
    if tk <= veh.tb && m(k) > (veh.m0 - veh.mp)
        T = veh.T;
        mdot = veh.mdot;
    else
        T = 0;
        mdot = 0;
    end

    % Atmosphere
    rho = env.rho0 * exp(-max(z(k),0)/env.H);

    % Gravity (variable with altitude)
    g = env.g0 * (env.Re/(env.Re + max(z(k),0)))^2;

    % Speed of sound (ISA-ish)
    a = isa_speed_of_sound(max(z(k),0), env);

    % Mach (force non-negative)
    Mach(k) = max(v(k)/max(a,1e-6), 0);

    % Cd (very simple Mach-dependent)
    Cd = Cd_of_Mach(Mach(k), veh);

    % Drag
    D = 0.5 * rho * v(k)^2 * Cd * veh.Aref;

    % Commanded gamma (pitch program) only during powered flight
    if tk <= veh.tb
        gam_cmd = gamma_cmd(tk);
        % Simple first-order "tracking" of commanded gamma
        tau_gam = 0.8; % s time constant (tune)
        gam_dot = (gam_cmd - gam(k))/tau_gam;
    else
        % Ballistic flight-path dynamics (no lift)
        % gamma_dot = -(g - v^2/(R+z)) * cos(gamma) / v  (approx)
        % Use a standard point-mass no-lift planar form:
        if v(k) > 1e-3
            gam_dot = -(g / v(k)) * cos(gam(k));
        else
            gam_dot = 0;
        end
    end

    % Along-track acceleration
    if m(k) > 1e-6
        v_dot = (T - D)/m(k) - g*sin(gam(k));
    else
        v_dot = -g*sin(gam(k));
    end

    % Update states (Euler-Cromer)
    v(k+1)   = max(v(k) + v_dot*dt, 0);          % do not go negative speed
    gam(k+1) = gam(k) + gam_dot*dt;

    x(k+1)   = x(k) + v(k+1)*cos(gam(k+1))*dt;
    z(k+1)   = z(k) + v(k+1)*sin(gam(k+1))*dt;

    m(k+1)   = max(m(k) - mdot*dt, veh.m0 - veh.mp);
end

% trim to impactIndex (or full length)
t    = t(1:impactIndex);
x    = x(1:impactIndex);
z    = z(1:impactIndex);
v    = v(1:impactIndex);
gam  = gam(1:impactIndex);
m    = m(1:impactIndex);
Mach = Mach(1:impactIndex);

% Compute Mach for last point
a_last = isa_speed_of_sound(max(z(end),0), env);
Mach(end) = max(v(end)/max(a_last,1e-6), 0);

%% -----------------------------
% 6) Metrics
% ------------------------------
maxMach     = max(Mach);
maxAlt_km   = max(z)/1000;
downrange_km= max(x)/1000;

fprintf("\n=== RAVEN BASELINE RESULTS ===\n");
fprintf("Max Mach:     %6.2f\n", maxMach);
fprintf("Max Altitude: %6.1f km\n", maxAlt_km);
fprintf("Downrange:    %6.1f km\n", downrange_km);
fprintf("Burn time:    %6.1f s\n", veh.tb);
fprintf("Final time:   %6.1f s\n\n", t(end));

%% -----------------------------
% 7) Plots
% ------------------------------
% Mach vs time
fig1 = figure('Name','Mach vs Time');
plot(t, Mach, 'LineWidth', 2);
grid on;
xlabel('Time (s)'); ylabel('Mach');
title('Mach vs Time');
saveas(fig1, fullfile(figDir, 'mach_vs_time.png'));

% Altitude vs downrange
fig2 = figure('Name','Trajectory (Altitude vs Downrange)');
plot(x/1000, z/1000, 'LineWidth', 2);
grid on;
xlabel('Downrange (km)'); ylabel('Altitude (km)');
title('Trajectory');
saveas(fig2, fullfile(figDir, 'altitude_vs_downrange.png'));

%% -----------------------------
% 8) 3D Animation
% ------------------------------
if MAKE_ANIMATION
    % Map into 3D with y=0 (later you can add heading / crossrange)
    X = x/1000;
    Y = zeros(size(X));
    Z = z/1000;

    fig3 = figure('Name','3D Trajectory Animation');
    ax = axes(fig3);
    grid(ax,'on');
    xlabel(ax,'Downrange X (km)');
    ylabel(ax,'Crossrange Y (km)');
    zlabel(ax,'Altitude Z (km)');
    title(ax,'RAVEN 3D Trajectory (Playback)');
    view(ax, 35, 20);
    axis(ax, 'equal');

    % Set nice bounds
    pad = 0.05;
    xlim(ax, [min(X) max(X)*(1+pad)]);
    ylim(ax, [-1 1]); % since Y=0 line; widen later if you add crossrange
    zlim(ax, [0 max(Z)*(1+pad)]);

    hold(ax,'on');
    track = plot3(ax, X(1), Y(1), Z(1), 'LineWidth', 2);
    marker = plot3(ax, X(1), Y(1), Z(1), 'o', 'MarkerSize', 8, 'MarkerFaceColor','w');

    % Optional: ground plane
    gpX = [min(X) max(X) max(X) min(X)];
    gpY = [-1 -1  1  1];
    gpZ = [0 0 0 0];
    patch(ax, gpX, gpY, gpZ, 'w', 'FaceAlpha', 0.05, 'EdgeColor', 'none');

    for k = 1:FRAME_STRIDE:length(X)
        set(track,  'XData', X(1:k), 'YData', Y(1:k), 'ZData', Z(1:k));
        set(marker, 'XData', X(k),   'YData', Y(k),   'ZData', Z(k));

        drawnow;

        if SAVE_FRAMES
            frameName = fullfile(animDir, sprintf('frame_%05d.png', k));
            exportgraphics(fig3, frameName);
        end
    end
end

%% =============================
% Local functions
% =============================

function gc = pitch_profile(tk, pitch)
% Piecewise pitch program: hold, then ramp to final by pitch.t_ramp_end
    if tk <= pitch.t_hold
        gc = pitch.gamma_init;
    else
        t0 = pitch.t_hold;
        t1 = max(pitch.t_ramp_end, t0 + 1e-6);
        s  = min(max((tk - t0)/(t1 - t0), 0), 1);
        % Smoothstep ramp (less jerk than linear)
        s  = s*s*(3 - 2*s);
        gc = (1 - s)*pitch.gamma_init + s*pitch.gamma_final;
    end
end

function Cd = Cd_of_Mach(M, veh)
% Extremely simple Cd model you can refine later
% - subsonic: veh.Cd0
% - transonic bump around M~1
% - supersonic: veh.Cd1
    if M < 0.8
        Cd = veh.Cd0;
    elseif M < 1.2
        % transonic bump
        bump = 0.15 * exp(-((M - 1.0)/0.18)^2);
        Cd = veh.Cd0 + bump;
    else
        Cd = veh.Cd1;
    end
end

function a = isa_speed_of_sound(h, env)
% Quick ISA-ish speed of sound
% Troposphere lapse to 11 km, then isothermal to 20 km (enough for baseline)
    if h < 11000
        T = 288.15 - 0.0065*h;
    else
        T = 216.65; % isothermal layer
    end
    a = sqrt(env.gamma * env.R * T);
end


