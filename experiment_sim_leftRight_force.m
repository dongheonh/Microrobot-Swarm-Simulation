clear; clc; close all;

%% TUNE
SCALE_EXTERNAL = 28;
SCALE_INTERNAL = 0.625 * 1e-06;
C_EFF          = 0.00028;

%% Parameters
MASS_DISK = 2.4e-7;

EXTERNAL_RANGE = 39.0;
INTERNAL_RANGE = 18.0;

DURATION     = 1.0;
N_DISKS      = 200;
spawn_radius = 30.0;
DISK_RADIUS  = 1.5;

SWARM_MARGIN_FACTOR = 1.5;
VELOCITY_HORIZON    = 7;

%% internal force

PAIR_A = 9318.0;
PAIR_N = 4.22;

%% external magnetic force

DISK_MOMENT = 1.7e-5;

BZ_TABLE = [
     357.3266,  30.3195,   7.1030,   2.6410;
    -30.1112,   11.3810,   3.9168,   2.1647;
     -6.6998,   -0.1569,   1.3462,   0.9944;
     -1.9077,   -0.7306,   0.2489,   0.2774;
];

BZ_RADIAL_PROFILE = 0.5 * (BZ_TABLE(:, 3) + BZ_TABLE(:, 4));

GAUSS_TO_TESLA = 1e-4;
MM_TO_M        = 1e-3;
RADIAL_STEP    = 10.0;

%% electromagnet arrangement

N_MAGS    = 32;
EM_RADIUS = 10.0;
EM_PITCH  = 20.0;

EM_X = (0 : N_MAGS-1)' * EM_PITCH;
EM_Y = zeros(N_MAGS, 1);

%% simulation

DT  = 0.001;
FPS = 15;

STEPS_PER_FRAME = max(1, round(1.0 / (FPS * DT)));

TOTAL_TIME   = N_MAGS * DURATION;   % no warm-up: magnet 1 on from t=0
TOTAL_FRAMES = floor(TOTAL_TIME / (DT * STEPS_PER_FRAME));

%% disk placement

rng(42);

disk_x = zeros(N_DISKS, 1);
disk_y = zeros(N_DISKS, 1);

disk_vx = zeros(N_DISKS, 1);
disk_vy = zeros(N_DISKS, 1);

attempts = 0;
placed   = 0;

MIN_DIST = 2 * DISK_RADIUS;

while placed < N_DISKS

    angle       = rand * 2 * pi;
    distance    = MIN_DIST + rand * (spawn_radius - MIN_DIST);
    candidate_x = EM_X(1) + distance * cos(angle);
    candidate_y = EM_Y(1) + distance * sin(angle);

    no_overlap = true;

    for k = 1:placed
        if hypot(candidate_x - disk_x(k), candidate_y - disk_y(k)) < MIN_DIST * 1.05
            no_overlap = false;
            break;
        end
    end

    if no_overlap
        placed         = placed + 1;
        disk_x(placed) = candidate_x;
        disk_y(placed) = candidate_y;
    end

    attempts = attempts + 1;

    if attempts > 50000
        spawn_radius = spawn_radius + 5.0;
        attempts     = 0;
    end
end

%% swarm radius — defined at t=0.5 s after magnet 1 activates

RADIUS_CAPTURE_TIME = 0.5;   % seconds after first magnet activates

t_sim = 0.0;
for s = 1:round(RADIUS_CAPTURE_TIME / DT)
    t_sim = t_sim + DT;
    active_magnet = min(floor(t_sim / DURATION) + 1, N_MAGS);
    [disk_x, disk_y, disk_vx, disk_vy] = physics_step( ...
        disk_x, disk_y, disk_vx, disk_vy, active_magnet, ...
        EM_X, EM_Y, MASS_DISK, C_EFF, DISK_MOMENT, ...
        SCALE_EXTERNAL, SCALE_INTERNAL, PAIR_A, PAIR_N, ...
        GAUSS_TO_TESLA, MM_TO_M, RADIAL_STEP, BZ_RADIAL_PROFILE, ...
        EXTERNAL_RANGE, INTERNAL_RANGE, MIN_DIST, DT);
end

initial_center_x  = mean(disk_x);
initial_center_y  = mean(disk_y);
dist_from_center  = hypot(disk_x - initial_center_x, disk_y - initial_center_y);
initial_radius_mm = max(dist_from_center + DISK_RADIUS);
swarm_radius_mm   = SWARM_MARGIN_FACTOR * initial_radius_mm;

fprintf('Swarm radius defined at t=%.1f s (magnet on):\n', RADIUS_CAPTURE_TIME);
fprintf('  raw radius      = %.3f mm\n', initial_radius_mm);
fprintf('  boundary (1.5x) = %.3f mm\n', swarm_radius_mm);

%% visualization

X_MARGIN = 100;
Y_MARGIN = 100;

X_MIN = EM_X(1)   - X_MARGIN;
X_MAX = EM_X(end) + X_MARGIN;
Y_MIN = -Y_MARGIN;
Y_MAX =  Y_MARGIN;

fig = figure('Color', 'w', 'Position', [50, 200, 1600, 400], ...
             'MenuBar', 'none', 'ToolBar', 'none');

set(fig, 'Renderer', 'opengl');

ax = axes(fig, 'Color', 'w', 'Position', [0.01 0.01 0.98 0.90]);

hold(ax, 'on');
axis(ax, 'equal', 'off');
xlim(ax, [X_MIN, X_MAX]);
ylim(ax, [Y_MIN, Y_MAX]);

theta_em = linspace(0, 2*pi, 60);

for i = 1:N_MAGS
    plot(ax, ...
        EM_X(i) + EM_RADIUS * cos(theta_em), ...
        EM_Y(i) + EM_RADIUS * sin(theta_em), ...
        'Color', [0.8 0.8 0.8], ...
        'LineWidth', 0.5);
end

theta_disk = linspace(0, 2*pi, 40)';
unit_circle_x = cos(theta_disk);
unit_circle_y = sin(theta_disk);

h_disks = gobjects(N_DISKS, 1);

for i = 1:N_DISKS
    h_disks(i) = patch(ax, ...
        disk_x(i) + DISK_RADIUS * unit_circle_x, ...
        disk_y(i) + DISK_RADIUS * unit_circle_y, ...
        'k', 'EdgeColor', 'k', 'FaceColor', 'k', 'LineWidth', 0.7);
end

theta_swarm = linspace(0, 2*pi, 200);

h_swarm_boundary = plot(ax, nan, nan, ...
    'Color', [0.0 0.7 0.2], ...
    'LineWidth', 2.0);

h_swarm_center = plot(ax, nan, nan, ...
    'o', ...
    'MarkerFaceColor', [0.0 0.7 0.2], ...
    'MarkerEdgeColor', [0.0 0.7 0.2], ...
    'MarkerSize', 6);

h_right_anchor = plot(ax, nan, nan, ...
    'o', ...
    'MarkerFaceColor', [1.0 0.0 0.0], ...
    'MarkerEdgeColor', [1.0 0.0 0.0], ...
    'MarkerSize', 5);

h_text = text(ax, X_MIN + 5, Y_MAX - 10, '', ...
    'FontSize', 11, ...
    'Color', 'k', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left');

drawnow;

% script_dir = fileparts(mfilename('fullpath'));
% output_path = fullfile(script_dir, 'microdisk_transport_newton_model.mp4');
% 
% vid = VideoWriter(output_path, 'MPEG-4');
% vid.FrameRate = FPS;
% vid.Quality   = 90;
% open(vid);

%% swarm evaluation variables

prev_anchor_x = nan;
prev_anchor_y = nan;
travel_anchor_mm = 0.0;

prev_frame_x = disk_x;
prev_frame_y = disk_y;

agent_speed_hist = nan(N_DISKS, VELOCITY_HORIZON);
speed_hist_col   = 1;

rec_t          = nan(TOTAL_FRAMES, 1);
rec_dlb        = nan(TOTAL_FRAMES, 1);
rec_left_frac  = nan(TOTAL_FRAMES, 1);
rec_anchor_spd = nan(TOTAL_FRAMES, 1);

%% animation loop
% t_sim already advanced by one frame (radius capture step above)

for frame = 1:TOTAL_FRAMES

    for s = 1:STEPS_PER_FRAME

        t_sim = t_sim + DT;

        % magnet 1 active from t=0 (no warm-up)
        active_magnet = min(floor(t_sim / DURATION) + 1, N_MAGS);

        [disk_x, disk_y, disk_vx, disk_vy] = physics_step( ...
            disk_x, disk_y, disk_vx, disk_vy, active_magnet, ...
            EM_X, EM_Y, ...
            MASS_DISK, C_EFF, DISK_MOMENT, ...
            SCALE_EXTERNAL, SCALE_INTERNAL, ...
            PAIR_A, PAIR_N, ...
            GAUSS_TO_TESLA, MM_TO_M, RADIAL_STEP, BZ_RADIAL_PROFILE, ...
            EXTERNAL_RANGE, INTERNAL_RANGE, ...
            MIN_DIST, DT);
    end

    %% evaluate current swarm using rightmost-anchor circle

    frame_dt = STEPS_PER_FRAME * DT;
    rec_t(frame) = t_sim;
    speed_now = hypot(disk_x - prev_frame_x, disk_y - prev_frame_y) / frame_dt;

    agent_speed_hist(:, speed_hist_col) = speed_now;
    speed_hist_col = speed_hist_col + 1;

    if speed_hist_col > VELOCITY_HORIZON
        speed_hist_col = 1;
    end

    prev_frame_x = disk_x;
    prev_frame_y = disk_y;

    [~, right_idx] = max(disk_x);

        right_x = disk_x(right_idx);
        right_y = disk_y(right_idx);

        anchor_center_x = right_x - swarm_radius_mm;
        anchor_center_y = right_y;

        anchor_frame_speed = nan;
        if ~isnan(prev_anchor_x)
            d_anchor = hypot(anchor_center_x - prev_anchor_x, anchor_center_y - prev_anchor_y);
            travel_anchor_mm   = travel_anchor_mm + d_anchor;
            anchor_frame_speed = d_anchor / frame_dt;
        end

        prev_anchor_x = anchor_center_x;
        prev_anchor_y = anchor_center_y;

        dist_current = hypot(disk_x - anchor_center_x, disk_y - anchor_center_y);
        inside_swarm = dist_current <= swarm_radius_mm;
        outside_count = sum(~inside_swarm);

        rec_dlb(frame)        = outside_count;
        rec_left_frac(frame)  = outside_count / N_DISKS;
        rec_anchor_spd(frame) = anchor_frame_speed;

        rolling_speed = zeros(N_DISKS, 1);
        has_speed = false(N_DISKS, 1);

        for i = 1:N_DISKS
            vals = agent_speed_hist(i, :);
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                rolling_speed(i) = mean(vals);
                has_speed(i) = true;
            end
        end

        valid_inside = inside_swarm & has_speed;

        if any(valid_inside)
            vel_mean = mean(rolling_speed(valid_inside));
            vel_max  = max(rolling_speed(valid_inside));
            vel_min  = min(rolling_speed(valid_inside));
        else
            vel_mean = nan;
            vel_max  = nan;
            vel_min  = nan;
        end

        set(h_swarm_boundary, ...
            'XData', anchor_center_x + swarm_radius_mm * cos(theta_swarm), ...
            'YData', anchor_center_y + swarm_radius_mm * sin(theta_swarm));

        set(h_swarm_center, ...
            'XData', anchor_center_x, ...
            'YData', anchor_center_y);

        set(h_right_anchor, ...
            'XData', right_x, ...
            'YData', right_y);

        text_str = sprintf([ ...
            'Initial swarm radius(mm): %.2f\n', ...
            'Boundary radius after 1.5x(mm): %.2f\n', ...
            'Anchor-center travel distance(mm): %.2f\n', ...
            'Robots outside boundary: %d\n', ...
            'Swarm agent velocity(mm/s): mean=%.2f max=%.2f min=%.2f\n', ...
            'Active magnet: %d'], ...
            initial_radius_mm, ...
            swarm_radius_mm, ...
            travel_anchor_mm, ...
            outside_count, ...
            vel_mean, vel_max, vel_min, ...
            active_magnet);


    set(h_text, 'String', text_str);

    %% update disks

    for i = 1:N_DISKS
        set(h_disks(i), ...
            'XData', disk_x(i) + DISK_RADIUS * unit_circle_x, ...
            'YData', disk_y(i) + DISK_RADIUS * unit_circle_y);
    end

    drawnow;
    pause(0.001);

    % frame_img = getframe(fig);
    % writeVideo(vid, frame_img);

    if mod(frame, 30) == 0
        fprintf('\r  %5.1f%%  frame %d / %d  active magnet = %d', ...
                100 * frame / TOTAL_FRAMES, frame, TOTAL_FRAMES, active_magnet);
    end
end

% close(vid);

%% Save simulation metrics to CSV

simDir = fileparts(mfilename('fullpath'));

valid_frames = ~isnan(rec_dlb);
t_out    = rec_t(valid_frames);
dlb_out  = rec_dlb(valid_frames);
frac_out = rec_left_frac(valid_frames);
spd_out  = rec_anchor_spd(valid_frames);

T1 = table(initial_radius_mm, swarm_radius_mm, ...
    'VariableNames', {'initial_radius_mm', 'check_radius_mm'});
writetable(T1, fullfile(simDir, 'sim_swarm_size.csv'));
fprintf('\nSaved sim_swarm_size.csv\n');

T2 = table(t_out, dlb_out, frac_out, ...
    'VariableNames', {'t', 'disks_left_behind', 'left_fraction'});
writetable(T2, fullfile(simDir, 'sim_left_behind_time.csv'));
fprintf('Saved sim_left_behind_time.csv  (%d rows)\n', height(T2));

T3 = table(dlb_out(end), frac_out(end), ...
    'VariableNames', {'final_disks_left_behind', 'final_left_fraction'});
writetable(T3, fullfile(simDir, 'sim_left_behind_final.csv'));
fprintf('Saved sim_left_behind_final.csv\n');

T4 = table(t_out, spd_out, ...
    'VariableNames', {'t', 'anchor_speed_mm_s'});
writetable(T4, fullfile(simDir, 'sim_anchor_speed_time.csv'));
fprintf('Saved sim_anchor_speed_time.csv  (%d rows)\n', height(T4));

%% Functions

function [new_x, new_y, new_vx, new_vy] = physics_step( ...
        disk_x, disk_y, disk_vx, disk_vy, active_magnet, ...
        EM_X, EM_Y, ...
        MASS_DISK, C_EFF, DISK_MOMENT, ...
        SCALE_EXTERNAL, SCALE_INTERNAL, ...
        PAIR_A, PAIR_N, ...
        GAUSS_TO_TESLA, MM_TO_M, RADIAL_STEP, BZ_RADIAL_PROFILE, ...
        EXTERNAL_RANGE, INTERNAL_RANGE, ...
        MIN_DIST, DT)

    if active_magnet == 0
        Fx_mag = zeros(size(disk_x));
        Fy_mag = zeros(size(disk_y));
    else
        [Fx_mag, Fy_mag] = calc_magnetic_force( ...
            disk_x, disk_y, EM_X(active_magnet), EM_Y(active_magnet), ...
            DISK_MOMENT, GAUSS_TO_TESLA, MM_TO_M, ...
            RADIAL_STEP, BZ_RADIAL_PROFILE, EXTERNAL_RANGE);
    end

    [Fx_pair, Fy_pair] = calc_pair_force( ...
        disk_x, disk_y, PAIR_A, PAIR_N, MIN_DIST, INTERNAL_RANGE);

    Fx_damp = -C_EFF * disk_vx;
    Fy_damp = -C_EFF * disk_vy;

    Fx_total = SCALE_EXTERNAL * Fx_mag + SCALE_INTERNAL * Fx_pair + Fx_damp;
    Fy_total = SCALE_EXTERNAL * Fy_mag + SCALE_INTERNAL * Fy_pair + Fy_damp;

    ax = Fx_total / MASS_DISK;
    ay = Fy_total / MASS_DISK;

    new_vx = disk_vx + ax * DT;
    new_vy = disk_vy + ay * DT;

    new_x = disk_x + new_vx * DT / MM_TO_M;
    new_y = disk_y + new_vy * DT / MM_TO_M;
end

function [Fx, Fy] = calc_magnetic_force( ...
        disk_x, disk_y, magnet_x, magnet_y, ...
        DISK_MOMENT, GAUSS_TO_TESLA, MM_TO_M, ...
        RADIAL_STEP, BZ_RADIAL_PROFILE, EXTERNAL_RANGE)

    dx = disk_x - magnet_x;
    dy = disk_y - magnet_y;

    r = hypot(dx, dy);

    in_range = double(r < EXTERNAL_RANGE);

    [~, dbz_dr_G_per_mm] = lookup_bz(r, BZ_RADIAL_PROFILE, RADIAL_STEP);

    r_safe = max(r, 1e-9);

    dbz_dr_T_per_m = dbz_dr_G_per_mm * GAUSS_TO_TESLA / MM_TO_M;

    grad_bz_x = dbz_dr_T_per_m .* (dx ./ r_safe) .* in_range;
    grad_bz_y = dbz_dr_T_per_m .* (dy ./ r_safe) .* in_range;

    Fx = DISK_MOMENT * grad_bz_x;
    Fy = DISK_MOMENT * grad_bz_y;
end

function [Fx_pair, Fy_pair] = calc_pair_force( ...
        disk_x, disk_y, PAIR_A, PAIR_N, MIN_DIST, INTERNAL_RANGE)

    n_disks = length(disk_x);

    Fx_pair = zeros(n_disks, 1);
    Fy_pair = zeros(n_disks, 1);

    for i = 1:n_disks-1

        dx = disk_x(i) - disk_x(i+1:end);
        dy = disk_y(i) - disk_y(i+1:end);

        r_raw = hypot(dx, dy);
        active = double(r_raw < INTERNAL_RANGE);

        r = max(r_raw, MIN_DIST);

        F_pair = active .* PAIR_A .* r .^ (-PAIR_N);

        ux = dx ./ r;
        uy = dy ./ r;

        Fx_pair(i)       = Fx_pair(i)       + sum(F_pair .* ux);
        Fy_pair(i)       = Fy_pair(i)       + sum(F_pair .* uy);

        Fx_pair(i+1:end) = Fx_pair(i+1:end) - F_pair .* ux;
        Fy_pair(i+1:end) = Fy_pair(i+1:end) - F_pair .* uy;
    end
end

function [bz_value, dbz_dr] = lookup_bz(r_array, BZ_RADIAL_PROFILE, RADIAL_STEP)

    idx_float = min(max(r_array / RADIAL_STEP, 0.0), 2.9999);

    idx_low  = floor(idx_float) + 1;
    idx_high = idx_low + 1;
    frac     = idx_float - floor(idx_float);

    bz_low  = BZ_RADIAL_PROFILE(idx_low);
    bz_high = BZ_RADIAL_PROFILE(idx_high);

    bz_value = bz_low .* (1 - frac) + frac .* bz_high;

    dbz_dr = (bz_high - bz_low) / RADIAL_STEP;
end