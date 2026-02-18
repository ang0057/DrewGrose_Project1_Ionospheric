clear;
close all;
clc;

% ==========================================
% setup general stuff

% uncomment the desired date
% My birthday: August 15, 2003. Using 7am for the time.
UT = [2003 8 15 7 0];        % UT - year, month, day, hour, minute
% THE NEXT DAY 24HR LATER:
% UT = [2003 8 16 7 0];        % UT - year, month, day, hour, minute
R12 = 100;                   % R12 index
speed_of_light = 2.99792458e8;

% uncomment the desired transmitter location
transmit = "VA";
% transmit = "TX";

% switch statement to make it easier for me to change transmitters
switch transmit
    case "VA"
        % Chesapeake, VA
        origin_lat = 36.77;          % latitude of the start point of ray
        origin_long = -76.29;         % longitude of the start point of ray
        fprintf('Transmitter Location: Chesapeake, VA\n');
    case "TX"
        % Corpus Christi, TX
        origin_lat = 27.80;          % latitude of the start point of ray
        origin_long = -97.40;         % longitude of the start point of ray
        fprintf('Transmitter Location: Corpus Christi, TX\n');
end

% Auburn, AL
receiver_lat = 32.61;          % latitude of the receiver
receiver_long = -85.48;         % longitude of the receiver

% Set Origin Height (Required for 3D raytrace)
origin_ht = 0;

% calculate bearing and range
[dist_m, ray_bear] = latlon2raz(receiver_lat,receiver_long,origin_lat,origin_long);
dist_km = dist_m/1000;
fprintf('Target Range: %.2f km | Bearing: %.2f deg\n', dist_km, ray_bear);

% these are leftover from ray_test4.m
% DOPPLER FLAG SET TO 0
doppler_flag = 0;            % generate ionosphere 5 minutes later so that
                             % Doppler shift can be calculated
irregs_flag = 0;             % no irregularities - not interested in 
                             % Doppler spread or field aligned irregularities
kp = 0;                      % kp not used as irregs_flag = 0. Set it to a 
                             % dummy value 

%
% generate ionospheric, geomagnetic and irregularity grids
%
max_range = dist_km + 1000;      % maximum range for sampling the ionosphere (km)
num_range = 201;        % number of ranges (must be < 2000)
range_inc = max_range ./ (num_range - 1);  % range cell size (km)

start_height = 0 ;      % start height for ionospheric grid (km)
height_inc = 3;         % height increment (km)
num_heights = 200;      % number of  heights (must be < 2000)

% ==========================================
% Generate 3D Grid

% 3D param setup

% Define Grid Boundaries (Must cover Tx and Rx with buffer)
% Path: Chesapeake (~36.8N, -76.3E) <-> Auburn (~32.6N, -85.5E)

lat_start = 25.0;      % Start south of both points
lat_step  = 0.5;       % 0.5 degree resolution
num_lats  = 40;        % Covers 25.0 to 45.0 degrees

lon_start = -100.0;    % Start west of both (covers TX and AL)
lon_step  = 1.0;       % 1.0 degree resolution
num_lons  = 40;        % Covers -100.0 to -60.0 degrees

ht_start  = 0;         % Start height (km)
ht_step   = 3;         % Height step (km)
num_hts   = 200;       % Covers 0 to 600 km

% (Required by PHaRLAP 3D)
% Vector Format: [lat_start, lat_step, num_lats, lon_start, lon_step, num_lons, h_start, h_step, num_hts]

iono_grid_parms = [lat_start; lat_step; num_lats; ...
                   lon_start; lon_step; num_lons; ...
                   ht_start; ht_step; num_hts];

geomag_grid_parms = iono_grid_parms; % Usually safe to use same grid for B-field


% actually generate it
fprintf('Generating 3D Ionosphere Grid...\n');
tic
[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d(UT, R12, iono_grid_parms, geomag_grid_parms, ...
                     doppler_flag);
toc

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

% freq and elevation angle setup
freqs = 1:0.1:10; % 1-10 MHz in 0.1 MHz increments
elevs = 3:0.2:85; % dense fan of elevation angles to hit the receiver
num_elevs = length(elevs);
tol = 1e-7; % ODE tolerance
nhops = 2;

% 3D requirement
ray_bear_vec = ones(size(elevs)) * ray_bear;

% ==========================================
% MAIN SIMULATION LOOP

% initialize results storage for all modes now
ionogram_data_2D = []; % No field (Mode 0)
ionogram_data_O  = []; % O-Mode (Mode 1)
ionogram_data_X  = []; % X-Mode (Mode -1

fprintf('Starting Frequency Sweep (%.1f - %.1f MHz)...\n', 1, 10);

for f = freqs

    % run raytrace_3d for the whole fan of angles
    nhops = 1;
    freq_vec = ones(size(elevs)) * f;

    % Loop through the 3 modes: 0 (No B-Field), 1 (O-Mode), -1 (X-Mode)
    % (raytrace_3d has the mode input formatted using these numbers)
    for mode = [0, 1, -1]

        % Run 3D Raytrace
        [ray_data, ~] = raytrace_3d(origin_lat, origin_long, origin_ht, ...
            elevs, ray_bear_vec, freq_vec, mode, nhops, tol, ...
            iono_en_grid, iono_en_grid_5, collision_freq, iono_grid_parms, ...
            Bx, By, Bz, geomag_grid_parms);

        % Organize data
        g_ranges = [ray_data.ground_range];
        p_ranges = [ray_data.group_range]; % Virtual path length (P')

        % Identify skip distance
        [~, min_idx] = min(g_ranges);

        % Isolate Low Ray
        g_subset = g_ranges(1:min_idx);
        p_subset = p_ranges(1:min_idx);
        elev_subset = elevs(1:min_idx);

        % Interpolate
        if length(g_subset) >= 2
            [g_sorted, sort_idx] = sort(g_subset);
            p_sorted = p_subset(sort_idx);
            elev_sorted = elev_subset(sort_idx);

            min_sim_range = min(g_sorted);
            max_sim_range = max(g_sorted);

            if dist_km >= min_sim_range && dist_km <= max_sim_range
                P_prime = interp1(g_sorted, p_sorted, dist_km, 'linear');
                target_elev = interp1(g_sorted, elev_sorted, dist_km, 'linear');
                h_virtual = sqrt( (P_prime/2)^2 - (dist_km/2)^2 );

                % Save result into the correct array based on the current mode
                if mode == 0
                    ionogram_data_2D = [ionogram_data_2D; f, h_virtual, target_elev];
                elseif mode == 1
                    ionogram_data_O = [ionogram_data_O; f, h_virtual, target_elev];
                elseif mode == -1
                    ionogram_data_X = [ionogram_data_X; f, h_virtual, target_elev];
                end
            end
        end
    end
end

% ==========================================
% Plotting
% ionogram
figure('Color', 'w');  % forcing to white bc I'm in dark mode
hold on;

% Plot each mode if data exists
if ~isempty(ionogram_data_2D)
    plot(ionogram_data_2D(:,1), ionogram_data_2D(:,2), 'k.-', 'LineWidth', 1.5, 'DisplayName', '2D Eqv. (No B-Field)');
end
if ~isempty(ionogram_data_O)
    plot(ionogram_data_O(:,1), ionogram_data_O(:,2), 'b.-', 'LineWidth', 1.5, 'DisplayName', '3D O-Mode');
end
if ~isempty(ionogram_data_X)
    plot(ionogram_data_X(:,1), ionogram_data_X(:,2), 'r.-', 'LineWidth', 1.5, 'DisplayName', '3D X-Mode');
end

grid on;
legend('Location', 'best');
xlabel('Frequency (MHz)');
ylabel('Virtual Reflection Height (km)');

% switch statement to easily update title based on date
switch UT(3)
    case 15
        % my birthday
        title(['8/15/2003 7:00am, 2D vs 3D Ionogram: ' num2str(dist_km, '%.1f') ' km path'], 'Color','k');
    case 16
        % the next day
        title(['8/16/2003 7:00am, 2D vs 3D Ionogram: ' num2str(dist_km, '%.1f') ' km path'], 'Color','k');
end
% switch statement to easily update subtitle based on transmitter
switch transmit
    case "VA"
        subtitle('Tx: Chesapeake -> Rx: Auburn', 'Color','k');
    case "TX"
        subtitle('Tx: Corpus Christi -> Rx: Auburn', 'Color','k');
end

xlim([1 10]);
ylim([0 650]); % Ionogram height range from instructions
set(gca, 'Color', 'w');  % Force axes background to white
set(gca, 'XColor', 'k', 'YColor', 'k'); % Force axis lines to black
