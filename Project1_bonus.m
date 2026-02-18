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

% Convert Plasma Frequency to Electron Density
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;

%%
% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

% freq and elevation angle setup
freqs = 1:0.1:10; % 1-10 MHz in 0.1 MHz increments
elevs = 3:0.2:85; % dense fan of elevation angles to hit the receiver
num_elevs = length(elevs);
tol = 1e-7; % ODE tolerance
nhops = 2;


% ==========================================
% MAIN SIMULATION LOOP

% initialize results storage: [frequency, virtual height, elev angle] eventually
ionogram_data = [];

fprintf('Starting Frequency Sweep (%.1f - %.1f MHz)...\n', 1, 10);

for f = freqs

    % run raytrace_2d for the whole fan of angles
    nhops = 1;
    freq_vec = ones(size(elevs)) * f;

    % Run 3D Raytrace
    [ray_data_3d, ray_path_data_3d] = raytrace_3d(origin_lat, origin_long, ...
        elevs, ray_bear, freq_vec, ...
        nhops, tol, iono_en_grid, ...
        B_field_grid);

    % Organize data
    g_ranges = [ray_data.ground_range];
    p_ranges = [ray_data.group_range]; % Virtual path length (P')

    % Identify the skip distance (Minimum Ground Range)
    % The Low Ray is on the left side of this minimum (lower elev angles)
    % The High Ray is on the right side (higher elev angles)
    [min_val, min_idx] = min(g_ranges);

    % Isolate the Low Ray (and associated virtual heights and elev angles)
    g_subset = g_ranges(1:min_idx);
    p_subset = p_ranges(1:min_idx);
    elev_subset = elevs(1:min_idx);

    % Interpolate (in case rays missed the receiver by landing around it)
    if length(g_subset) >= 2
        
        % Sort by ground range (required for interp1 function)
        [g_sorted, sort_idx] = sort(g_subset);
        p_sorted = p_subset(sort_idx);
        elev_sorted = elev_subset(sort_idx);
        
        % Check bounds
        min_sim_range = min(g_sorted);
        max_sim_range = max(g_sorted);
        
        % Interpolation logic
        if dist_km >= min_sim_range && dist_km <= max_sim_range
            
            % Interpolate values at the receiver distance
            P_prime = interp1(g_sorted, p_sorted, dist_km, 'linear');
            target_elev = interp1(g_sorted, elev_sorted, dist_km, 'linear');
            
            % Calculate virtual height from time-of-flight
            % (Breit-Tuve Theorem)
            h_virtual = sqrt( (P_prime/2)^2 - (dist_km/2)^2 );
            
            % Save result: [frequency, virtual height, elev angle]
            % also storing elevation angle to satisfy 1) in the
            % instructions and later plot it
            ionogram_data = [ionogram_data; f, h_virtual, target_elev];
        end
    end
end

% ==========================================
% Plotting

% Ionogram
figure('Color', 'w'); % forcing to white bc I'm in dark mode
if ~isempty(ionogram_data)
    plot(ionogram_data(:,1), ionogram_data(:,2), 'b.-', 'LineWidth', 1.5);
    grid on;
    xlabel('Frequency (MHz)');
    ylabel('Virtual Reflection Height (km)');
    % title(['Oblique Ionogram: ' num2str(dist_km, '%.1f') ' km path'], 'Color','k');

    % switch statement to easily update title based on date
    switch UT(3)
        case 15
            % my birthday
            title(['8/15/2003 7:00am, Oblique Ionogram: ' num2str(dist_km, '%.1f') ' km path'], 'Color','k');
        case 16
            % the next day
            title(['8/16/2003 7:00am, Oblique Ionogram: ' num2str(dist_km, '%.1f') ' km path'], 'Color','k');
    end

    % switch statement to easily update subtitle based on transmitter
    switch transmit
        case "VA"
            % Chesapeake, VA
            subtitle('Tx: Chesapeake -> Rx: Auburn', 'Color','k');
        case "TX"
            % Corpus Christi, TX
            subtitle('Tx: Corpus Christi -> Rx: Auburn', 'Color','k');
    end

    xlim([1 10]);
    ylim([0 650]); % Ionogram height range from instructions
    set(gca, 'Color', 'w'); % Force axes background to white
    set(gca, 'XColor', 'k', 'YColor', 'k'); % Force axis lines to black
else
    fprintf('No rays reached the receiver! Check geometry or frequencies.\n');
end

% Elevation angle of the ray that reaches the receiver for each frequency
figure('Color', 'w'); % forcing to white bc I'm in dark mode
plot(ionogram_data(:,1), ionogram_data(:,3), 'b.-', 'LineWidth', 1.5);
grid on;
set(gca, 'Color', 'w'); % Force axes background to white
set(gca, 'XColor', 'k', 'YColor', 'k'); % Force axis lines to black
xlabel('Frequency (MHz)');
ylabel('Elevation Angle (deg)');
title('Elevation angle of the ray that reaches the receiver for each frequency', 'Color','k');

% ==========================================
% Extra code for checking the "rule of thumb"
% Commented out so it doesn't run every time

% Find maximum electron density in the grid
% N_max = max(iono_en_grid(:)); 
% fprintf('Max Electron Density (N_max): %.2e electrons/cm^3\n', N_max);

% Calculate f_c using the "rule of thumb" (approx 9*sqrt(N))
% convert N to SI (m^-3) for the standard formula f = 9*sqrt(N)
% N_cm3 * 1e6 = N_m3
% fc_Hz = 9 * sqrt(N_max * 1e6); 
% fc_MHz = fc_Hz / 1e6;
% fprintf('Theoretical Critical Frequency (fc): %.2f MHz\n', fc_MHz);