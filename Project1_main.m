clear;
close all;
clc;

% setup general stuff

% My birthday: August 15, 2003. Using 7am for the time.
UT = [2003 8 15 7 0];        % UT - year, month, day, hour, minute
R12 = 100;                   % R12 index
speed_of_light = 2.99792458e8;

% ray_bear = 324.7; % from example

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

% clear iri_options
% iri_options.Ne_B0B1_model = 'Bil-2000'; % this is a non-standard setting for 
%                                         % IRI but is used as an example

% Generate the ionospheric grid
tic
fprintf('Generating ionospheric grid...\n')
[iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
    gen_iono_grid_2d(origin_lat, origin_long, R12, UT, ray_bear, ...
    max_range, num_range, range_inc, start_height, ...
    height_inc, num_heights, kp, doppler_flag, 'iri2020');
toc
 

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;


%
% Example 1 - rays of different frequency with same launch elevation
% Changing from example to do project

% first call to raytrace so pass in the ionospheric and geomagnetic grids 
freqs = 1:0.1:10; % 1-10 MHz in 0.1 MHz increments
elevs = 3:0.2:85; % dense fan of elevation angles to hit the receiver
num_elevs = length(elevs);
tol = 1e-7;          % ODE tolerance
nhops = 2;

% I will call raytrace_2d in a loop instead
% tic
% fprintf('Generating %d 2D NRT rays ...', num_elevs);
% [ray_data, ray_path_data] = ...
%     raytrace_2d(origin_lat, origin_long, elevs, ray_bear, freqs, nhops, ...
%              tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
% 	     collision_freq, start_height, height_inc, range_inc, irreg);
% toc

% ==========================================
% MAIN SIMULATION LOOP

% initialize results storage: [frequency, virtual height, elev angle] eventually
ionogram_data = [];

fprintf('Starting Frequency Sweep (%.1f - %.1f MHz)...\n', 1, 10);

for f = freqs

    % run raytrace_2d for the whole fan of angles
    nhops = 1;
    freq_vec = ones(size(elevs)) * f;

    [ray_data, ~] = raytrace_2d(origin_lat, origin_long, elevs, ...
        ray_bear, freq_vec, nhops, tol, irregs_flag, ...
        iono_en_grid, iono_en_grid_5, ...
        collision_freq, start_height, height_inc, ...
        range_inc, irreg);

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
    plot(ionogram_data(:,1), ionogram_data(:,2), 'b.', 'LineWidth', 1.5);
    grid on;
    xlabel('Frequency (MHz)');
    ylabel('Virtual Reflection Height (km)');
    title(['Oblique Ionogram: ' num2str(dist_km, '%.1f') ' km path'], 'Color','k');

    % switch statement to easily update subtitle
    switch transmit
        case "VA"
            % Chesapeake, VA
            subtitle('Tx: Chesapeake -> Rx: Auburn', 'Color','k');
        case "TX"
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
