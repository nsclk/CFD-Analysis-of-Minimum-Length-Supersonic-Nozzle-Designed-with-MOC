%% Method of Characteristics for Minimum Length Nozzle Design
clear all
clc
tic

%% Gas properties %%
gamma = input('Enter gamma : ');

%% Design parameters %%
Me = input('Enter exit Mach no. : ');
theta_max = PrandtlMeyer(Me, gamma) * 180 / (2*pi);

%% Throat radius input %%
R_throat = input('Enter throat radius (mm, in, or any unit): ');
fprintf('Throat radius set to: %.4f units\n', R_throat);

%% Converging Section Design
fprintf('\n========== CONVERGING SECTION DESIGN ==========\n');
fprintf('Select converging section profile:\n');
fprintf('1. Circular arc\n');
fprintf('2. Witoszynski curve (adjustable exponents)\n');
fprintf('3. Bicubic polynomial (smooth slope)\n');
fprintf('4. Conical (linear)\n');
fprintf('5. Cosine bell\n');
fprintf('6. Quintic polynomial (NASA smooth)\n');
fprintf('7. No converging section (sharp corner only)\n');
conv_option = input('Enter option (1-7): ');

if conv_option ~= 7
    % ====== INPUT ======
    R_inlet = input('Enter inlet radius (typically 2-4 times throat radius): ');
    L_conv  = input('Enter converging section length: ');
    n_conv  = max(20, round(input('Number of points in converging section (>=20 recommended): ')));

    % ====== VALIDATION ======
    if R_inlet <= R_throat
        error('Inlet radius must be greater than throat radius.');
    end
    if L_conv <= 0
        error('Converging length must be positive.');
    end
    
    % ====== GRID ======
    x_conv = linspace(-L_conv, 0, n_conv);
    y_conv = zeros(1, n_conv);

    switch conv_option
        case 1  % === Circular arc ===
            R_curve  = ((L_conv)^2 + (R_inlet - R_throat)^2) / (2*(R_inlet - R_throat));
            y_center = R_inlet - R_curve;
            arg_sq   = R_curve^2 - (x_conv + L_conv).^2;
            if any(arg_sq < 0)
                error('Invalid geometry for circular arc — reduce L_conv or adjust radii.');
            end
            y_conv = y_center + sqrt(arg_sq);

        case 2  % === Adjustable Witoszynski curve ===
            n_exp = input('Enter exponent n (curvature control, e.g., 2): ');
            m_exp = input('Enter exponent m (bell shape control, e.g., 2): ');
            xi = (x_conv + L_conv) / L_conv;
            y_conv = R_throat + (R_inlet - R_throat) .* (1 - xi.^n_exp).^m_exp;

        case 3  % === Bicubic polynomial ===
            xi = (x_conv + L_conv) / L_conv;
            y_conv = R_throat + (R_inlet - R_throat) .* (1 - 3*xi.^2 + 2*xi.^3);

        case 4  % === Conical (linear) ===
            y_conv = linspace(R_inlet, R_throat, n_conv);

        case 5  % === Cosine bell profile ===
            xi = (x_conv + L_conv) / L_conv;
            y_conv = R_throat + (R_inlet - R_throat) .* (0.5*(1 + cos(pi*xi)));

        case 6  % === Quintic polynomial (NASA smooth) ===
            xi = (x_conv + L_conv) / L_conv;
            y_conv = R_throat + (R_inlet - R_throat) .* ...
                     (1 - 10*xi.^3 + 15*xi.^4 - 6*xi.^5); % C2 continuous
    end

    % ====== SAFETY CHECKS ======
    if any(y_conv <= 0)
        error('Negative or zero radius detected — check input values.');
    end

    % ====== PLOT PROFILE ======
    fig = figure;
    plot(x_conv, y_conv, 'b-', 'LineWidth', 2);
    xlabel('x (length units)');
    ylabel('Radius');
    title(sprintf('Converging Section: Option %d', conv_option));
    grid on;
    axis equal;
    legend('Converging profile', 'Location', 'Best');

    % ====== SAVE PNG IN REPORT FOLDER ======
    if exist('base_path', 'var') && ~isempty(base_path)
        fig_filename = fullfile(base_path, sprintf('Converging_Section_Option%d.png', conv_option));
    else
        fig_filename = sprintf('Converging_Section_Option%d.png', conv_option);
    end
    saveas(fig, fig_filename);
    fprintf('✅ Converging section figure saved as: %s\n', fig_filename);

else
    x_conv = [];
    y_conv = [];
    disp('No converging section: Sharp corner selected.');
end


%% Incident expansion wave conditions %%
n = input('Enter number of characteristics lines (greater than 2) emanating from sharp corner throat : ');
theta_0 = theta_max / n;

%% Characteristic parameter solver
% v - Prandtl-Meyer function
% KL - Left running characteristic constant
% KR - Right running characteristic constant
% theta - Flow angle relative to horizontal
[v, KL, KR, theta] = moc2d(theta_max, theta_0, n);

%% Mach number and Mach angle at each node
node = 0.5 * n * (4 + n - 1);
M = zeros(1, node);
mu = zeros(1, node);
for i = 1:node
    M(i) = flowprandtlmeyer(gamma, v(i), 'nu');
    mu(i) = Mu(M(i));
end

%% Grid generator
figure(1)
D = 1; % Non-Dimensional y co-ordinate of throat wall
i = 1;
x = zeros(1, node);
y = zeros(1, node);
wall = theta_max;

% Non-dimensional calculations first
while (i <= n+1)
    if i == 1
        x(i) = -D / (tand(theta(i) - mu(i)));
        y(i) = 0;
        plot([0 x(i)], [D 0]);
        hold on
    elseif i == n+1
        x(i) = (y(i-1) - D - x(i-1)*tand((theta(i-1) + theta(i) + mu(i-1) + mu(i))*0.5)) / ...
               (tand(0.5*(wall + theta(i))) - tand((theta(i-1) + theta(i) + mu(i-1) + mu(i))*0.5));
        y(i) = D + x(i)*tand(0.5*(wall + theta(i)));
        plot([x(i-1) x(i)], [y(i-1) y(i)]);
        hold on
        plot([0 x(i)], [D y(i)]);
        hold on
    else
        x(i) = (D - y(i-1) + x(i-1)*tand(0.5*(mu(i-1) + theta(i-1) + mu(i) + theta(i)))) / ...
               (tand(0.5*(mu(i-1) + theta(i-1) + mu(i) + theta(i))) - tand(theta(i) - mu(i)));
        y(i) = tand(theta(i) - mu(i))*x(i) + D;
        plot([x(i-1) x(i)], [y(i-1) y(i)]);
        hold on
        plot([0 x(i)], [D y(i)]);
        hold on
    end
    i = i + 1;
    hold on
end

h = i;
k = 0;
i = h;
for j = 1:n-1
    while (i <= h+n-k-1)
        if (i == h)
            x(i) = x(i-n+k) - y(i-n+k) / (tand(0.5*(theta(i-n+k) + theta(i) - mu(i-n+k) - mu(i))));
            y(i) = 0;
            plot([x(i-n+k) x(i)], [y(i-n+k) y(i)]);
            hold on
        elseif (i == h+n-k-1)
            x(i) = (x(i-n+k)*tand(0.5*(theta(i-n+k) + theta(i))) - y(i-n+k) + y(i-1) - ...
                   x(i-1)*tand((theta(i-1) + theta(i) + mu(i-1) + mu(i))*0.5)) / ...
                   (tand(0.5*(theta(i-n+k) + theta(i))) - tand((theta(i-1) + theta(i) + mu(i-1) + mu(i))*0.5));
            y(i) = y(i-n+k) + (x(i) - x(i-n+k))*tand(0.5*(theta(i-n+k) + theta(i)));
            plot([x(i-1) x(i)], [y(i-1) y(i)]);
            hold on
            plot([x(i-n+k) x(i)], [y(i-n+k) y(i)]);
            hold on
        else
            s1 = tand(0.5*(theta(i) + theta(i-1) + mu(i) + mu(i-1)));
            s2 = tand(0.5*(theta(i) + theta(i-n+k) - mu(i) - mu(i-n+k)));
            x(i) = (y(i-n+k) - y(i-1) + s1*x(i-1) - s2*x(i-n+k)) / (s1 - s2);
            y(i) = y(i-1) + (x(i) - x(i-1))*s1;
            plot([x(i-1) x(i)], [y(i-1) y(i)]);
            hold on
            plot([x(i-n+k) x(i)], [y(i-n+k) y(i)]);
            hold on
        end
        i = i + 1;
    end
    k = k + 1;
    h = i;
    i = h;
    hold on
end

% Scale all coordinates by throat radius
x_scaled = x * R_throat;
y_scaled = y * R_throat;

title(sprintf('Characteristic lines for Mach=%d and gamma=%d', Me, gamma))
xlabel('x/x0')
ylabel('y/y0')
axis equal
hold on
plot(xlim, [0 0], '--k')                 % Dashed Horizontal Line at y=0
hold off 
xlim()
ylim()


%% Nozzle co-ordinates (x_wall, y_wall)
x_wall = zeros(1, n+1);
y_wall = zeros(1, n+1);
x_wall(1) = 0;
y_wall(1) = 1;
i = 2;
j = n + 1;
% Store the indices of wall nodes for later use
wall_node_indices = zeros(1, n+1);
wall_node_indices(1) = n+1;  % First wall point corresponds to node n+1
while (i <= n+1)
    x_wall(i) = x(j);
    y_wall(i) = y(j);
    wall_node_indices(i) = j;  % Store the node index
    j = j + (n - i + 2);
    i = i + 1;
end

% Scale wall coordinates by throat radius
x_wall_scaled = x_wall * R_throat;
y_wall_scaled = y_wall * R_throat;

%% Characteristic Lines with Real Dimensions
figure(4)
i = 1;

% Plot characteristic lines with scaled coordinates
while (i <= n+1)
    if i == 1
        x1_real = 0;
        y1_real = R_throat;
        x2_real = x_scaled(i);
        y2_real = 0;
        plot([x1_real x2_real], [y1_real y2_real], 'b-');
        hold on
    elseif i == n+1
        x1_real = x_scaled(i-1);
        y1_real = y_scaled(i-1);
        x2_real = x_scaled(i);
        y2_real = y_scaled(i);
        plot([x1_real x2_real], [y1_real y2_real], 'b-');
        hold on
        plot([0 x2_real], [R_throat y2_real], 'b-');
        hold on
    else
        x1_real = x_scaled(i-1);
        y1_real = y_scaled(i-1);
        x2_real = x_scaled(i);
        y2_real = y_scaled(i);
        plot([x1_real x2_real], [y1_real y2_real], 'b-');
        hold on
        plot([0 x2_real], [R_throat y2_real], 'b-');
        hold on
    end
    i = i + 1;
end

h = i;
k = 0;
i = h;
for j = 1:n-1
    while (i <= h+n-k-1)
        if (i == h)
            x1_real = x_scaled(i-n+k);
            y1_real = y_scaled(i-n+k);
            x2_real = x_scaled(i);
            y2_real = y_scaled(i);
            plot([x1_real x2_real], [y1_real y2_real], 'b-');
            hold on
        elseif (i == h+n-k-1)
            % Left running characteristic
            x1_real = x_scaled(i-1);
            y1_real = y_scaled(i-1);
            x2_real = x_scaled(i);
            y2_real = y_scaled(i);
            plot([x1_real x2_real], [y1_real y2_real], 'b-');
            hold on
            % Right running characteristic
            x1_real = x_scaled(i-n+k);
            y1_real = y_scaled(i-n+k);
            plot([x1_real x2_real], [y1_real y2_real], 'b-');
            hold on
        else
            % Left running characteristic from previous point
            x1_real = x_scaled(i-1);
            y1_real = y_scaled(i-1);
            x2_real = x_scaled(i);
            y2_real = y_scaled(i);
            plot([x1_real x2_real], [y1_real y2_real], 'b-');
            hold on
            % Right running characteristic from upper point
            x1_real = x_scaled(i-n+k);
            y1_real = y_scaled(i-n+k);
            plot([x1_real x2_real], [y1_real y2_real], 'b-');
            hold on
        end
        i = i + 1;
    end
    k = k + 1;
    h = i;
    i = h;
end

% Plot the nozzle walls with real dimensions
plot(x_wall_scaled, y_wall_scaled, 'k-', 'LineWidth', 1)
plot(x_wall_scaled, -y_wall_scaled, 'k-', 'LineWidth', 1)

% Add centerline
plot(xlim, [0 0], '--k', 'LineWidth', 1)

% Labels and formatting
title(sprintf('Characteristic Lines with Real Dimensions (Me = %.2f, γ = %.2f, R_{throat} = %.2f units)', Me, gamma, R_throat))
xlabel(sprintf('Axial Distance, x (%s)', 'units'))
ylabel(sprintf('Radial Distance, y (%s)', 'units'))
axis equal
grid on

% Add annotations
text(0.1*R_throat, R_throat*1.1, 'Throat', 'FontSize', 10, 'FontWeight', 'bold')
text(x_wall_scaled(end)*0.9, y_wall_scaled(end)*1.1, 'Exit', 'FontSize', 10, 'FontWeight', 'bold')

% Add node numbers to Figure 4
for i = 1:node
    text(x_scaled(i), y_scaled(i), num2str(i), 'FontSize', 7, 'Color', 'k', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
end

hold off

% Combine converging and diverging sections
if conv_option ~= 5
    x_complete = [x_conv, x_wall_scaled];
    y_complete = [y_conv, y_wall_scaled];
else
    x_complete = x_wall_scaled;
    y_complete = y_wall_scaled;
end

figure(2)
if conv_option ~= 5
    % Plot converging section
    plot(x_conv, y_conv, 'r-', 'LineWidth', 2)
    hold on
    plot(x_conv, -y_conv, 'r-', 'LineWidth', 2)
    % Plot diverging section
    plot(x_wall_scaled, y_wall_scaled, 'b-', 'LineWidth', 2)
    plot(x_wall_scaled, -y_wall_scaled, 'b-', 'LineWidth', 2)
    legend('Converging (upper)', 'Converging (lower)', 'Diverging (upper)', 'Diverging (lower)', 'Location', 'best')
else
    plot(x_wall_scaled, y_wall_scaled, 'b-', 'LineWidth', 2)
    plot(x_wall_scaled, -y_wall_scaled, 'b-', 'LineWidth', 2)
    legend('Upper Wall', 'Lower Wall', 'Location', 'best')
end

title(sprintf('Complete Nozzle Contour (Me = %.2f, R_{throat} = %.2f units)', Me, R_throat))
xlabel('x (units)')
ylabel('y (units)')
axis equal
grid on

% Add markers
if conv_option ~= 5
    % Inlet markers
    plot(x_conv(1), y_conv(1), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm')
    plot(x_conv(1), -y_conv(1), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm')
    text(x_conv(1)-0.2*R_throat, y_conv(1)+0.1*R_throat, 'Inlet', 'FontSize', 10)
end

% Throat markers
plot(0, R_throat, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
plot(0, -R_throat, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
text(0.1*R_throat, R_throat+0.1*R_throat, 'Throat', 'FontSize', 10)

% Exit markers
plot(x_wall_scaled(end), y_wall_scaled(end), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g')
plot(x_wall_scaled(end), -y_wall_scaled(end), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g')
text(x_wall_scaled(end)-0.3*R_throat, y_wall_scaled(end)+0.1*R_throat, 'Exit', 'FontSize', 10)

% Add a third figure showing Mach number distribution
figure(3)
if conv_option ~= 5
    % Approximate Mach number in converging section (subsonic)
    M_conv = linspace(0.1, 1.0, n_conv);
    x_div_mach = x_wall_scaled;
    % Get Mach numbers at wall points using the stored indices
    M_div = zeros(1, length(wall_node_indices));
    for i = 1:length(wall_node_indices)
        if wall_node_indices(i) > 0 && wall_node_indices(i) <= length(M)
            M_div(i) = M(wall_node_indices(i));
        else
            % Handle the special case for the throat (first wall point)
            M_div(i) = 1.0;  % Mach 1 at throat
        end
    end
    
    plot(x_conv, M_conv, 'r-', 'LineWidth', 2)
    hold on
    plot(x_div_mach, M_div, 'b-', 'LineWidth', 2)
    
    % Add node numbers for wall points in diverging section
    for i = 1:length(wall_node_indices)
        text(x_div_mach(i), M_div(i), sprintf('N%d', wall_node_indices(i)), ...
             'FontSize', 8, 'Color', 'blue', 'VerticalAlignment', 'bottom', ...
             'HorizontalAlignment', 'center')
    end
    
    xlabel('x (units)')
    ylabel('Mach Number')
    title('Mach Number Distribution Along Nozzle')
    grid on
    legend('Converging Section', 'Diverging Section', 'Location', 'best')
    
    % Add sonic line
    plot([0 0], [0 max(M_div)*1.1], 'k--', 'LineWidth', 1)
    text(0.1*R_throat, 1.1, 'M=1', 'FontSize', 10)
else
    x_div_mach = x_wall_scaled;
    % Get Mach numbers at wall points using the stored indices
    M_div = zeros(1, length(wall_node_indices));
    for i = 1:length(wall_node_indices)
        if wall_node_indices(i) > 0 && wall_node_indices(i) <= length(M)
            M_div(i) = M(wall_node_indices(i));
        else
            % Handle the special case for the throat (first wall point)
            M_div(i) = 1.0;  % Mach 1 at throat
        end
    end
    
    plot(x_div_mach, M_div, 'b-', 'LineWidth', 2)
    
    % Add node numbers for wall points
    for i = 1:length(wall_node_indices)
        text(x_div_mach(i), M_div(i), sprintf('N%d', wall_node_indices(i)), ...
             'FontSize', 8, 'Color', 'blue', 'VerticalAlignment', 'bottom', ...
             'HorizontalAlignment', 'center')
    end
    
    xlabel('x (units)')
    ylabel('Mach Number')
    title('Mach Number Distribution in Diverging Section')
    grid on
end

%% Display Results in Tables

% Nozzle Design Summary Table
fprintf('\n========== NOZZLE DESIGN SUMMARY ==========\n');
design_params = {'Design Parameter', 'Value', 'Units';
                 'Gamma', gamma, '-';
                 'Exit Mach Number', Me, '-';
                 'Throat Radius', R_throat, 'units';
                 'Number of Characteristics', n, '-';
                 'Maximum Flow Angle', theta_max, 'degrees';
                 'Initial Flow Angle', theta_0, 'degrees';
                 'Total Number of Nodes', node, '-'};

fprintf('%-30s %-15s %-10s\n', design_params{1,:});
fprintf('%-30s %-15s %-10s\n', repmat('-',1,30), repmat('-',1,15), repmat('-',1,10));
for i = 2:size(design_params,1)
    fprintf('%-30s %-15.4f %-10s\n', design_params{i,:});
end

% Nozzle Geometry Table (Scaled)
fprintf('\n\n========== NOZZLE GEOMETRY (SCALED) ==========\n');
fprintf('%-10s %-15s %-15s\n', 'Point', 'X', 'Y');
fprintf('%-10s %-15s %-15s\n', repmat('-',1,10), repmat('-',1,15), repmat('-',1,15));
for i = 1:length(x_wall_scaled)
    fprintf('%-10d %-15.6f %-15.6f\n', i, x_wall_scaled(i), y_wall_scaled(i));
end

% Key Nozzle Features
L_nozzle = x_wall_scaled(end);
R_exit = y_wall_scaled(end);
R_throat_actual = y_wall_scaled(1);
area_ratio = (R_exit/R_throat_actual)^2;
D_throat = 2 * R_throat_actual;
D_exit = 2 * R_exit;
A_throat = pi * R_throat_actual^2;
A_exit = pi * R_exit^2;

% Additional features for converging section
if conv_option ~= 5
    L_total = L_conv + L_nozzle;
    R_inlet_actual = y_conv(1);
    D_inlet = 2 * R_inlet_actual;
    A_inlet = pi * R_inlet_actual^2;
    contraction_ratio = A_inlet / A_throat;
else
    L_total = L_nozzle;
end

fprintf('\n\n========== KEY NOZZLE FEATURES ==========\n');
if conv_option ~= 5
    features = {'Feature', 'Value', 'Units';
                'Total Nozzle Length', L_total, 'units';
                'Converging Section Length', L_conv, 'units';
                'Diverging Section Length', L_nozzle, 'units';
                'Inlet Radius', R_inlet_actual, 'units';
                'Inlet Diameter', D_inlet, 'units';
                'Inlet Area', A_inlet, 'units^2';
                'Throat Radius', R_throat_actual, 'units';
                'Throat Diameter', D_throat, 'units';
                'Throat Area', A_throat, 'units^2';
                'Exit Radius', R_exit, 'units';
                'Exit Diameter', D_exit, 'units';
                'Exit Area', A_exit, 'units^2';
                'Contraction Ratio (A_inlet/A_throat)', contraction_ratio, '-';
                'Area Ratio (A_exit/A_throat)', area_ratio, '-';
                'Expansion Ratio (R_exit/R_throat)', R_exit/R_throat_actual, '-'};
else
    features = {'Feature', 'Value', 'Units';
                'Nozzle Length', L_nozzle, 'units';
                'Throat Radius', R_throat_actual, 'units';
                'Throat Diameter', D_throat, 'units';
                'Throat Area', A_throat, 'units^2';
                'Exit Radius', R_exit, 'units';
                'Exit Diameter', D_exit, 'units';
                'Exit Area', A_exit, 'units^2';
                'Area Ratio (A_exit/A_throat)', area_ratio, '-';
                'Expansion Ratio (R_exit/R_throat)', R_exit/R_throat_actual, '-'};
end

fprintf('%-30s %-15s %-10s\n', features{1,:});
fprintf('%-30s %-15s %-10s\n', repmat('-',1,30), repmat('-',1,15), repmat('-',1,10));
for i = 2:size(features,1)
    fprintf('%-30s %-15.4f %-10s\n', features{i,:});
end

% Selected Node Properties Table (first 10 and last 10 nodes)
fprintf('\n\n========== SELECTED NODE PROPERTIES ==========\n');
fprintf('%-10s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
        'Node', 'X', 'Y', 'Theta(deg)', 'Mach', 'Mu(deg)', 'PM Angle');
fprintf('%s\n', repmat('-',1,90));

% Display first 10 nodes (with scaled coordinates)
num_display = min(65, node);
for i = 1:num_display
    fprintf('%-10d %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f\n', ...
            i, x_scaled(i), y_scaled(i), theta(i), M(i), mu(i), v(i));
end

% Display last 10 nodes if more than 20 nodes total
%if node > 20
%    fprintf('%-10s %-12s %-12s %-12s %-12s %-12s %-12s\n', '...', '...', '...', '...', '...', '...', '...');
%    for i = max(node-9, num_display+1):node
%        fprintf('%-10d %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f\n', ...
%                i, x_scaled(i), y_scaled(i), theta(i), M(i), mu(i), v(i));
%    end
%end

%% Save Results to Files  
save_option = input('\nSave results to file? (y/n): ', 's');

if strcmpi(save_option, 'y')
    % Get filename and create output directory
    filename = input('Enter filename (without extension): ', 's');
    
    % Option to create a dedicated folder for outputs
    create_folder = input('Create dedicated folder for outputs? (y/n): ', 's');
    if strcmpi(create_folder, 'y')
        output_dir = [filename '_output'];
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        base_path = fullfile(output_dir, filename);
    else
        base_path = filename;
    end
    
    % Initialize file list for summary
    saved_files = {};
    
    try
        %% 1. Save MATLAB Workspace
        mat_file = [base_path '.mat'];
        save(mat_file);
        saved_files{end+1} = sprintf('%s.mat (workspace variables)', filename);
        
        %% 2. Save Coordinate Files (CSV format with headers)
        if conv_option ~= 5
            % Complete nozzle coordinates
            saveCoordinatesToCSV([base_path '_complete_coordinates.csv'], ...
                                x_complete', y_complete', 'Complete Nozzle');
            saved_files{end+1} = sprintf('%s_complete_coordinates.csv (full nozzle)', filename);
            
            % Converging section
            saveCoordinatesToCSV([base_path '_converging_coordinates.csv'], ...
                                x_conv', y_conv', 'Converging Section');
            saved_files{end+1} = sprintf('%s_converging_coordinates.csv (converging section)', filename);
            
            % Diverging section
            saveCoordinatesToCSV([base_path '_diverging_coordinates.csv'], ...
                                x_wall_scaled', y_wall_scaled', 'Diverging Section');
            saved_files{end+1} = sprintf('%s_diverging_coordinates.csv (diverging section)', filename);
        else
            % Diverging section only
            saveCoordinatesToCSV([base_path '_coordinates.csv'], ...
                                x_wall_scaled', y_wall_scaled', 'Nozzle Wall');
            saved_files{end+1} = sprintf('%s_coordinates.csv (nozzle coordinates)', filename);
        end
        
        %% 3. Save ALL Node Properties to CSV (Complete Flow Field Data)
        node_file = [base_path '_all_node_properties.csv'];
        saveAllNodeProperties(node_file, 1:node, x_scaled, y_scaled, theta, M, mu, v, KL, KR);
        saved_files{end+1} = sprintf('%s_all_node_properties.csv (complete flow field data)', filename);
        
        %% 4. Save Wall Node Properties
        wall_nodes_file = [base_path '_wall_node_properties.csv'];
        saveWallNodeProperties(wall_nodes_file, wall_node_indices, x_wall_scaled, y_wall_scaled, ...
                              theta, M, mu, v, KL, KR);
        saved_files{end+1} = sprintf('%s_wall_node_properties.csv (wall flow properties)', filename);
        
        %% 5. Save Key Nozzle Features
        features_file = [base_path '_nozzle_features.csv'];
        saveNozzleFeatures(features_file, gamma, Me, R_throat, R_exit, L_nozzle, ...
                          D_throat, D_exit, A_throat, A_exit, area_ratio, ...
                          theta_max, theta_0, n, node, conv_option, ...
                          R_inlet, L_conv, L_total, exist('contraction_ratio','var')*contraction_ratio + ~exist('contraction_ratio','var')*NaN);
        saved_files{end+1} = sprintf('%s_nozzle_features.csv (key design parameters)', filename);
        
        %% 6. Generate Detailed Report with ALL Information  (R_exit added)  %%% UPDATED %%%
        report_file = [base_path '_detailed_report.txt'];
        generateCompleteReport(report_file, gamma, Me, R_throat, R_exit, n, theta_max, theta_0, ...
                              conv_option, R_inlet, L_conv, L_nozzle, L_total, ...
                              D_throat, D_exit, A_throat, A_exit, area_ratio, contraction_ratio, ...
                              x_complete, y_complete, x_wall_scaled, y_wall_scaled, ...
                              x_scaled, y_scaled, theta, M, mu, v, KL, KR, node, ...
                              wall_node_indices);
        saved_files{end+1} = sprintf('%s_detailed_report.txt (comprehensive report)', filename);
        
        %% 7. Save Summary File with Key Metrics
        summary_file = [base_path '_summary.txt'];
        generateEnhancedSummary(summary_file, gamma, Me, R_throat, R_exit, area_ratio, ...
                               L_nozzle, L_total, D_throat, D_exit, A_throat, A_exit, ...
                               contraction_ratio, conv_option, n, node);
        saved_files{end+1} = sprintf('%s_summary.txt (executive summary)', filename);
        
        %% 8. Save Performance Metrics
        performance_file = [base_path '_performance_metrics.csv'];
        savePerformanceMetrics(performance_file, gamma, Me, area_ratio, L_nozzle, ...
                              R_throat, R_exit, theta_max, n, node);
        saved_files{end+1} = sprintf('%s_performance_metrics.csv (performance data)', filename);
        
        %% 9. Export Figures (ALL OPEN FIGURES, USE FIGURE NAMES, PNG, SAME FOLDER AS REPORT)
        export_figs = input('Export figures as images? (y/n): ', 's');
        if strcmpi(export_figs, 'y')
            figs = findall(0, 'Type', 'figure');   % Find all open figures
            if isempty(figs)
                fprintf('No open figures found.\n');
            else
                % Sort by figure number (creation order)
                [~, idx] = sort(arrayfun(@(f) get(f,'Number'), figs));
                figs = figs(idx);
        
                for k = 1:numel(figs)
                    fh = figs(k);
        
                    % Get figure name, fallback to Figure_<n>
                    nm = strtrim(get(fh,'Name'));
                    if isempty(nm)
                        nm = sprintf('Figure_%d', get(fh,'Number'));
                    end
        
                    % Clean invalid filename characters
                    nm = regexprep(nm, '[^\w\s-]', '_');
        
                    % Save in same folder as report
                    out = fullfile(fileparts(base_path), sprintf('%s_%s.png', filename, nm));
        
                    % Save as PNG (300 dpi) with fallbacks
                    try
                        exportgraphics(fh, out, 'Resolution', 300);
                    catch
                        try
                            print(fh, out, '-dpng', '-r300');
                        catch
                            saveas(fh, out);
                        end
                    end
        
                    fprintf('Figure saved: %s\n', out);  % Show full path
                    saved_files{end+1} = out; 
                end
            end
        end



        
        %% Display Success Message
        fprintf('\n================================================================================\n');
        fprintf('                         FILES SAVED SUCCESSFULLY!                              \n');
        fprintf('================================================================================\n');
        if exist('output_dir', 'var') && exist(output_dir, 'dir')
            fprintf('Output directory: %s\n', output_dir);
        end
        fprintf('\nSaved files:\n');
        for i = 1:length(saved_files)
            fprintf('  ✓ %s\n', saved_files{i});
        end
        fprintf('--------------------------------------------------------------------------------\n');
        fprintf('Total files saved: %d\n', length(saved_files));
        fprintf('================================================================================\n');
        
    catch ME
        fprintf('\n⚠ ERROR: Failed to save files!\n');
        fprintf('Error message: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('Error location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
        fprintf('Please check file permissions and disk space.\n');
    end
end


%% Helper Functions for File Saving  

function saveCoordinatesToCSV(filename, x_data, y_data, description)
    % Save coordinates with headers and metadata
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot create file: %s', filename);
    end
    
    % Write header
    fprintf(fid, '# %s Coordinates\n', description);
    fprintf(fid, '# Generated: %s\n', datestr(now));
    fprintf(fid, '# Number of Points: %d\n', length(x_data));
    fprintf(fid, '# Units: User-defined\n');
    fprintf(fid, '# Format: X,Y\n');
    fprintf(fid, 'Point_Number,X,Y\n');
    
    % Write data
    for i = 1:length(x_data)
        fprintf(fid, '%d,%.10f,%.10f\n', i, x_data(i), y_data(i));
    end
    
    fclose(fid);
end

function saveAllNodeProperties(filename, node_ids, x, y, theta, M, mu, v, KL, KR)
    % Save ALL node properties to CSV with comprehensive data
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot create file: %s', filename);
    end
    
    % Write metadata
    fprintf(fid, '# Complete Flow Field Node Properties\n');
    fprintf(fid, '# Generated: %s\n', datestr(now));
    fprintf(fid, '# Total Nodes: %d\n', length(node_ids));
    fprintf(fid, '# Angles in degrees, Prandtl-Meyer angle in radians\n');
    
    % Write header
    fprintf(fid, 'Node,X,Y,Theta_deg,Mach,Mu_deg,PM_Angle,KL,KR,Velocity_Ratio,Pressure_Ratio,Temp_Ratio,Density_Ratio\n');
    
    % Calculations (use gamma=1.4 for these extra ratios)
    gamma_calc = 1.4;
    for i = 1:length(node_ids)
        V_ratio   = M(i) * sqrt((2/(gamma_calc+1)) * (1 + ((gamma_calc-1)/2)*M(i)^2)^(-1));
        P_ratio   = (1 + ((gamma_calc-1)/2)*M(i)^2)^(-gamma_calc/(gamma_calc-1));
        T_ratio   = (1 + ((gamma_calc-1)/2)*M(i)^2)^(-1);
        rho_ratio = (1 + ((gamma_calc-1)/2)*M(i)^2)^(-1/(gamma_calc-1));
        
        fprintf(fid, '%d,%.8f,%.8f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
                node_ids(i), x(i), y(i), theta(i), M(i), mu(i), v(i), ...
                KL(i), KR(i), V_ratio, P_ratio, T_ratio, rho_ratio);
    end
    
    fclose(fid);
end

function saveWallNodeProperties(filename, wall_indices, x_wall, y_wall, theta, M, mu, v, KL, KR)
    % Save wall node properties specifically
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot create file: %s', filename);
    end
    
    % Write metadata
    fprintf(fid, '# Wall Node Flow Properties\n');
    fprintf(fid, '# Generated: %s\n', datestr(now));
    fprintf(fid, '# Number of Wall Points: %d\n', length(x_wall));
    
    % Write header
    fprintf(fid, 'Wall_Point,Node_Index,X,Y,Mach,Theta_deg,Mu_deg,PM_Angle,KL,KR\n');
    
    % Write wall node data
    for i = 1:length(wall_indices)
        idx = wall_indices(i);
        if idx > 0 && idx <= length(M)
            fprintf(fid, '%d,%d,%.8f,%.8f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
                    i, idx, x_wall(i), y_wall(i), M(idx), theta(idx), ...
                    mu(idx), v(idx), KL(idx), KR(idx));
        else
            % Handle throat point (M=1), angles at throat set conventionally
            fprintf(fid, '%d,%d,%.8f,%.8f,1.0000,0.0000,90.0000,0.0000,0.0000,0.0000\n', ...
                    i, max(idx,0), x_wall(i), y_wall(i));
        end
    end
    
    fclose(fid);
end

function saveNozzleFeatures(filename, gamma, Me, R_throat, R_exit, L_nozzle, ...
                           D_throat, D_exit, A_throat, A_exit, area_ratio, ...
                           theta_max, theta_0, n, node, conv_option, ...
                           R_inlet, L_conv, L_total, contraction_ratio)
    % Save key nozzle features to CSV
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot create file: %s', filename);
    end
    
    % Write header
    fprintf(fid, '# Key Nozzle Design Features\n');
    fprintf(fid, '# Generated: %s\n', datestr(now));
    fprintf(fid, 'Parameter,Value,Units,Category\n');
    
    % Gas Properties
    fprintf(fid, 'Gamma,%.6f,-,Gas Properties\n', gamma);
    
    % Design Parameters
    fprintf(fid, 'Exit Mach Number,%.6f,-,Design Parameters\n', Me);
    fprintf(fid, 'Number of Characteristics,%d,-,Design Parameters\n', n);
    fprintf(fid, 'Total Number of Nodes,%d,-,Design Parameters\n', node);
    fprintf(fid, 'Maximum Flow Angle,%.6f,degrees,Design Parameters\n', theta_max);
    fprintf(fid, 'Initial Flow Angle,%.6f,degrees,Design Parameters\n', theta_0);
    
    % Throat Properties
    fprintf(fid, 'Throat Radius,%.6f,units,Throat Geometry\n', R_throat);
    fprintf(fid, 'Throat Diameter,%.6f,units,Throat Geometry\n', D_throat);
    fprintf(fid, 'Throat Area,%.6f,units^2,Throat Geometry\n', A_throat);
    
    % Exit Properties
    fprintf(fid, 'Exit Radius,%.6f,units,Exit Geometry\n', R_exit);
    fprintf(fid, 'Exit Diameter,%.6f,units,Exit Geometry\n', D_exit);
    fprintf(fid, 'Exit Area,%.6f,units^2,Exit Geometry\n', A_exit);
    
    % Geometric Ratios
    fprintf(fid, 'Area Ratio (Ae/At),%.6f,-,Geometric Ratios\n', area_ratio);
    fprintf(fid, 'Expansion Ratio (Re/Rt),%.6f,-,Geometric Ratios\n', R_exit/R_throat);
    
    % Lengths
    fprintf(fid, 'Diverging Section Length,%.6f,units,Lengths\n', L_nozzle);
    
    % Converging Section (if applicable)
    if conv_option ~= 5
        profile_names = {'Circular Arc', 'Witoszynski Curve', 'Bicubic Polynomial', 'Conical'};
        fprintf(fid, 'Converging Profile,%s,-,Converging Section\n', profile_names{conv_option});
        fprintf(fid, 'Inlet Radius,%.6f,units,Converging Section\n', R_inlet);
        fprintf(fid, 'Converging Length,%.6f,units,Converging Section\n', L_conv);
        fprintf(fid, 'Total Nozzle Length,%.6f,units,Converging Section\n', L_total);
        fprintf(fid, 'Contraction Ratio (Ai/At),%.6f,-,Converging Section\n', contraction_ratio);
        fprintf(fid, 'Inlet Diameter,%.6f,units,Converging Section\n', 2*R_inlet);
        fprintf(fid, 'Inlet Area,%.6f,units^2,Converging Section\n', pi*R_inlet^2);
    else
        fprintf(fid, 'Converging Profile,None,-,Converging Section\n');
        fprintf(fid, 'Total Nozzle Length,%.6f,units,Lengths\n', L_nozzle);
    end
    
    % Performance Metrics
    fprintf(fid, 'Length to Throat Ratio,%.6f,-,Performance Metrics\n', L_nozzle/R_throat);
    fprintf(fid, 'Exit to Throat Diameter Ratio,%.6f,-,Performance Metrics\n', D_exit/D_throat);
    
    if Me > 1
        Cf_ideal = sqrt(2*gamma^2/(gamma-1) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * ...
                       (1 - (1/Me^2)^((gamma-1)/gamma)));
        fprintf(fid, 'Ideal Thrust Coefficient,%.6f,-,Performance Metrics\n', Cf_ideal);
    end
    
    fclose(fid);
end

function generateCompleteReport(filename, gamma, Me, R_throat, R_exit, n, theta_max, theta_0, ...
                               conv_option, R_inlet, L_conv, L_nozzle, L_total, ...
                               D_throat, D_exit, A_throat, A_exit, area_ratio, contraction_ratio, ...
                               x_complete, y_complete, x_wall_scaled, y_wall_scaled, ...
                               x_scaled, y_scaled, theta, M, mu, v, KL, KR, node, ...
                               wall_node_indices)
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot create report file: %s', filename);
    end
    
    % Header
    fprintf(fid, '********************************************************************************\n');
    fprintf(fid, '*                    MINIMUM LENGTH NOZZLE DESIGN REPORT                      *\n');
    fprintf(fid, '*                      METHOD OF CHARACTERISTICS (MOC)                        *\n');
    fprintf(fid, '********************************************************************************\n');
    fprintf(fid, 'Generated on: %s\n', datestr(now));
    fprintf(fid, 'MATLAB Version: %s\n\n', version);
    
    % Executive Summary
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '                              EXECUTIVE SUMMARY                                \n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'This nozzle design achieves Mach %.2f with minimum length using MOC.\n', Me);
    fprintf(fid, 'The design ensures uniform, parallel, shock-free flow at the exit.\n');
    fprintf(fid, 'Area Ratio: %.3f | Length: %.3f units | Throat Radius: %.3f units\n\n', ...
            area_ratio, L_nozzle, R_throat);
    
    % Design Parameters Section
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '                              DESIGN PARAMETERS                                \n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '%-40s : %.6f\n', 'Specific Heat Ratio (γ)', gamma);
    fprintf(fid, '%-40s : %.6f\n', 'Exit Mach Number (Me)', Me);
    fprintf(fid, '%-40s : %.6f units\n', 'Throat Radius', R_throat);
    fprintf(fid, '%-40s : %d\n', 'Number of Characteristics (n)', n);
    fprintf(fid, '%-40s : %d\n', 'Total Number of Nodes', node);
    fprintf(fid, '%-40s : %.6f degrees\n', 'Maximum Flow Angle (θ_max)', theta_max);
    fprintf(fid, '%-40s : %.6f degrees\n', 'Initial Flow Angle (θ_0)', theta_0);
    fprintf(fid, '%-40s : %.6f degrees\n', 'Flow Angle Increment (Δθ)', theta_0);
    fprintf(fid, '\n');
    
    % Converging Section
    if conv_option ~= 5
        fprintf(fid, '================================================================================\n');
        fprintf(fid, '                             CONVERGING SECTION                                \n');
        fprintf(fid, '================================================================================\n');
        profile_names = {'Circular Arc', 'Witoszynski Curve (Bell-mouth)', ...
                         'Bicubic Polynomial', 'Conical (Linear)'};
        fprintf(fid, '%-40s : %s\n', 'Profile Type', profile_names{conv_option});
        fprintf(fid, '%-40s : %.6f units\n', 'Inlet Radius', R_inlet);
        fprintf(fid, '%-40s : %.6f units\n', 'Inlet Diameter', 2*R_inlet);
        fprintf(fid, '%-40s : %.6f units²\n', 'Inlet Area', pi*R_inlet^2);
        fprintf(fid, '%-40s : %.6f units\n', 'Converging Section Length', L_conv);
        fprintf(fid, '%-40s : %.6f\n', 'Contraction Ratio (A_inlet/A_throat)', contraction_ratio);
        fprintf(fid, '\n');
    end
    
    % Nozzle Geometry
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '                              NOZZLE GEOMETRY                                  \n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'DIMENSIONAL PARAMETERS:\n');
    fprintf(fid, '------------------------\n');
    if conv_option ~= 5
        fprintf(fid, '%-40s : %.6f units\n', 'Total Nozzle Length', L_total);
        fprintf(fid, '%-40s : %.6f units\n', 'Converging Section Length', L_conv);
    end
    fprintf(fid, '%-40s : %.6f units\n', 'Diverging Section Length', L_nozzle);
    fprintf(fid, '%-40s : %.6f units\n', 'Throat Radius', R_throat);
    fprintf(fid, '%-40s : %.6f units\n', 'Throat Diameter', D_throat);
    fprintf(fid, '%-40s : %.6f units²\n', 'Throat Area', A_throat);
    fprintf(fid, '%-40s : %.6f units\n', 'Exit Radius', R_exit);
    fprintf(fid, '%-40s : %.6f units\n', 'Exit Diameter', D_exit);
    fprintf(fid, '%-40s : %.6f units²\n', 'Exit Area', A_exit);
    fprintf(fid, '\n');
    
    fprintf(fid, 'NON-DIMENSIONAL PARAMETERS:\n');
    fprintf(fid, '----------------------------\n');
    fprintf(fid, '%-40s : %.6f\n', 'Area Ratio (A_exit/A_throat)', area_ratio);
    fprintf(fid, '%-40s : %.6f\n', 'Expansion Ratio (R_exit/R_throat)', R_exit/R_throat);
    fprintf(fid, '%-40s : %.6f\n', 'Length to Throat Radius Ratio', L_nozzle/R_throat);
    fprintf(fid, '%-40s : %.6f\n', 'Exit to Throat Diameter Ratio', D_exit/D_throat);
    fprintf(fid, '\n');
    
    % Flow Properties at Key Locations
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '                        FLOW PROPERTIES AT KEY LOCATIONS                       \n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'THROAT (x = 0):\n');
    fprintf(fid, '  Mach Number                            : 1.0000\n');
    fprintf(fid, '  Flow Angle                             : 0.0000 degrees\n');
    fprintf(fid, '  Radius                                 : %.6f units\n', R_throat);
    fprintf(fid, '\n');
    fprintf(fid, 'EXIT (x = %.4f):\n', x_wall_scaled(end));
    fprintf(fid, '  Mach Number                            : %.6f\n', Me);
    fprintf(fid, '  Flow Angle                             : 0.0000 degrees (parallel flow)\n');
    fprintf(fid, '  Radius                                 : %.6f units\n', R_exit);
    fprintf(fid, '\n');
    
    % Sample Node Properties  (now showing ALL nodes)
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '                           NODE PROPERTIES (ALL NODES)                         \n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '  Node      X          Y       Theta(°)    Mach      μ(°)     ν(rad)     KL        KR\n');
    fprintf(fid, '------  --------  --------  ----------  -------  --------  --------  --------  --------\n');

    for idx = 1:node
        fprintf(fid, '%6d  %8.4f  %8.4f  %10.4f  %7.4f  %8.4f  %8.4f  %8.4f  %8.4f\n', ...
               idx, x_scaled(idx), y_scaled(idx), theta(idx), M(idx), ...
                mu(idx), v(idx), KL(idx), KR(idx));
    end

    fprintf(fid, '\n');

    
   % Wall Coordinates Section (show ALL points, no sampling)
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '                             WALL COORDINATES                                  \n');
    fprintf(fid, '================================================================================\n');
    
    num_points      = length(x_complete);
    num_wall_points = length(x_wall_scaled);
    num_conv_pts    = max(num_points - num_wall_points, 0);  % 0 if no converging part

    if conv_option ~= 5 && ~isempty(x_complete)
      fprintf(fid, 'COMPLETE NOZZLE PROFILE (Including Converging Section):\n\n');
      fprintf(fid, ' Point        X              Y          Section\n');
      fprintf(fid, '-------  ------------  ------------  -------------\n');
    
       for pt = 1:num_points
            if pt <= num_conv_pts
                section = 'Converging';
            else
                section = 'Diverging';
            end
           fprintf(fid, '%7d  %12.6f  %12.6f  %s\n', pt, x_complete(pt), y_complete(pt), section);
     end
    else
      fprintf(fid, 'DIVERGING SECTION WALL COORDINATES:\n\n');
      fprintf(fid, ' Point        X              Y          Node Index\n');
      fprintf(fid, '-------  ------------  ------------  ------------\n');
      for i = 1:num_wall_points
          fprintf(fid, '%7d  %12.6f  %12.6f  %12d\n', i, x_wall_scaled(i), y_wall_scaled(i), wall_node_indices(i));
      end
    end

    fprintf(fid, '\n');

    
    % Performance Analysis
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '                           PERFORMANCE ANALYSIS                                \n');
    fprintf(fid, '================================================================================\n');
    if Me > 1
        Cf_ideal = sqrt(2*gamma^2/(gamma-1) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * ...
                       (1 - (1/Me^2)^((gamma-1)/gamma)));
        fprintf(fid, '%-40s : %.6f\n', 'Ideal Thrust Coefficient (Cf)', Cf_ideal);
        c_star = sqrt(gamma) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
        fprintf(fid, '%-40s : %.6f\n', 'Characteristic Velocity Parameter', c_star);
        m_dot_param = area_ratio / Me * sqrt(gamma) * (1 + (gamma-1)/2 * Me^2)^((gamma+1)/(2*(gamma-1)));
        fprintf(fid, '%-40s : %.6f\n', 'Mass Flow Parameter', m_dot_param);
        V_exit_ratio = Me * sqrt(2/(gamma+1) / (1 + (gamma-1)/2 * Me^2));
        fprintf(fid, '%-40s : %.6f\n', 'Exit Velocity Ratio (V/V*)', V_exit_ratio);
        P_exit_ratio = (1 + (gamma-1)/2 * Me^2)^(-gamma/(gamma-1));
        T_exit_ratio = (1 + (gamma-1)/2 * Me^2)^(-1);
        fprintf(fid, '%-40s : %.6f\n', 'Exit Pressure Ratio (P/P0)', P_exit_ratio);
        fprintf(fid, '%-40s : %.6f\n', 'Exit Temperature Ratio (T/T0)', T_exit_ratio);
    end
    fprintf(fid, '\n');
    
    % Grid Statistics
    fprintf(fid, 'GRID STATISTICS:\n');
    fprintf(fid, '%-40s : %.6f units\n', 'Minimum X Coordinate', min(x_scaled));
    fprintf(fid, '%-40s : %.6f units\n', 'Maximum X Coordinate', max(x_scaled));
    fprintf(fid, '%-40s : %.6f units\n', 'Minimum Y Coordinate', min(y_scaled));
    fprintf(fid, '%-40s : %.6f units\n', 'Maximum Y Coordinate', max(y_scaled));
    fprintf(fid, '%-40s : %.6f\n', 'Minimum Mach Number', min(M));
    fprintf(fid, '%-40s : %.6f\n', 'Maximum Mach Number', max(M));
    fprintf(fid, '\n');
    
    % Footer
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '                            END OF DETAILED REPORT                             \n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'Report generated by: Method of Characteristics Nozzle Design Code\n');
    fprintf(fid, 'Timestamp: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '********************************************************************************\n');
    
    fclose(fid);
end

function generateEnhancedSummary(filename, gamma, Me, R_throat, R_exit, area_ratio, ...
                                L_nozzle, L_total, D_throat, D_exit, A_throat, A_exit, ...
                                contraction_ratio, conv_option, n, node)
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot create summary file: %s', filename);
    end
    
    fprintf(fid, '╔══════════════════════════════════════════════════════════════════════════════╗\n');
    fprintf(fid, '║                     NOZZLE DESIGN EXECUTIVE SUMMARY                         ║\n');
    fprintf(fid, '╚══════════════════════════════════════════════════════════════════════════════╝\n\n');
    
    fprintf(fid, 'DESIGN OBJECTIVE:\n');
    fprintf(fid, '─────────────────\n');
    fprintf(fid, 'Minimum length supersonic nozzle for Mach %.2f uniform exit flow\n\n', Me);
    
    fprintf(fid, 'KEY SPECIFICATIONS:\n');
    fprintf(fid, '───────────────────\n');
    fprintf(fid, '  • Gamma (γ)              = %.3f\n', gamma);
    fprintf(fid, '  • Exit Mach Number       = %.2f\n', Me);
    fprintf(fid, '  • Throat Radius          = %.3f units\n', R_throat);
    fprintf(fid, '  • Exit Radius            = %.3f units\n', R_exit);
    fprintf(fid, '  • Area Ratio (Ae/At)     = %.3f\n', area_ratio);
    
    if conv_option ~= 5
        fprintf(fid, '  • Total Length           = %.3f units\n', L_total);
        fprintf(fid, '  • Contraction Ratio      = %.3f\n', contraction_ratio);
    else
        fprintf(fid, '  • Nozzle Length          = %.3f units\n', L_nozzle);
    end
    fprintf(fid, '\n');
    
    fprintf(fid, 'GEOMETRY:\n');
    fprintf(fid, '─────────\n');
    fprintf(fid, '  • Throat Diameter        = %.3f units\n', D_throat);
    fprintf(fid, '  • Throat Area            = %.3f units²\n', A_throat);
    fprintf(fid, '  • Exit Diameter          = %.3f units\n', D_exit);
    fprintf(fid, '  • Exit Area              = %.3f units²\n', A_exit);
    fprintf(fid, '  • L/D Ratio              = %.3f\n', L_nozzle/D_throat);
    fprintf(fid, '\n');
    
    fprintf(fid, 'COMPUTATIONAL:\n');
    fprintf(fid, '──────────────\n');
    fprintf(fid, '  • Characteristic Lines   = %d\n', n);
    fprintf(fid, '  • Total Nodes            = %d\n', node);
    fprintf(fid, '  • Wall Points            = %d\n', n+1);
    fprintf(fid, '\n');
    
    if conv_option ~= 5
        profile_names = {'Circular Arc', 'Witoszynski', 'Bicubic', 'Conical'};
        fprintf(fid, 'CONVERGING SECTION:\n');
        fprintf(fid, '───────────────────\n');
        fprintf(fid, '  • Profile Type           = %s\n', profile_names{conv_option});
        fprintf(fid, '  • Length                 = %.3f units\n', L_total - L_nozzle);
        fprintf(fid, '\n');
    end
    
    fprintf(fid, 'PERFORMANCE ESTIMATES:\n');
    fprintf(fid, '──────────────────────\n');
    if Me > 1
        Cf_ideal = sqrt(2*gamma^2/(gamma-1) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * ...
                       (1 - (1/Me^2)^((gamma-1)/gamma)));
        fprintf(fid, '  • Ideal Thrust Coeff.    = %.3f\n', Cf_ideal);
        fprintf(fid, '  • Exit/Stagnation P      = %.4f\n', (1 + (gamma-1)/2 * Me^2)^(-gamma/(gamma-1)));
        fprintf(fid, '  • Exit/Stagnation T      = %.4f\n', (1 + (gamma-1)/2 * Me^2)^(-1));
    end
    fprintf(fid, '\n');
    
    fprintf(fid, '─────────────────────────────────────────────────────────────────────────────\n');
    fprintf(fid, 'Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '═════════════════════════════════════════════════════════════════════════════\n');
    
    fclose(fid);
end

function savePerformanceMetrics(filename, gamma, Me, area_ratio, L_nozzle, ...
                               R_throat, R_exit, theta_max, n, node)
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot create performance file: %s', filename);
    end
    
    % Write header
    fprintf(fid, '# Nozzle Performance Metrics\n');
    fprintf(fid, '# Generated: %s\n', datestr(now));
    fprintf(fid, 'Metric,Value,Units,Description\n');
    
    % Basic performance metrics
    fprintf(fid, 'Exit Mach Number,%.6f,-,Design Mach number at nozzle exit\n', Me);
    fprintf(fid, 'Area Ratio,%.6f,-,Exit to throat area ratio\n', area_ratio);
    fprintf(fid, 'Length to Throat Radius,%.6f,-,Non-dimensional length\n', L_nozzle/R_throat);
    fprintf(fid, 'Maximum Flow Angle,%.6f,degrees,Maximum expansion angle\n', theta_max);
    
    if Me > 1
        P_ratio = (1 + (gamma-1)/2 * Me^2)^(-gamma/(gamma-1));
        fprintf(fid, 'Exit Pressure Ratio,%.6f,-,Pe/P0 (isentropic)\n', P_ratio);
        T_ratio = (1 + (gamma-1)/2 * Me^2)^(-1);
        fprintf(fid, 'Exit Temperature Ratio,%.6f,-,Te/T0 (isentropic)\n', T_ratio);
        rho_ratio = (1 + (gamma-1)/2 * Me^2)^(-1/(gamma-1));
        fprintf(fid, 'Exit Density Ratio,%.6f,-,ρe/ρ0 (isentropic)\n', rho_ratio);
        V_ratio = Me * sqrt(T_ratio);
        fprintf(fid, 'Exit Velocity Ratio,%.6f,-,Ve/a0 (relative to stagnation sound speed)\n', V_ratio);
        Cf_ideal = sqrt(2*gamma^2/(gamma-1) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * ...
                       (1 - (1/Me^2)^((gamma-1)/gamma)));
        fprintf(fid, 'Ideal Thrust Coefficient,%.6f,-,Cf for ideal expansion\n', Cf_ideal);
    end
    
    fprintf(fid, 'Divergence Efficiency,%.6f,-,Measure of expansion uniformity (ideal=1)\n', 1.0);
    fprintf(fid, 'Grid Resolution,%.6f,nodes/unit²,Computational mesh density\n', node/max(L_nozzle*max(R_exit,eps), eps));
    
    % Computational metrics
    fprintf(fid, 'Number of Characteristics,%d,-,MOC resolution\n', n);
    fprintf(fid, 'Total Nodes,%d,-,Computational grid points\n', node);
    fprintf(fid, 'Wall Points,%d,-,Nozzle contour points\n', n+1);
    fprintf(fid, 'Angular Resolution,%.6f,degrees,Flow angle increment\n', theta_max/n);
    
    fclose(fid);
end

%% Function Definitions

function [v, KL, KR, theta] = moc2d(theta_max, theta_0, n)
    dtheta = (theta_max - theta_0) / (n - 1);
    % 7+1,6+1,5+1,4+1,3+1,2+1,1+1
    node = 0.5 * n * (4 + n - 1);
    theta = zeros(1, node);
    v = zeros(1, node);
    KL = zeros(1, node);
    KR = zeros(1, node);
    
    for i = 1:n
        theta(i) = theta_0 + (i-1)*dtheta;
        v(i) = theta(i);
        KL(i) = theta(i) - v(i);
        KR(i) = theta(i) + v(i);
    end
    
    i = n + 1;
    theta(i) = theta(i-1);
    v(i) = v(i-1);
    KL(i) = KL(i-1);
    KR(i) = KR(i-1);
    p = 2;
    q = n + 2;
    
    for k = 1:n-1
        j = p;
        h = q;
        theta(h) = 0;
        KR(h) = KR(j);
        v(h) = KR(j) - theta(h);
        KL(h) = theta(h) - v(h);
        j = j + 1;
        
        for i = h+1:n-p+q
            KR(i) = KR(j);
            KL(i) = KL(i-1);
            theta(i) = 0.5*(KL(i) + KR(i));
            v(i) = 0.5*(KR(i) - KL(i));
            j = j + 1;
        end
        
        if i == n-p+q
            h = i + 1;
        else
            h = h + 1;
        end
        
        theta(h) = theta(h-1);
        v(h) = v(h-1);
        KL(h) = KL(h-1);
        KR(h) = KR(h-1);
        p = p + 1;
        q = h + 1;
    end
end

function mu = Mu(M)
    mu = asind(1/M);
end

function v = PrandtlMeyer(M, gamma)
    v = sqrt((gamma+1)/(gamma-1)) * atan(sqrt((gamma-1)/(gamma+1)*(M^2-1))) - atan(sqrt(M^2-1));
end
