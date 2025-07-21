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

%% Converging section design %%
fprintf('\n========== CONVERGING SECTION DESIGN ==========\n');
fprintf('Select converging section profile:\n');
fprintf('1. Circular arc\n');
fprintf('2. Witoszynski curve (bell-mouth)\n');
fprintf('3. Bicubic (smooth polynomial)\n');
fprintf('4. Conical (linear)\n');
fprintf('5. No converging section (sharp corner only)\n');
conv_option = input('Enter option (1-5): ');

if conv_option ~= 5
    R_inlet = input('Enter inlet radius (typically 2-4 times throat radius): ');
    L_conv = input('Enter converging section length: ');
    
    % Generate converging section coordinates
    n_conv = 50; % Number of points in converging section
    x_conv = linspace(-L_conv, 0, n_conv);
    y_conv = zeros(1, n_conv);
    
    switch conv_option
        case 1 % Circular arc
            % Calculate radius of curvature for smooth transition
            R_curve = ((L_conv)^2 + (R_inlet - R_throat)^2) / (2*(R_inlet - R_throat));
            y_center = R_inlet - R_curve;
            for i = 1:n_conv
                y_conv(i) = y_center + sqrt(R_curve^2 - (x_conv(i) + L_conv)^2);
            end
            
        case 2 % Witoszynski curve
            % Bell-mouth profile for minimum pressure loss
            for i = 1:n_conv
                xi = (x_conv(i) + L_conv) / L_conv; % Normalized position 0 to 1
                y_conv(i) = R_throat + (R_inlet - R_throat) * (1 - xi^2)^2;
            end
            
        case 3 % Bicubic polynomial
            % Smooth polynomial with zero slope at inlet
            for i = 1:n_conv
                xi = (x_conv(i) + L_conv) / L_conv; % Normalized position 0 to 1
                y_conv(i) = R_throat + (R_inlet - R_throat) * (1 - 3*xi^2 + 2*xi^3);
            end
            
        case 4 % Conical (linear)
            y_conv = linspace(R_inlet, R_throat, n_conv);
    end
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
xlabel('x (units)')
ylabel('y (units)')
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
num_display = min(10, node);
for i = 1:num_display
    fprintf('%-10d %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f\n', ...
            i, x_scaled(i), y_scaled(i), theta(i), M(i), mu(i), v(i));
end

% Display last 10 nodes if more than 20 nodes total
if node > 20
    fprintf('%-10s %-12s %-12s %-12s %-12s %-12s %-12s\n', '...', '...', '...', '...', '...', '...', '...');
    for i = max(node-9, num_display+1):node
        fprintf('%-10d %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f %-12.4f\n', ...
                i, x_scaled(i), y_scaled(i), theta(i), M(i), mu(i), v(i));
    end
end

% Save results to file
save_option = input('\nSave results to file? (y/n): ', 's');
if strcmpi(save_option, 'y')
    filename = input('Enter filename (without extension): ', 's');
    
    % Save workspace variables
    save([filename '.mat']);
    
    % Save nozzle coordinates to CSV (scaled values)
    if conv_option ~= 5
        nozzle_data = [x_complete', y_complete'];
        csvwrite([filename '_complete_coordinates.csv'], nozzle_data);
        
        % Also save separate files for converging and diverging sections
        conv_data = [x_conv', y_conv'];
        div_data = [x_wall_scaled', y_wall_scaled'];
        csvwrite([filename '_converging_coordinates.csv'], conv_data);
        csvwrite([filename '_diverging_coordinates.csv'], div_data);
    else
        nozzle_data = [x_wall_scaled', y_wall_scaled'];
        csvwrite([filename '_coordinates.csv'], nozzle_data);
    end
    
    % Save detailed report to text file
    fid = fopen([filename '_report.txt'], 'w');
    
    fprintf(fid, 'MINIMUM LENGTH NOZZLE DESIGN REPORT\n');
    fprintf(fid, 'Generated on: %s\n\n', datestr(now));
    
    fprintf(fid, 'DESIGN PARAMETERS:\n');
    fprintf(fid, 'Gamma = %.4f\n', gamma);
    fprintf(fid, 'Exit Mach Number = %.4f\n', Me);
    fprintf(fid, 'Throat Radius = %.4f units\n', R_throat);
    fprintf(fid, 'Number of Characteristics = %d\n', n);
    fprintf(fid, 'Maximum Flow Angle = %.4f degrees\n', theta_max);
    
    if conv_option ~= 5
        fprintf(fid, '\nCONVERGING SECTION:\n');
        switch conv_option
            case 1
                fprintf(fid, 'Profile Type: Circular Arc\n');
            case 2
                fprintf(fid, 'Profile Type: Witoszynski Curve (Bell-mouth)\n');
            case 3
                fprintf(fid, 'Profile Type: Bicubic Polynomial\n');
            case 4
                fprintf(fid, 'Profile Type: Conical (Linear)\n');
        end
        fprintf(fid, 'Inlet Radius = %.4f units\n', R_inlet);
        fprintf(fid, 'Converging Length = %.4f units\n', L_conv);
    end
    
    fprintf(fid, '\nNOZZLE FEATURES:\n');
    if conv_option ~= 5
        fprintf(fid, 'Total Length = %.4f units\n', L_total);
        fprintf(fid, 'Contraction Ratio = %.4f\n', contraction_ratio);
    else
        fprintf(fid, 'Nozzle Length = %.4f units\n', L_nozzle);
    end
    fprintf(fid, 'Throat Diameter = %.4f units\n', D_throat);
    fprintf(fid, 'Exit Diameter = %.4f units\n', D_exit);
    fprintf(fid, 'Area Ratio = %.4f\n', area_ratio);
    
    if conv_option ~= 5
        fprintf(fid, '\nCOMPLETE NOZZLE COORDINATES:\n');
        fprintf(fid, 'X\t\tY\n');
        for i = 1:length(x_complete)
            fprintf(fid, '%.6f\t%.6f\n', x_complete(i), y_complete(i));
        end
    else
        fprintf(fid, '\nNOZZLE WALL COORDINATES:\n');
        fprintf(fid, 'X\t\tY\n');
        for i = 1:length(x_wall_scaled)
            fprintf(fid, '%.6f\t%.6f\n', x_wall_scaled(i), y_wall_scaled(i));
        end
    end
    
    fclose(fid);
    fprintf('\nResults saved to:\n');
    fprintf('  - %s.mat (workspace variables)\n', filename);
    if conv_option ~= 5
        fprintf('  - %s_complete_coordinates.csv (full nozzle coordinates)\n', filename);
        fprintf('  - %s_converging_coordinates.csv (converging section only)\n', filename);
        fprintf('  - %s_diverging_coordinates.csv (diverging section only)\n', filename);
    else
        fprintf('  - %s_coordinates.csv (nozzle coordinates)\n', filename);
    end
    fprintf('  - %s_report.txt (detailed report)\n', filename);
end

toc

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
