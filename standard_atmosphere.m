
%% LVD I --- HW4 --- Nicolas Sarmiento

%% Standard Atmosphere Code

% The following code will perform the calculations for properties like
% temperature, pressure, density, and viscosity at different altitudes
% in the atmosphere.

% Note that this file will export the tabular numbers as a csv file under
% the name 'Standard_Atmosphere_Nicolas_Sarmiento.csv' and will also export
% on a separate pdf file the corresponding plots under the name
% 'Standard_Atmosphere_Nicolas_Sarmiento.pdf' 

%% Definition and Initiation of Arrays, Constants, and Variables

% Altitude Array (km)
altitude = 0:1:86;

% Preallocate Result Arrays
viscosity = zeros(size(altitude));
pressure = zeros(size(altitude));
density = zeros(size(altitude));
temperature = zeros(size(altitude));

% Constant Values
g = 9.80665; % m/s^2
R_air = 287.05; % J/(Kg-K)

% ------------------------------------------------------------------ %
% Troposphere, Tropopause, Stratosphere, Stratosphere 2, Stratopause,
% Mesosphere, Mesosphere 2, Mesopause 
% ------------------------------------------------------------------ %

% Layer Altitude Array (m)
layer_boundary_altitude = [0, 11000, 20000, 32000, 47000, 51000, 71000, 86000]; % m

% Temperature Lapse Rate Array (ºC/m)
temp_lapse_rate = [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002];

% Density, Temperature, and Pressure Base Values for Each Layer
dens_base_values = [1.225, 0.3639, 0.088, 0.0132, 0.0014, 0.0009, 0.0001, 0]; % kg/m^3
pres_base_values = [101325, 22632.1, 5474.89, 868.02, 110.91, 66.94, 3.96]; % Pa
temp_base_values = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65]; % K

%% Data Calculations

% Interation Through Altitude Array
for i = 1:length(altitude)
    h = altitude(i) * 1000;

    % If Statement to Determine Location of Altitude Within Atmosphere Regions
    if h >= layer_boundary_altitude(end)
        layer = length(layer_boundary_altitude) - 1;
    else
        layer = find(h >= layer_boundary_altitude, 1, 'last');
    end

    % Pressure and Temperature Calculations
    if temp_lapse_rate(layer) == 0
        % Altitude lies within an isothermic layer
        temperature(i) = temp_base_values(layer);
        pressure(i) = pres_base_values(layer) * exp(-g * (h - layer_boundary_altitude(layer)) / (R_air * temp_base_values(layer)));
    else
        % Altitude lies within a gradient layer
        temperature(i) = temp_base_values(layer) + temp_lapse_rate(layer) * (h - layer_boundary_altitude(layer));
        pressure(i) = pres_base_values(layer) * (temperature(i) / temp_base_values(layer))^(-g / (R_air * temp_lapse_rate(layer)));
    end

    % Density Calculation
    density(i) = pressure(i) / (R_air * temperature(i));

    % Dynamic Viscosity Calculation
    C1 = 1.458e-6; % kg/(m-s-K^(1/2))
    S = 110.4; % K
    viscosity(i) = C1 * temperature(i)^(3/2) / (temperature(i) + S);
end

%% Organize Results into a Table and Export as a csv File

output = table(altitude', temperature', pressure', density', viscosity', ...
    'VariableNames', {'Altitude', 'Temperature', 'Pressure', 'Density', 'Viscosity'});

[filename, pathname] = uiputfile('Standard_Atmosphere_Nicolas_Sarmiento.csv', 'Save as');

if isequal(filename, 0) || isequal(pathname, 0)
    disp(" ");
    disp('User canceled file saving');
    disp(" ");
else
    fullFilePath = fullfile(pathname, filename);
    writetable(output, fullFilePath);

    disp(" ");
    fprintf('File saved successfully to : %s\n', fullFilePath);
    disp(" ");
end

%% Plot Graphs and Export to a pdf File

[filename, pathname] = uiputfile('Standard_Atmosphere_Nicolas_Sarmiento.pdf', 'Save as');

if isequal(filename, 0) || isequal(pathname, 0)
    disp(" ");
    disp('User canceled file saving');
    disp(" ");
else
    pdf_file_name = fullfile(pathname, filename);
end

% Temperature vs Altitude
figure;
plot(temperature, altitude, 'r', 'LineWidth', 2);
ylabel('Altitude (km)');
xlabel('Temperature (K)');
title('Temperature vs Altitude Graph');
grid on;
exportgraphics(gcf, pdf_file_name, 'Append', true);

% Pressure vs Altitude
figure;
plot(pressure, altitude, 'r', 'LineWidth', 2);
ylabel('Altitude (km)');
xlabel('Pressure (Pa)');
title('Pressure vs Altitude Graph');
grid on;
exportgraphics(gcf, pdf_file_name, 'Append', true);

% Density vs Altitude
figure;
plot(density, altitude, 'r', 'LineWidth', 2);
ylabel('Altitude (km)');
xlabel('Density (kg/m^3)');
title('Density vs Altitude Graph');
grid on;
exportgraphics(gcf, pdf_file_name, 'Append', true);

% Viscosity vs Altitude
figure;
plot(viscosity, altitude, 'r', 'LineWidth', 2);
ylabel('Altitude (km)');
xlabel('Viscosity (kg/(m-s))');
title('Viscosity vs Altitude Graph');
grid on;
exportgraphics(gcf, pdf_file_name, 'Append', true);

disp(" ");
fprintf('File saved successfully to: %s\n', pdf_file_name);
disp(" ");   