clear; close all; clc;

% DETUMBLING SECTION

% Time parameters
dt = 1;
duration_det = 1200;
time_steps = duration_det/dt;

% Orbital parameters
a = 6740e3; 
e = 0.0003;
i = deg2rad(30);
omega_periapsis = deg2rad(0);
omega_ascNode = deg2rad(0);
mu = 3.986004418*10^14; 
Period = 2*pi*sqrt((a^3)/mu); 

% Inertial parameter
Inertial = diag([5, 5, 5]);

wgs84 = wgs84Ellipsoid('kilometer');

% Time parameters 
currYear = 2024;
currMonth = 11;
currDay = 15;
currHour = 0;
currMin = 0;
currSec = 0;

% Initial angular velocity
angular_velocity = [2; 4; -5]; 
angular_velocity_data = zeros(3, time_steps);

% Initialize quaternion for attitude (starting at identity rotation)
q = [1; 0; 0; 0];  
q_history = zeros(4, time_steps);
q_history(:,1) = q;

% Initialize body frame axes
body_x = [1; 0; 0];
body_y = [0; 1; 0];
body_z = [0; 0; 1];

% Body axes history
body_x_history = zeros(3, time_steps);
body_y_history = zeros(3, time_steps);
body_z_history = zeros(3, time_steps);

body_x_history(:,1) = body_x;
body_y_history(:,1) = body_y;
body_z_history(:,1) = body_z;

for n = 1:time_steps
    % Original orbit and magnetic field calculations
    theta = TrueAnomaly(n, Period);
    r = a * (1 - e^2) / (1 + e * cos(theta));
    x = r * (cos(omega_periapsis + theta) * cos(omega_ascNode) - sin(omega_periapsis + theta) * cos(i) * sin(omega_ascNode));
    y = r * (cos(omega_periapsis + theta) * sin(omega_ascNode) + sin(omega_periapsis + theta) * cos(i) * cos(omega_ascNode));
    z = r * sin(omega_periapsis + theta) * sin(i);

    % Covert longitude, latitude, altitude to Bdot information
    [lat, long, alt] = ecef2geodetic(wgs84, x, y, z);
    [NEDx, NEDy, NEDz] = igrf('15-Nov-2024', lat, long, alt, 'geodetic'); 
    [ECEFx, ECEFy, ECEFz] = ned2ecef(NEDx, NEDy, NEDz, lat, long, alt, wgs84, "radians");
    
    utc = [currYear currMonth currDay currHour currMin currSec];
    ECEFPos = [ECEFx; ECEFy; ECEFz];
    B_eci = ecef2eci(utc, ECEFPos);
    b = B_eci / norm(B_eci);
    
    k = 0.1; % Changing this changes end direction and how long it takes to stabilize
    Torque = cross(k*(cross(angular_velocity, b)), b);
    
    % Angular velocity update
    angular_acceleration = Inertial \ Torque;
    angular_velocity = angular_velocity + angular_acceleration * dt;
    angular_velocity_data(:,n) = angular_velocity;
    
    % Update quaternion using angular velocity
    omega_norm = norm(angular_velocity);
    if omega_norm > 0
        axis = angular_velocity / omega_norm;
        angle = omega_norm * dt;
        
        % Rotation quaternion
        dq = [cos(angle/2); 
              axis(1)*sin(angle/2);
              axis(2)*sin(angle/2);
              axis(3)*sin(angle/2)];
        
        % Quaternion multiplication
        q = quatmultiply(q', dq')';
        q = q / norm(q);  % Normalize 
    end
    
    % Store quaternion history
    q_history(:,n) = q;
    
    % Update body frame vectors using current quaternion
    R = quat2rotm(q');  % Convert quaternion to rotation matrix
    
    % Rotate body frame vectors
    body_x = R * [1; 0; 0];
    body_y = R * [0; 1; 0];
    body_z = R * [0; 0; 1];
    
    % Store body frame vectors
    body_x_history(:,n) = body_x;
    body_y_history(:,n) = body_y;
    body_z_history(:,n) = body_z;
end

% Plotting
time_vector_det = (1:time_steps);

% Original angular velocity plot
figure(1);
set(gcf, 'Name', 'Detumbling Plots', 'NumberTitle', 'off');
set(gcf, 'Position', [50, 150, 700, 600]); % Position: [x y width height]

% Plot 1: Angular Velocity
subplot(2,1,1);
plot(time_vector_det, angular_velocity_data(1,:), 'b-', ...
    time_vector_det, angular_velocity_data(2,:), 'r-', ...
    time_vector_det, angular_velocity_data(3,:), 'g-');
ylabel('Angular Velocity Vector Components (rad/s)');
xlabel('Time (minutes)')
title('Angular Velocity Vector Components');
legend('\omega_x', '\omega_y', '\omega_z');
grid on;

% Plot 2: X, Y, Z directions 
subplot(2,1,2);
hold on;
plot(time_vector_det, body_x_history(3,:), 'b-', ...
     time_vector_det, body_y_history(3,:), 'r-', ...
     time_vector_det, body_z_history(3,:), 'g-');
ylabel('Direction Vector Components');
xlabel('Time (minutes)');
legend('X', 'Y', 'Z');
title('Satellite Direction Vector Components');
grid on;
hold off;


%% 

% SUN ACQUISITION SECTION

% From Detumbling Section
final_x = body_x_history(3,end);
final_y = body_y_history(3,end);
final_z = body_z_history(3,end);

A = [final_x, final_y, final_z];  % Initial direction vector

% Adjust accordingly
Kp = 0.15;  
Ki = 0.005;
Kd = 0.01; 

integral_error = 0;
previous_error = 0;
dt_sa = 0.1; 
duration_sa = 600; 

% History arrays - all vectors are in ECI frame
theta_history = zeros(1, duration_sa);
A_history = zeros(3, duration_sa);  % Actual direction vector
desired_history = zeros(3, duration_sa);  % Desired Sun vector
sun_pos_history = zeros(3, duration_sa);  % Sun position vector 
satellite_pos_history = zeros(3, duration_sa);  % Satellite position vector 


for k = 1:duration_sa
    % Calculate position
    th = TrueAnomaly(k, Period);
    r = a * (1 - e^2) / (1 + e * cos(th));
    x = r * (cos(omega_periapsis + th) * cos(omega_ascNode) - sin(omega_periapsis + th) * cos(i) * sin(omega_ascNode));
    y = r * (cos(omega_periapsis + th) * sin(omega_ascNode) + sin(omega_periapsis + th) * cos(i) * cos(omega_ascNode));
    z = r * sin(omega_periapsis + th) * sin(i);
    
    % Satellite position vector 
    EarthSatelliteVector = [x; y; z];
    satellite_pos_history(:,k) = EarthSatelliteVector / norm(EarthSatelliteVector);
    
    % Sun position vector
    B = CalcSunVector(currYear, currMonth, currDay, currHour, currMin, currSec);
    sun_pos_history(:,k) = B / norm(B); % Normalized
    
    % Desired Sun vector 
    SatelliteSunVector = B - EarthSatelliteVector;
    desired_direction = SatelliteSunVector / norm(SatelliteSunVector);
    
    % Store current directions
    A_history(:, k) = A / norm(A);  % Normalized 
    desired_history(:, k) = desired_direction;

    % Angle between current and desired vector
    costheta = dot(A, SatelliteSunVector) / (norm(A) * norm(SatelliteSunVector));
    theta = acos(costheta);
    theta_history(:, k) = theta;
    
    % Prevents fluctuating after low error is achieved
    if abs(theta) < 0.02
        integral_error = 0;
    end

    % PID Control and Torque
    integral_error = integral_error + theta * dt_sa;
    derivative_error = (theta - previous_error) / dt_sa;
    
    Torque = Kp * theta + Ki * integral_error + Kd * derivative_error;
    max_torque = 0.5; % Can adjust 
    Torque = min(max(Torque, -max_torque), max_torque);

    rotation_axis = cross(A, SatelliteSunVector);  
    if norm(rotation_axis) > 0
        rotation_axis = rotation_axis / norm(rotation_axis);
    end
    
    q = [cos(Torque/2), rotation_axis * sin(Torque/2)];  % Rotation axis in quaternion form, using Torque
    A_quat = [0, A]; % Current direction vector in quaternion form
    q_inv = [q(1), -q(2), -q(3), -q(4)]; % For quaternion multiplication

    % Calculate and update new direction vector after torque rotation is applied
    A_new = quatmultiply(quatmultiply(q, A_quat), q_inv); 
    A = A_new(2:4); 
    
    previous_error = theta;

    currSec = currSec + dt_sa * 60;  
    if currSec >= 60
        currSec = currSec - 60;
        currMin = currMin + 1;
        if currMin == 60
            currMin = 0;
            currHour = currHour + 1;
            if currHour == 24
                currHour = 0;
                currDay = currDay + 1;
            end
        end
    end
end

% Plotting 
time_vector_sa = (1:duration_sa);

figure(2);
set(gcf, 'Name', 'Sun Acquisition Plots', 'NumberTitle', 'off');
set(gcf, 'Position', [750, 150, 700, 600]); % Position: [x y width height]

% Plot 1: Angular Error
subplot(3,1,1);
plot(time_vector_sa, theta_history, 'o', 'MarkerFaceColor', 'k', 'Color', 'k');
xlabel('Time (minutes)');
ylabel('Angle Error (radians)');
title('Angle Error Over Time');
grid on;

% Plot 2: Position Vectors
subplot(3,1,2);
hold on;
plot(time_vector_sa, sun_pos_history(1,:), 'b-', 'DisplayName', 'Sun X');
plot(time_vector_sa, sun_pos_history(2,:), 'r-', 'DisplayName', 'Sun Y');
plot(time_vector_sa, sun_pos_history(3,:), 'g-', 'DisplayName', 'Sun Z');
plot(time_vector_sa, satellite_pos_history(1,:), 'b--', 'DisplayName', 'Sat X');
plot(time_vector_sa, satellite_pos_history(2,:), 'r--', 'DisplayName', 'Sat Y');
plot(time_vector_sa, satellite_pos_history(3,:), 'g--', 'DisplayName', 'Sat Z');
xlabel('Time (minutes)');
ylabel('Position Vector Components');
title('Earth-Sun vs Earth-Satellite Position Vector Components');
legend('show');
grid on;

% Plot 3: Direction Vectors
subplot(3,1,3);
hold on;
plot(time_vector_sa, A_history(1,:), 'b-', 'DisplayName', 'Actual X');
plot(time_vector_sa, desired_history(1,:), 'b--', 'DisplayName', 'Desired X');
plot(time_vector_sa, A_history(2,:), 'r-', 'DisplayName', 'Actual Y');
plot(time_vector_sa, desired_history(2,:), 'r--', 'DisplayName', 'Desired Y');
plot(time_vector_sa, A_history(3,:), 'g-', 'DisplayName', 'Actual Z');
plot(time_vector_sa, desired_history(3,:), 'g--', 'DisplayName', 'Desired Z');
xlabel('Time (minutes)');
ylabel('Direction Vector Components');
title('Current Satellite Direction vs Desired Sun Direction Vector Components');
legend('show');
grid on;

%% 

% FUNCTIONS SECTION

% Function to calculate EarthSun Vector
function sv = CalcSunVector(y, mon, d, h, min, s)
    JD = 367*y - floor(7*(y + floor(mon + 9) / 12)) + floor(275 * mon / 9) + d + 1721013.5 + ((((s/60) + min)/60) + h)/24;
    T = (JD - 2451545.0) / 36525;
    lambdaMo = 280.2606184 + 36000.77005361 * T;
    Mo = 357.5277233 + 359.05034 * T;
    lambdaEcliptical = lambdaMo + 1.1914666471 * sin(Mo) + 0.019994643 * sin(Mo);
    epsilon = 23.439291 - 0.0130042 * T;
    sv = [cos(lambdaEcliptical); cos(epsilon) * sin(lambdaEcliptical); sin(epsilon) * sin(lambdaEcliptical)];
end

% Function to calculate true anomaly from time and period
function theta = TrueAnomaly(t, period)
    % t is in minutes, period is in seconds
    t_seconds = t * 60;  % Convert to seconds
    meanMotion = 2*pi/period;  % radians per second
    theta = meanMotion * t_seconds;  % Returns angle in radians
end


