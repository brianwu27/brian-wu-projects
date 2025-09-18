
a = 10000; % Semi-major axis
e = 0.5;   % Eccentricity
i = 45;    % Inclination in degrees
omega_periapsis = 45; % Argument of periapsis in degrees
omega_ascNode = 45; % Longitude of ascending node in degrees

i = deg2rad(i);
omega_periapsis = deg2rad(omega_periapsis);
omega_ascNode = deg2rad(omega_ascNode);

mu = 1.137*10^11; % km^3/s^2

theta = linspace(0, 2*pi, 360);

x = zeros(1, length(theta));
y = zeros(1, length(theta));
z = zeros(1, length(theta));

for k = 1:length(theta)
    r = a * (1 - e^2) / (1 + e * cos(theta(k)));
    x(k) = r * (cos(omega_periapsis + theta(k)) * cos(omega_ascNode) - sin(omega_periapsis + theta(k)) * cos(i) * sin(omega_ascNode));
    y(k) = r * (cos(omega_periapsis + theta(k)) * sin(omega_ascNode) + sin(omega_periapsis + theta(k)) * cos(i) * cos(omega_ascNode));
    z(k) = r * sin(omega_periapsis + theta(k)) * sin(i);
end


figure;
plot3(x, y, z);
hold on; 
fill3(x, y, z, 'c', 'FaceAlpha', 0.5);
grid on; 
title('3D Plot of the Orbit with the Sun');
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');

[sx, sy, sz] = sphere(50); 
sunRadius = 3000; % Sun's radius km
surf(sunRadius * sx, sunRadius * sy, sunRadius * sz, 'FaceColor', 'y', 'EdgeColor', 'none'); 

% Reference direction line (x-axis)
refLength = 1.5 * max(x); % Length
plot3([-refLength, refLength], [0, 0], [0, 0], 'k--', 'LineWidth', 1.5); 

% Equatorial plane 
[xp, yp] = meshgrid(linspace(-refLength, refLength, 50));
zp = zeros(size(xp));
surf(xp, yp, zp, 'FaceColor', 'blue', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); 

% Specific true anomaly for satellite
theta_celestial = deg2rad(165); 

% Calculate position for satellite
r_celestial = a * (1 - e^2) / (1 + e * cos(theta_celestial));
x_celestial = r_celestial * (cos(omega_periapsis + theta_celestial) * cos(omega_ascNode) - sin(omega_periapsis + theta_celestial) * cos(i) * sin(omega_ascNode));
y_celestial = r_celestial * (cos(omega_periapsis + theta_celestial) * sin(omega_ascNode) + sin(omega_periapsis + theta_celestial) * cos(i) * cos(omega_ascNode));
z_celestial = r_celestial * sin(omega_periapsis + theta_celestial) * sin(i);

plot3(x_celestial, y_celestial, z_celestial, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Red dot

% Calculate position for periapsis
theta_periapsis = 0; % Periapsis occurs at true anomaly = 0
r_periapsis = a * (1 - e^2) / (1 + e * cos(theta_periapsis));
x_periapsis = r_periapsis * (cos(omega_periapsis + theta_periapsis) * cos(omega_ascNode) - sin(omega_periapsis + theta_periapsis) * cos(i) * sin(omega_ascNode));
y_periapsis = r_periapsis * (cos(omega_periapsis + theta_periapsis) * sin(omega_ascNode) + sin(omega_periapsis + theta_periapsis) * cos(i) * cos(omega_ascNode));
z_periapsis = r_periapsis * sin(omega_periapsis + theta_periapsis) * sin(i);

plot3(x_periapsis, y_periapsis, z_periapsis, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'black'); % Green dot


axis equal; 
view(3); 

hold off; 