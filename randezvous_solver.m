clc; clear; clear all;

%% Constants
pi = 4*atan(1);
R_Earth = 6378.137; % [km]
mu_Earth = 3.986004418e5; % [km^3/s^2]
lower_orb_altitude = 205; % [km]
R_low = R_Earth + lower_orb_altitude;
theta = 0:1:360; % true anomaly
a_ISS = 6790.69; % semi-major axis [ISS]
e_ISS = 0.0008389; % eccentricity [ISS]
b_ISS = a_ISS * sqrt(1-e_ISS^2); % semi-minor axis [ISS]
r_p_ISS = a_ISS * (1 + e_ISS); % perigee radius
r_a_ISS = a_ISS * (1 - e_ISS); % apogee radius
r_p_transfer = R_low; % perigee radius

%% Matrix initials
x_low = ones(1,length(theta));
y_low = ones(1,length(theta));
x_Earth = ones(1,length(theta));
y_Earth = ones(1,length(theta));
x_ISS = ones(1,length(theta));
y_ISS = ones(1,length(theta));
r_ISS = ones(1,length(theta));

%% Lower orbit
for i = 1:length(theta)    
    x_low(i) = R_low * cosd(theta(i));
    y_low(i) = R_low * sind(theta(i));
end

%% Earth
for j = 1:length(theta)    
    x_Earth(j) = R_Earth * cosd(theta(j));
    y_Earth(j) = R_Earth * sind(theta(j));
end

%% ISS orbit
for k = 1:length(theta)
    r_ISS(k) = (a_ISS * (1 - e_ISS^2)) / (1 + e_ISS * cosd(theta(k)));
    x_ISS(k) = r_ISS(k) .* cosd(theta(k));
    y_ISS(k) = r_ISS(k) .* sind(theta(k));
end  
T_ISS = 2*pi*sqrt((a_ISS^3)/mu_Earth)/60;

%% Transfer orbit
random_value1 = randi([1, length(theta)]);
x1 = x_low(random_value1); % Transfer orbit entry
y1 = y_low(random_value1); % Transfer orbit entry
theta1 = theta(random_value1);
if (theta1 < 180)
    theta2 = theta1 + 180;
else
    theta2 = theta1 - 180;
end
index1 = find(theta == theta2);
x2 = x_ISS(index1); % Transfer orbit exit
y2 = y_ISS(index1); % Transfer orbit exit
r_a_transfer = r_ISS(index1); % apogee radius
a_transfer = (r_a_transfer + r_p_transfer ) / 2; % semi-major axis [transfer]
T_transfer = 2*pi*sqrt((a_transfer^3)/mu_Earth)/60;

%% ISS orbit future point location
random_value2 = randi([1, length(theta)]);
x3 = x_ISS(random_value2);
y3 = y_ISS(random_value2);
theta3 = theta(random_value2);
E_anomaly_ISS = (e_ISS+cosd(theta3))/(1+e_ISS*cosd(theta3));
E_anomaly_ISS = acosd(E_anomaly_ISS); % Eccentric anomaly of ISS
if (theta3 > 180)
    E_anomaly_ISS = 360 - E_anomaly_ISS;       
end
E_anomaly_ISS_radian = E_anomaly_ISS*pi/180; 
t1 = T_ISS * (E_anomaly_ISS_radian-e_ISS*sin(E_anomaly_ISS_radian)) / (2*pi);
t2 = T_transfer/2;
t3 = t1 + t2;

M_ISS = 2*pi*t3/T_ISS;

E_future_initial = 1;
E_future_ISS = E_future_initial;
tol = 1e-10;                   % Tolerance
max_iter = 1000;              % Maximum iterations
i = 0; % # of iterations

while true
    % Newton-Raphson
    f = E_future_ISS - e_ISS * sin(E_future_ISS) - M_ISS;  % Kepler's equation
    f_prime = 1 - e_ISS * cos(E_future_ISS);               % Derivative
    E_next = E_future_ISS - f / f_prime;

    if abs(E_next - E_future_ISS) < tol
        break;
    end
    if i > max_iter
        error('Newton-Raphson failed to converge');
    end

    E_future_ISS = E_next;
    i = i + 1;
    
end

E_future_ISS = E_future_ISS*180/pi;
theta4 = (cosd(E_future_ISS) - e_ISS) / (1 - e_ISS*cosd(E_future_ISS));
theta4 = round(acosd(theta4));
index2 = find(theta == theta4);
x4 = x_ISS(index2); % future ISS coordinate
y4 = y_ISS(index2); % future ISS coordinate

%% Plotting
figure;
plot(x_low, y_low, 'g', LineWidth=1.5)
hold on;
plot(x_ISS, y_ISS, 'b', LineWidth=1.5)
earthBlue = [0, 0.2, 0.6]; % Earth's blue color
fill(x_Earth, y_Earth, earthBlue, 'EdgeColor', 'black');
axis equal;
ylim([-9000 9000])
xlim([-9000 9000])
xlabel('x-axis [km]')
ylabel('y-axis [km]')
set(gca, 'Color', 'k');
set(gca, 'FontSize', 24);

tran_entry = plot(x1,y1, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'c'); 
tran_exit = plot(x2,y2, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
initial_ISS = plot(x3,y3, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'm'); 
fut_ISS = plot(x4,y4, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
legend([tran_entry, tran_exit, initial_ISS, fut_ISS], 'Transfer orbit entry',...
    'Transfer orbit exit', 'ISS initial coordinate', 'ISS future coordinate', ...
    'Location', 'northeast', 'TextColor', [1, 1, 0], 'FontSize', 24);

%% Time vs Theta for ISS
% for p=1:length(theta)
% E_anomaly_ISS(p) = (e_ISS+cosd(theta(p)))/(1+e_ISS*cosd(theta(p)));
% E_anomaly_ISS(p) = acosd(E_anomaly_ISS(p)); % Eccentric anomaly of ISS
% if (theta(p) > 180)
%     E_anomaly_ISS(p) = 360 - E_anomaly_ISS(p);       
% end
% E_anomaly_ISS_radian(p) = E_anomaly_ISS(p)*pi/180; 
% t1(p) = T_ISS * (E_anomaly_ISS_radian(p)-e_ISS*sin(E_anomaly_ISS_radian(p))) / (2*pi);
% end
% 
% plot(t1,theta)

