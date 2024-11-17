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
meet_senarios = [];
x_low = ones(1,length(theta));
y_low = ones(1,length(theta));
x_Earth = ones(1,length(theta));
y_Earth = ones(1,length(theta));
x_ISS = ones(1,length(theta));
y_ISS = ones(1,length(theta));
r_ISS = ones(1,length(theta));

%% Lower orbit
for u = 1:length(theta)    
    x_low(u) = R_low * cosd(theta(u));
    y_low(u) = R_low * sind(theta(u));
end

%% Earth
for q = 1:length(theta)    
    x_Earth(q) = R_Earth * cosd(theta(q));
    y_Earth(q) = R_Earth * sind(theta(q));
end

%% ISS orbit
for k = 1:length(theta)
    r_ISS(k) = (a_ISS * (1 - e_ISS^2)) / (1 + e_ISS * cosd(theta(k)));
    x_ISS(k) = r_ISS(k) .* cosd(theta(k));
    y_ISS(k) = r_ISS(k) .* sind(theta(k));
end  
T_ISS = 2*pi*sqrt((a_ISS^3)/mu_Earth)/60;

%% The loop
for i = 1:length(theta)
    for j = 1:length(theta)
    
        % Transfer orbit
        x1(i,j) = x_low(j); % Transfer orbit entry
        y1(i,j) = y_low(j); % Transfer orbit entry
        theta1 = theta(j);
        if (theta1 < 180)
            theta2 = theta1 + 180;
        else
            theta2 = theta1 - 180;
        end
        index1 = find(theta == theta2);
        x2(i,j) = x_ISS(index1); % Transfer orbit exit
        y2(i,j) = y_ISS(index1); % Transfer orbit exit
        r_a_transfer = r_ISS(index1); % apogee radius
        a_transfer = (r_a_transfer + r_p_transfer ) / 2; % semi-major axis [transfer]
        T_transfer = 2*pi*sqrt((a_transfer^3)/mu_Earth)/60;
    
        % ISS orbit future point location
        
        x3(i,j) = x_ISS(i);
        y3(i,j) = y_ISS(i);
        theta3 = theta(i);
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
        E_future_ISS(i,j) = E_future_initial;
        tol = 1e-10;                   % Tolerance
        max_iter = 1000;              % Maximum iterations
        iter = 0; % # of iterations
        
        while true
            % Newton-Raphson
            f = E_future_ISS(i,j) - e_ISS * sin(E_future_ISS(i,j)) - M_ISS;  % Kepler's equation
            f_prime = 1 - e_ISS * cos(E_future_ISS(i,j));               % Derivative
            E_next(i,j) = E_future_ISS(i,j) - f / f_prime;
        
            if abs(E_next(i,j) - E_future_ISS(i,j)) < tol
                break;
            end
            if iter > max_iter
                error('Newton-Raphson failed to converge');
            end
        
            E_future_ISS(i,j) = E_next(i,j);
            iter = iter + 1;
            
        end
        
        E_future_ISS(i,j) = E_future_ISS(i,j)*180/pi;
        if E_future_ISS(i,j) >= 360
            E_future_ISS(i,j) = E_future_ISS(i,j) - 360;
            theta4(i,j) = (cosd(E_future_ISS(i,j)) - e_ISS) / (1 - e_ISS*cosd(E_future_ISS(i,j)));
            theta4(i,j) = round(acosd(theta4(i,j)));
        else
            theta4(i,j) = (cosd(E_future_ISS(i,j)) - e_ISS) / (1 - e_ISS*cosd(E_future_ISS(i,j)));
            theta4(i,j) = round(acosd(theta4(i,j)));
            theta4(i,j) = 360 - theta4(i,j);
        end    
        theta4(1,:) = theta4(end,1);

        index2(i,j) = find(theta == theta4(i,j));
        x4(i,j) = x_ISS(index2(i,j)); % future ISS coordinate
        y4(i,j) = y_ISS(index2(i,j)); % future ISS coordinate

        distance(i,j) = sqrt(((x4(i,j)-x2(i,j))^2) + (y4(i,j)-y2(i,j))^2);
        if (distance(i,j) < 0.1) 
            meet_senarios = [meet_senarios; x1(i,j), y1(i,j), x2(i,j), y2(i,j), x3(i,j), y3(i,j), x4(i,j), y4(i,j)];
        end
    end
end    

%% Plotting
figure;
% Orbits
plot(x_low, y_low, 'g', LineWidth=1.5)
hold on;
plot(x_ISS, y_ISS, 'b', LineWidth=1.5)
tran_entry = plot(meet_senarios(:,1),meet_senarios(:,2), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'c');
randevous = plot(meet_senarios(:,7),meet_senarios(:,8), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'r');

% Earth
earthBlue = [0, 0.2, 0.6]; % Earth's blue color
fill(x_Earth, y_Earth, earthBlue);
axis equal;
ylim([-9000 9000])
xlim([-9000 9000])
xlabel('x-axis [km]')
ylabel('y-axis [km]')
set(gca, 'Color', 'k'); % Set plot background to black
set(gca, 'FontSize', 24); % Increase font size

legendObj = legend([tran_entry, randevous], ...
    'Transfer orbit entries for randezvous', 'Possible randezvous locations', ...
    'Location', 'northeast', 'TextColor', [1, 1, 0], 'FontSize', 24);
set(legendObj, 'Color', 'k');

% tran_entry = plot(x1,y1, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'c'); 
% tran_exit = plot(x2,y2, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% initial_ISS = plot(x3,y3, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'm'); 
% fut_ISS = plot(x4,y4, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
% legend([tran_entry, tran_exit, initial_ISS, fut_ISS], 'Transfer orbit entry',...
%     'Transfer orbit exit', 'ISS initial coordinate', 'ISS future coordinate', ...
%     'Location', 'northeast', 'TextColor', [1, 1, 0], 'FontSize', 24);

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

