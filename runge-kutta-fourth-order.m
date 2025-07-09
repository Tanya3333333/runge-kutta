clear;   clc;  close all;  

% Constants
g = 9.81; 
l = 0.12;
m = 0.0242;
I_G = 2.9E-5;
h = 0.1;  % time step
t = 0:h:8;  % (0 to 8 seconds)
N = length(t);  

% Define initial angles for each set
theta0_list = [10.01, 18.2, 41.86, 58.422, 71.162];  % Initial angles in degrees
theta0_list_rad = deg2rad(theta0_list);  % Convert to radians

% Loop through each initial condition
for j = 1:length(theta0_list_rad)
    theta1 = zeros(1, N);  % angle (radians)
    theta2 = zeros(1, N);  % angular velocity (rad/s)
    theta1(1) = theta0_list_rad(j);  % initial angle (radians)
    theta2(1) = 0;  % initial angular velocity

    % Runge-Kutta 4th-order method
    for i = 1:N-1   
        
        % k for angle and angular velocity
        k1_1 = h * theta2(i); 
        k1_2 = h * (-g/l * sin(theta1(i)));
        
        k2_1 = h * (theta2(i) + 0.5 * k1_2); 
        k2_2 = h * (-g/l * sin(theta1(i) + 0.5 * k1_1));
        
        k3_1 = h * (theta2(i) + 0.5 * k2_2);
        k3_2 = h * (-g/l * sin(theta1(i) + 0.5 * k2_1));
        
        k4_1 = h * (theta2(i) + k3_2);
        k4_2 = h * (-g/l * sin(theta1(i) + k3_1)); 
                  
          

        % Update the values 
        theta1(i+1) = theta1(i) + (1/6) * (k1_1 + 2*k2_1 + 2*k3_1 + k4_1);  % angle
        theta2(i+1) = theta2(i) + (1/6) * (k1_2 + 2*k2_2 + 2*k3_2 + k4_2);  % angular velocity


        disp(['Initial Angle: ', num2str(theta0_list(j)), 'degree']);
        disp(['Step ', num2str(i), ':']);
        disp(['  k1_1 = ', num2str(k1_1), ', k1_2 = ', num2str(k1_2)]);
        disp(['  k2_1 = ', num2str(k2_1), ', k2_2 = ', num2str(k2_2)]);
        disp(['  k3_1 = ', num2str(k3_1), ', k3_2 = ', num2str(k3_2)]);
        disp(['  k4_1 = ', num2str(k4_1), ', k4_2 = ', num2str(k4_2)]);
    end

    % Convert back to degrees for plots
    angle_deg = rad2deg(theta1); 
    velocity_deg = rad2deg(theta2); 

    % Plots
    figure;
    subplot(2,1,1);
    plot(t, angle_deg);  % angle vs time
    xlabel('Time (s)');
    ylabel('Angle (°)');
    title(['Pendulum Angle vs Time (Initial Angle: ', num2str(theta0_list(j)), '°)']);

    subplot(2,1,2);
    plot(t, velocity_deg);  % angular velocity vs time
    xlabel('Time (s)');
    ylabel('Angular Velocity (deg/s)');
    title(['Pendulum Angular Velocity vs Time (Initial Angle: ', num2str(theta0_list(j)), '°)']);
end
