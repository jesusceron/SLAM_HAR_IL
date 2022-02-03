function [positions]=ZUPT_KF(acc, gyr_unbiased, g, acc_mean, fs, Step_events, angle, participant)
%% Read data from file.
% Data should include timestamps (seconds), 3 axis accelerations (m/s^2), 3
% axis gyroscopic rates of turn (rad/s).
% load('.\log_files\timestamps_KF_evaluation.mat');
% load('.\log_files\acc_KF_evaluation.mat');
% load('.\log_files\gyr_KF_evaluation.mat');
% load('.\log_files\stancePhase_KF_evaluation.mat');

data_size = size(acc,1);
Acc = acc';
Gyro = gyr_unbiased';


%% Initialise parameters.
% Orientation from accelerometers. Sensor is assumed to be stationary.
pitch =-atan2(acc_mean(1),sqrt(acc_mean(2)^2+acc_mean(3)^2));
roll = atan2(acc_mean(2),acc_mean(3));  
yaw = 0;

C = [cos(pitch)*cos(yaw) (sin(roll)*sin(pitch)*cos(yaw))-(cos(roll)*sin(yaw)) (cos(roll)*sin(pitch)*cos(yaw))+(sin(roll)*sin(yaw));
    cos(pitch)*sin(yaw)  (sin(roll)*sin(pitch)*sin(yaw))+(cos(roll)*cos(yaw))  (cos(roll)*sin(pitch)*sin(yaw))-(sin(roll)*cos(yaw));
    -sin(pitch) sin(roll)*cos(pitch) cos(roll)*cos(pitch)];
C_prev = C;

% Preallocate storage for heading estimate. Different from direction of
% travel, the heading indicates the direction that the sensor, and therefore
% the pedestrian, is facing.
heading = nan(1, data_size);
heading(1) = yaw;

% Gyroscope bias, to be determined for each sensor.
%  -- Defined above so we don't forget to change for each dataset. --

% Preallocate storage for accelerations in navigation frame.
acc_n = nan(3, data_size);
acc_n(:,1) = C*Acc(:,1);


% Preallocate storage for velocity (in navigation frame).
% Initial velocity assumed to be zero.
vel_n = nan(3, data_size);
vel_n(:,1) = [0 0 0]';

% Preallocate storage for position (in navigation frame).
% Initial position arbitrarily set to the origin.
pos_n = nan(3, data_size);
pos_n(:,1) = [0   0   0]';

% Preallocate storage for distance travelled used for altitude plots.
distance = nan(1,data_size-1);
distance(1) = 0;


% Error covariance matrix.
P = zeros(9);

% Process noise parameter, gyroscope and accelerometer noise.
sigma_omega = 0.02;
sigma_a = 0.03;

% ZUPT measurement matrix.
H = [zeros(3) zeros(3) eye(3)];

% ZUPT measurement noise covariance matrix.
sigma_v = 0.005;
R = diag([sigma_v sigma_v sigma_v]).^2;

% Gyroscope stance phase detection threshold.
gyro_threshold = 0.6;

%% Main Loop
for t = 2:data_size
    %%% Start INS (transformation, double integration) %%%
    dt = 1/fs;
    
    % Remove bias from gyro measurements.
    gyro_s1 = Gyro(:,t); %- gyro_bias;
    
    % Skew-symmetric matrix for angular rates
    ang_rate_matrix = [0   -gyro_s1(3)   gyro_s1(2);
        gyro_s1(3)  0   -gyro_s1(1);
        -gyro_s1(2)  gyro_s1(1)  0];
    
    % orientation estimation
    C = C_prev*(2*eye(3)+(ang_rate_matrix*dt))/(2*eye(3)-(ang_rate_matrix*dt));
    
    % Transforming the acceleration from sensor frame to navigation frame.
    acc_n(:,t) = 0.5*(C + C_prev)*Acc(:,t);
    
    % Velocity and position estimation using trapeze integration.
    vel_n(:,t) = vel_n(:,t-1) + ((acc_n(:,t) - [0; 0; g] )+(acc_n(:,t-1) - [0; 0; g]))*dt/2;
    pos_n(:,t) = pos_n(:,t-1) + (vel_n(:,t) + vel_n(:,t-1))*dt/2;
    
    % Skew-symmetric cross-product operator matrix formed from the n-frame accelerations.
    S = [0  -acc_n(3,t)  acc_n(2,t);
        acc_n(3,t)  0  -acc_n(1,t);
        -acc_n(2,t) acc_n(1,t) 0];
    
    % State transition matrix.
    F = [eye(3)  zeros(3,3)    zeros(3,3);
        zeros(3,3)   eye(3)  dt*eye(3);
        -dt*S  zeros(3,3)    eye(3) ];
    
    % Compute the process noise covariance Q.
    Q = diag([sigma_omega sigma_omega sigma_omega 0 0 0 sigma_a sigma_a sigma_a]*dt).^2;
    
    % Propagate the error covariance matrix.
    P = F*P*F' + Q;
    %%% End INS %%%
    
    % Stance phase detection and zero-velocity updates.
    if norm(Gyro(:,t)) < gyro_threshold
        %%% Start Kalman filter zero-velocity update %%%
        % Kalman gain.
        K = (P*(H)')/((H)*P*(H)' + R);
        
        % Update the filter state.
        delta_x = K*vel_n(:,t);
        
        % Update the error covariance matrix.
        P = (eye(9) - K*H)*P;
        
        % Extract errors from the KF state.
        attitude_error = delta_x(1:3);
        pos_error = delta_x(4:6);
        vel_error = delta_x(7:9);
        %%% End Kalman filter zero-velocity update %%%
        
        %%% Apply corrections to INS estimates. %%%
        % Skew-symmetric matrix for small angles to correct orientation.
        ang_matrix = -[0   -attitude_error(3,1)   attitude_error(2,1);
            attitude_error(3,1)  0   -attitude_error(1,1);
            -attitude_error(2,1)  attitude_error(1,1)  0];
        
        % Correct orientation.
        C = (2*eye(3)+(ang_matrix))/(2*eye(3)-(ang_matrix))*C;
        
        % Correct position and velocity based on Kalman error estimates.
        vel_n(:,t)=vel_n(:,t)-vel_error;
        pos_n(:,t)=pos_n(:,t)-pos_error;
    end
    heading(t) = atan2(C(2,1), C(1,1)); % Estimate and save the yaw of the sensor (different from the direction of travel). Unused here but potentially useful for orienting a GUI correctly.
    C_prev = C; % Save orientation estimate, required at start of main loop.
    
    % Compute horizontal distance.
    distance(1,t) = distance(1,t-1) + sqrt((pos_n(1,t)-pos_n(1,t-1))^2 + (pos_n(2,t)-pos_n(2,t-1))^2);
end



%Rotate position estimates and plot.
rotation_matrix = [cosd(angle) -sind(angle);  % Rotation angle required to achieve an aesthetic alignment of the figure.
    sind(angle) cosd(angle)];
pos_r = zeros(2,length(pos_n));
for idx = 1:length(pos_n)
    pos_r(:,idx) = rotation_matrix*[pos_n(1,idx) pos_n(2,idx)]';
end
pos_r(1,:) = pos_r(1,:) + 3.2;
pos_r(2,:) = pos_r(2,:) + 7;
pos_r(3,:) = pos_n(3,:)';
positions_kf = [pos_r(1,:); pos_r(2,:); pos_r(3,:)]'; % minus to pass from west to east

% Plot in step basis
number_of_steps = length(Step_events);
positions = zeros(number_of_steps+1,3);
positions(1,:) = [3.2 7 0]; % initial position
distances = zeros(number_of_steps,1);
for i_step=1:number_of_steps
    sample_step = Step_events(i_step);
    positions(i_step+1,:) = [positions_kf(sample_step,1), positions_kf(sample_step,2), positions_kf(sample_step,3)]; % minus to pass from west to east
    distances(i_step,1) = distance(sample_step);
end


% % plot figure
%     %openfig('map.fig','new'); % open figure
%     %map2
%     %hold on
%     %plot(Thetas(:,1),Positions(:,1),'g');
    
%     % Plot in 2d in step basis
%     figure(participant+30)
%     new_map
%     plot(positions(:,1),positions(:,2),'bo-');
%     start = plot(positions(1,1),positions(1,2),'bs','MarkerSize',8,'MarkerFaceColor',[0 0 1]); % marcar posicion inicial
%     stop = plot(positions(end,1),positions(end,2),'bo','MarkerSize',8,'MarkerFaceColor',[0 0 1]); % marcar paso final
%     hold off
    
% %     %Plot in a continous way (complete KF)
% %     plot(positions_kf(:,1),positions_kf(:,2),'r');
% %     %plot(Positions(:,1),Positions(:,2),'LineWidth',1,'Color','b');
% %     start = plot3(positions_kf(1,1),positions_kf(1,2),positions_kf(1,3),'bs');
% %     stop = plot3(positions_kf(end,1),positions_kf(end,2),positions_kf(end,3),'bo');
% %     box on;
% %     hold off;


%% Plot altitude estimates.
% figure;
% box on;
% hold on;
% plot(distances,positions(:,3),'Linewidth',2, 'Color','b');
% xlabel('Distance Travelled (m)');
% ylabel('z (m)');
% title('Estimated altitude');
% grid;
% 
% % Display lines representing true altitudes of each floor.
% floor_colour = [0 0.5 0]; % Colour for lines representing floors.
% floor_heights = [0 3.6 7.2 10.8]; % Altitude of each floor measured from the ground floor.
% floor_names = {'A' 'B' 'C' 'D'};
% lim = xlim;
% for floor_idx = 1:length(floor_heights)
%     line(lim, [floor_heights(floor_idx) floor_heights(floor_idx)], 'LineWidth', 2, 'LineStyle', '--', 'Color', floor_colour);
% end
% ax1=gca; % Save handle to main axes.
% axes('YAxisLocation','right','Color','none','YTickLabel', floor_names, 'YTick', floor_heights,'XTickLabel', {});
% ylim(ylim(ax1));
% ylabel('Floor');
% hold off;
end