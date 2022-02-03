function [positions, Step_events, beac_rssi_fixed_filtered, beac_rssi_activity_filtered, beac_motion] = Load_data(participant, acc_gyr_gravity_angle)


    %% Load data   
    %Acc_Gyr_Beac_data= readtable(strcat('E:\Google Drive\matlab_slam3\dataset_synchronized\',int2str(participant_id),'.csv'));
    load(strcat('P',int2str(participant),'.mat'));
    
%     T_p=Acc_Gyr_Beac_data(:, {'motion_state_beacon_5', 'motion_state_beacon_6', 'motion_state_beacon_7', 'motion_state_beacon_8', 'motion_state_beacon_9', 'motion_state_beacon_10'});
%     a = sum(T_p{:,1:end});
    
%     plot(Acc_Gyr_Beac_data.Shimmer_93B6_Gyro_X_CAL,'DisplayName','Acc_Gyr_Beac_data.Shimmer_93B6_Gyro_X_CAL');hold on;plot(Acc_Gyr_Beac_data.Shimmer_93B6_Gyro_Y_CAL,'DisplayName','Acc_Gyr_Beac_data.Shimmer_93B6_Gyro_Y_CAL');plot(Acc_Gyr_Beac_data.Shimmer_93B6_Gyro_Z_CAL,'DisplayName','Acc_Gyr_Beac_data.Shimmer_93B6_Gyro_Z_CAL');plot(Acc_Gyr_Beac_data.motion_state_beacon_5*500,'DisplayName','Acc_Gyr_Beac_data.motion_state_beacon_5');hold off;
%     Acc_Gyr_Beac_data.motion_state_beacon_5(66069:end)=0;
%     
%     plot(Acc_Gyr_Beac_data.Shimmer_93B6_Gyro_X_CAL,'DisplayName','Acc_Gyr_Beac_data.Shimmer_93B6_Gyro_X_CAL');hold on;plot(Acc_Gyr_Beac_data.Shimmer_93B6_Gyro_Y_CAL,'DisplayName','Acc_Gyr_Beac_data.Shimmer_93B6_Gyro_Y_CAL');plot(Acc_Gyr_Beac_data.Shimmer_93B6_Gyro_Z_CAL,'DisplayName','Acc_Gyr_Beac_data.Shimmer_93B6_Gyro_Z_CAL');plot(Acc_Gyr_Beac_data.motion_state_beacon_5*500,'DisplayName','Acc_Gyr_Beac_data.motion_state_beacon_5');hold off;

    acc_mean = acc_gyr_gravity_angle(1:3);
    bias_gyr = acc_gyr_gravity_angle(4:6);
    gravity = acc_gyr_gravity_angle(7);
    angle = acc_gyr_gravity_angle(8);
    fs = 204.8; % IMU sample rate

    %% Load data.

    % Load IMU data
    acc_complete = table2array(Acc_Gyr_Beac_data(:,2:4));
    gyro_complete = table2array(Acc_Gyr_Beac_data(:,5:7));

    % Load BEACONS data.
    beac_rssi_fixed = table2array(Acc_Gyr_Beac_data(:,...
        {'RSSIs_beacon_1',...   % Coconut:      Pasillo
        'RSSIs_beacon_2',...    % Mint:         Sala
        'RSSIs_beacon_3',...    % Ice:          Cocina
        'RSSIs_beacon_4',...    % Blueberry:    Comedor
        'RSSIs_beacon_6',...    % P1:           Chair
        'RSSIs_beacon_9'}));    % G2:           Door

    beac_rssi_activity = table2array(Acc_Gyr_Beac_data(:,...
        {'RSSIs_beacon_7',...    % B2:          Toilet lid
        'RSSIs_beacon_8',...     % B1:          Water tap
        'RSSIs_beacon_5',...     % P2:          Pitcher
        'RSSIs_beacon_10'}));    % G1:          Broom

    beac_motion = table2array(Acc_Gyr_Beac_data(:,{'motion_state_beacon_5',...
        'motion_state_beacon_7',...
        'motion_state_beacon_8',...
        'motion_state_beacon_10'}));

    %% Data pre-processing
    acc = zeros(length(acc_complete),3);
    acc(:,1) = acc_complete(:,2); 
    acc(:,2) = acc_complete(:,3); % y = z
    acc(:,3) = acc_complete(:,1); 

    gyr = zeros(length(gyro_complete),3);
    gyr(:,1) = gyro_complete(:,2); 
    gyr(:,2) = gyro_complete(:,3); 
    gyr(:,3) = gyro_complete(:,1); 
    gyr = deg2rad(gyr);

    gyr_unbiased=[gyr(:,1) - bias_gyr(1),...
        gyr(:,2) - bias_gyr(2),...
        gyr(:,3) - bias_gyr(3)];


    % Filtering the beacons data in Python
    [beac_rssi_fixed_filtered, beac_rssi_activity_filtered] = rssiKF(beac_rssi_fixed,beac_rssi_activity);


    %% -----Step detection---------------------------
    idx_fig=20;
    [~,Step_events,~,~] = StepDetection_Acel(acc, 0, idx_fig);

    %stance_phase = load(strcat('E:\Google Drive\IL_HAR_app_Matlab\dataset_synchronized\stance_phase_p',int2str(participant),'.mat'));
    %Step_events = stance_phase.stance_phase(:,1)';


    %% -------------------- Trajectory reconstruction ZUPT KF -----------------------%
    [positions]=ZUPT_KF(acc, gyr_unbiased, gravity, acc_mean, fs, Step_events, angle, participant);
end

    