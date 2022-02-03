clc;
clear all;
close all;

%participant = 1; % 
distance_errors_total = cell(1,24);

for participant=1:1
    clearvars -except  participant distance_errors_total;
    if participant== 10 || participant==15
        continue
    end

    load('participants_acc_gyr_bias_angle.mat')
    participant_id = participants_acc_gyr_bias_angle(participant, 1);
    acc_gyr_gravity_angle = participants_acc_gyr_bias_angle(participant, 2:end);

    [positions, step_events, beac_rssi_fixed_filtered, beac_rssi_activity_filtered, beac_motion] = Load_data(participant, acc_gyr_gravity_angle);
    step_events = [1, step_events, length(beac_motion)]; % Add '1' and the lenght of the dataset as points of support

    % Get strides where the beacons were moved. That indicates that it is the begining of a new trajectory
    [trajectory_parts] = GetTrajectoryParts(beac_motion,step_events);

    figure(participant)
    new_map

    N_ITERATIONS = 1;
    N_PARTICLES = 50;
    N_LM = 7;                               % Number of landmarks (beacons used)
    LM_SIZE = 2;                            % each landmark has x and y coordinates
    THRESHOLD_RESAMPLE = N_PARTICLES / 2;
    STD_PERSON_POSITION = 0.1;
    distance_errors = zeros(N_ITERATIONS, size(trajectory_parts,1));
    
    % code snippet to get the average of the 10 trajectories
    % All trajectories are composed by 6 sub-trajectories (x,y) axis
    
    sub_trajectories = cell(6,2);  
    
    % x-axis
    sub_trajectories{1,1} = zeros(N_ITERATIONS, trajectory_parts(1,2) - trajectory_parts(1,1) + 2);
    sub_trajectories{2,1} = zeros(N_ITERATIONS, trajectory_parts(2,2) - trajectory_parts(2,1) + 2);
    sub_trajectories{3,1} = zeros(N_ITERATIONS, trajectory_parts(3,2) - trajectory_parts(3,1) + 2);
    sub_trajectories{4,1} = zeros(N_ITERATIONS, trajectory_parts(4,2) - trajectory_parts(4,1) + 2);
    sub_trajectories{5,1} = zeros(N_ITERATIONS, trajectory_parts(5,2) - trajectory_parts(5,1) + 2);
    sub_trajectories{6,1} = zeros(N_ITERATIONS, trajectory_parts(6,2) - trajectory_parts(6,1) + 2);
    
    % y-axis
    sub_trajectories{1,2} = zeros(N_ITERATIONS, trajectory_parts(1,2) - trajectory_parts(1,1) + 2);
    sub_trajectories{2,2} = zeros(N_ITERATIONS, trajectory_parts(2,2) - trajectory_parts(2,1) + 2);
    sub_trajectories{3,2} = zeros(N_ITERATIONS, trajectory_parts(3,2) - trajectory_parts(3,1) + 2);
    sub_trajectories{4,2} = zeros(N_ITERATIONS, trajectory_parts(4,2) - trajectory_parts(4,1) + 2);
    sub_trajectories{5,2} = zeros(N_ITERATIONS, trajectory_parts(5,2) - trajectory_parts(5,1) + 2);
    sub_trajectories{6,2} = zeros(N_ITERATIONS, trajectory_parts(6,2) - trajectory_parts(6,1) + 2);

    %% MAIN LOOP
    for i_iteration=1:N_ITERATIONS
        disp(i_iteration)

        errors=zeros(1,size(trajectory_parts,1));
        for i_trajectory_part=1:size(trajectory_parts,1)

            initial_step = trajectory_parts(i_trajectory_part,1);
            final_step = trajectory_parts(i_trajectory_part,2);
            pos_trajectory_part = positions(initial_step: final_step,:);

            x_initial = trajectory_parts(i_trajectory_part,3);
            y_initial = trajectory_parts(i_trajectory_part,4);
            x_final =  trajectory_parts(i_trajectory_part,5);
            y_final =  trajectory_parts(i_trajectory_part,6);

            % To correct error position for the current trajectory.       
            correction_x = pos_trajectory_part(1, 1) - x_initial;
            correction_y = pos_trajectory_part(1, 2) - y_initial;        
            pos_trajectory_part(:, 1:2) = pos_trajectory_part(: ,1:2) - [correction_x correction_y];


            % CREATION OF PARTICLES
            particles = cell(1,N_PARTICLES);
            trajectories = cell(1,N_PARTICLES);
            landmarks_figures = cell(1,N_LM);
            particles_weights = ones(1,N_PARTICLES) / N_PARTICLES;
            for i_particle=1:N_PARTICLES
                particles{i_particle} = Particle(N_PARTICLES,N_LM,LM_SIZE, x_initial, y_initial);
                trajectory = plot(nan, nan, 'color', [.9 .9 .1]); 
                trajectories{i_particle} = trajectory;
            end
            for i_landmark=1:N_LM
                landmark_figure = plot(nan, nan,'d');
                landmarks_figures{i_landmark} = landmark_figure;
            end

            particles = init_beacons_position(particles);
            current_best_particle = N_PARTICLES;    

            x_prev = x_initial;
            y_prev = y_initial;

            for i_step=2:size(pos_trajectory_part,1)

                % PREDICTION
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         step_lenght = sqrt((pos_trajectory_part(i_step, 1) - pos_trajectory_part(i_step-1, 1))^2 +...
        %                       (pos_trajectory_part(i_step, 2) - pos_trajectory_part(i_step-1, 2))^2);
        %                   
        %         step_angle = atan((pos_trajectory_part(i_step, 2) - pos_trajectory_part(i_step-1, 2)) /...
        %                      (pos_trajectory_part(i_step, 1) - pos_trajectory_part(i_step-1, 1)));   

                x = pos_trajectory_part(i_step, 1);
                y = pos_trajectory_part(i_step, 2);
                for i_particle=1:N_PARTICLES

    %                 x = particles{i_particle}.X + (pos_trajectory_part(i_step, 1) - pos_trajectory_part(i_step-1, 1));
    %                 y = particles{i_particle}.Y + (pos_trajectory_part(i_step, 2) - pos_trajectory_part(i_step-1, 2));

                    dist_x = (x - x_prev) + randn() * STD_PERSON_POSITION;
                    dist_y = (y - y_prev) + randn() * STD_PERSON_POSITION;
                    particles{i_particle}.X = particles{i_particle}.X + dist_x;
                    particles{i_particle}.Y = particles{i_particle}.Y + dist_y;

                    particles{i_particle}.T_x = [particles{i_particle}.T_x, particles{i_particle}.X];
                    particles{i_particle}.T_y = [particles{i_particle}.T_y, particles{i_particle}.Y];

                    trajectory_x_data = [particles{i_particle}.T_x];
                    trajectory_y_data = [particles{i_particle}.T_y];

                    set(trajectories{i_particle},'XData',trajectory_x_data, 'YData',trajectory_y_data, 'color', [0.31, 0.31, 0.31]);  %[0.7, 0.7, 0.7]

                end
                

                x_prev = x;
                y_prev = y;
                pause(1)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % UPDATE
                initial_sample = step_events(i_step + initial_step -1);
                final_sample = step_events(i_step + initial_step);
                % Do not include RSSI measurements while person is still % EVALUAR
                if (final_sample - initial_sample) > 1024
                    initial_sample = final_sample - 1024;
                end

                rssi_corridor = beac_rssi_fixed_filtered(initial_sample:final_sample,1);
                rssi_living = beac_rssi_fixed_filtered(initial_sample:final_sample,2);
                rssi_kitchen = beac_rssi_fixed_filtered(initial_sample:final_sample,3);
                rssi_table = beac_rssi_fixed_filtered(initial_sample:final_sample,4);
                rssi_chair = beac_rssi_fixed_filtered(initial_sample:final_sample,5);
                rssi_door = beac_rssi_fixed_filtered(initial_sample:final_sample,6);

                rssi_toilet = beac_rssi_activity_filtered(initial_sample:final_sample,1);
                rssi_watertap = beac_rssi_activity_filtered(initial_sample:final_sample,2);
                rssi_pitcher = beac_rssi_activity_filtered(initial_sample:final_sample,3);
                rssi_broom = beac_rssi_activity_filtered(initial_sample:final_sample,4);

                rssis = [rssi_corridor, rssi_living, rssi_kitchen, rssi_table, rssi_chair, rssi_door, rssi_toilet, rssi_watertap];

                for i_rssis=1:size(rssis,1)

                    non_zero_indexes = find(rssis(i_rssis,:));
                    if ~isempty(non_zero_indexes)

                        for j_rssi=non_zero_indexes
                            rssi = rssis(i_rssis,j_rssi);

                            if j_rssi <9 %j_rssi ~= 2 %Stationary beacons
                                if rssi >= -88
                                    [particles, particles_weights, trajectories, current_best_particle] = update_resample(particles, rssi, j_rssi, particles_weights, trajectories, THRESHOLD_RESAMPLE);
                                end
                            end                   
                        end
                    end
                end

            end
            
            % Correction of orientation
            best_trajectory_positions = zeros(size(trajectories{current_best_particle}.XData,2),2);
            best_trajectory_positions(:,1) = trajectories{current_best_particle}.XData;
            best_trajectory_positions(:,2) = trajectories{current_best_particle}.YData;

            Pos_delta = zeros(size(best_trajectory_positions,1),2);
            StrideLengths = zeros(size(best_trajectory_positions,1),1);
            for i_step=1:size(best_trajectory_positions,1)-1
                Pos_delta(i_step,:) = best_trajectory_positions(i_step+1,1:2)-best_trajectory_positions(i_step,1:2);
                StrideLengths(i_step)=sqrt(Pos_delta(i_step,1)^2+Pos_delta(i_step,2)^2); % stride lenght in X-Y plane
            end

            Thetas = zeros(size(best_trajectory_positions,1),1);
            for k=1:size(best_trajectory_positions,1) % Thetas (detected angle from position increments)
                Thetas(k)=atan2(Pos_delta(k,2),Pos_delta(k,1)); % orientations orientacion with respect to the north
            end

            [best_trajectory_positions_corrected, distance_error] = getDegrees(StrideLengths, Thetas, [x_initial, y_initial], [x_final, y_final]);
            errors(i_trajectory_part) = distance_error;
            plot(best_trajectory_positions_corrected(:,1),best_trajectory_positions_corrected(:,2), 'color', [0.7, 0.7, 0.7]);
            
            % to erase particle trajectories
            for i_particle=1:N_PARTICLES
                if i_particle ~= current_best_particle
                set(trajectories{i_particle},'XData',[], 'YData',[]);
                end
            end        
            set(trajectories{current_best_particle},'XData',[], 'YData',[]);

            % Save the positions (x,y) of the sub-trajectory
            sub_trajectories{i_trajectory_part,1}(i_iteration,:) = best_trajectory_positions_corrected(:,1).'; % x-axis
            sub_trajectories{i_trajectory_part,2}(i_iteration,:) = best_trajectory_positions_corrected(:,2).'; % y-axis
            
%             % Plot the landmarks of the best particle.
%             for i_lm=1:N_LM
%                 set(landmarks_figures{i_lm},'XData',particles{current_best_particle}.Lm(i_lm,1), 'YData',particles{current_best_particle}.Lm(i_lm,2));
%             end

        end
        distance_errors(i_iteration,1:length(errors)) = errors;
    end
    distance_errors_total{1,participant} = distance_errors;
    
    % Plot the average trajectory of the ten iterations
    for i_sub_trajectory=1:6
        plot(mean(sub_trajectories{i_sub_trajectory,1}(:,:)),mean(sub_trajectories{i_sub_trajectory,2}(:,:)),'color','b','LineWidth',1);
    end
    % to erase particle trajectories
    for i_particle=1:N_PARTICLES
        set(trajectories{i_particle},'XData',[], 'YData',[]);
    end   
end



function [particles, particles_weights, trajectories, current_best_particle] = ...
    update_resample(particles, rssi, j_rssi, particles_weights, trajectories, THRESHOLD_RESAMPLE)
    
    % To update
    [particles, particles_weights, current_best_particle] = update(particles, rssi, j_rssi, particles_weights);

    % To resample
    neff = 1 / sum(sqrt(particles_weights));
    if neff < THRESHOLD_RESAMPLE
        %indexes = double(py.mifunc.stratified_resample(particles_weights))+1;
        indexes = double(py.mifunc.systematic_resample(particles_weights))+1;
        [particles, particles_weights, trajectories] = resample_from_index(particles, particles_weights, indexes, trajectories);
        %particles_weights = ones(1,N_PARTICLES) / N_PARTICLES; % PARTICLES WEIGHTS ARE INITILIZED TO AND EQUAL VALUE

    end
end

function [particles, particles_weights, current_best_particle] = update(particles, rssi, lm_id, particles_weights)
    
    % different models of beacons have different propagation parameters
    if lm_id <= 4 % Big beacons
        distance_rssi = 10 ^ ((rssi + 80) / -10 * 0.6);  % Calculate distance to beacon
        R = distance_rssi / 10;
    else
        distance_rssi = 10 ^ ((rssi + 85) / -10 * 0.6);  % Calculate distance to beacon
        R = distance_rssi / 10;
    end
    
    N_PARTICLES = size(particles,2);
    for i_particle=1:N_PARTICLES
        % Get the distance between the person and the landmark
        dx = particles{i_particle}.X - particles{i_particle}.Lm(lm_id,1);
        dy = particles{i_particle}.Y - particles{i_particle}.Lm(lm_id,2);
        dist = sqrt((dx * dx) + (dy * dy));
        residual = distance_rssi - dist;

        % Compute Jacobians
        H = [-dx / dist, -dy / dist];

        % Compute covariance of the residual
        % covV = H * Cov_s * H^T + error
        HxCov = [particles{i_particle}.LmP{1,lm_id}(1, 1) * H(1) + particles{i_particle}.LmP{1,lm_id}(1, 2) * H(2),...
                 particles{i_particle}.LmP{1,lm_id}(2, 1) * H(1) + particles{i_particle}.LmP{1,lm_id}(2, 2) * H(2)];

        covV = (HxCov(1) * H(1)) + (HxCov(2) * H(2)) + R;

        % Calculate Kalman gain
        K_gain = [HxCov(1) * (1 / covV), HxCov(2) * (1.0 / covV)];

        % Calculate the new landmark position
        lm_x = particles{i_particle}.Lm(lm_id,1) + (K_gain(1) * residual);
        lm_y = particles{i_particle}.Lm(lm_id,2) + (K_gain(2) * residual);

        % Calculate the new covariance matrix of the landmark
        % cov_t = cov_t-1 - K * covV * K^T
        lm_P_aux = [[K_gain(1) * K_gain(1) * covV, K_gain(1) * K_gain(2) * covV];...
                    [K_gain(2) * K_gain(1) * covV, K_gain(2) * K_gain(2) * covV]];

        lm_P = [[particles{i_particle}.LmP{1,lm_id}(1, 1) - lm_P_aux(1,1),...
                 particles{i_particle}.LmP{1,lm_id}(1, 2) - lm_P_aux(1,2)];...
                [particles{i_particle}.LmP{1,lm_id}(2, 1) - lm_P_aux(2,1),...
                 particles{i_particle}.LmP{1,lm_id}(2, 2) - lm_P_aux(2,2)]];

        % Update landmark in particle
        particles{i_particle}.Lm(lm_id, 1) = lm_x;
        particles{i_particle}.Lm(lm_id, 2) = lm_y;
        particles{i_particle}.LmP{1,lm_id} = lm_P;

        % Update particles weight
        particles_weights(i_particle) = py.mifunc.calculate_weight(particles_weights(i_particle), dist, covV, distance_rssi);
        %particles_weights(i_particle) = (particles_weights(i_particle) * (1 / (covV * sqrt(2*pi))) * exp((-(distance_rssi - dist)^2) / (2 * covV^2))) + 1*exp(-300);
    end
    
    %Save particle weights in all the particles
    particles_weights = particles_weights / sum(particles_weights);
    biggest_weight = -1;
    for i_particle=1:N_PARTICLES
        particles{i_particle}.W = particles_weights(i_particle);
        
        % Get the best particle (particle with the biggest weight)
        if particles_weights(i_particle) >= biggest_weight
            biggest_weight = particles_weights(i_particle);
            current_best_particle = i_particle;
        end

    end
end

function [new_particles, new_particles_weights, trajectories] = resample_from_index(particles, particles_weights, indexes, trajectories)
    
    N_PARTICLES = size(particles,2);    
    % Resample particles
    new_particles_weights = zeros(1,N_PARTICLES);
    new_particles = particles;
    for i_particle = 1:N_PARTICLES
        new_particles_weights(i_particle) = particles_weights(indexes(i_particle));
        new_particles{i_particle} = particles{indexes(i_particle)};
    end
    
    % Update particle weights and plot trajectories
    new_particles_weights = new_particles_weights / sum(new_particles_weights);
    for i_particle=1:N_PARTICLES    

        new_particles{i_particle}.W = new_particles_weights(i_particle);        
        
        if ~ismember(i_particle,indexes)
            set(trajectories{i_particle},'XData',[], 'YData',[]);
        else
            trajectory_x_data = [new_particles{i_particle}.T_x];
            trajectory_y_data = [new_particles{i_particle}.T_y];
            set(trajectories{i_particle},'XData',trajectory_x_data, 'YData',trajectory_y_data);
        end

    end
end