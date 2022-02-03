function particles = init_beacons_position(particles)
% Initialization of beacons position
N_PARTICLES = size(particles,2);

% Initialize Lm[1] (beacon in the corridor)
    for i_particle=1:N_PARTICLES
        particles{i_particle}.Lm(1,:) = [10.03, 4];
        particles{i_particle}.LmP{1} = [[0.1, 0.1]; [0.1, 0.1]];
    end

    % Initialize Lm[2] (beacon in the living room)
    for i_particle=1:N_PARTICLES
        particles{i_particle}.Lm(2,:) = [13, 7.5585];
        particles{i_particle}.LmP{2} = [[0.1, 0.1]; [0.1, 0.1]];
    end

    % Initialize Lm[3] (beacon in the kitchen)
    for i_particle=1:N_PARTICLES
        particles{i_particle}.Lm(3,:) = [20.12, 11.1];
        particles{i_particle}.LmP{3} = [[0.1, 0.1]; [0.1, 0.1]];
    end

    % Initialize Lm[4] (beacon in the dining table)
    for i_particle=1:N_PARTICLES
        particles{i_particle}.Lm(4,:) = [17.3, 9.8];
        particles{i_particle}.LmP{4} = [[0.1, 0.1]; [0.1, 0.1]];
    end
    
    % Initialize Lm[5] (beacon in the chair)
    for i_particle=1:N_PARTICLES
        particles{i_particle}.Lm(5,:) = [15, 8.3];
        particles{i_particle}.LmP{5} = [[0.1, 0.1]; [0.1, 0.1]];
    end
    
    % Initialize Lm[6] (beacon in the door)
    for i_particle=1:N_PARTICLES
        particles{i_particle}.Lm(6,:) = [3.3, 5];
        particles{i_particle}.LmP{6} = [[0.1, 0.1]; [0.1, 0.1]];
    end
    
    % Initialize Lm[7] (beacon in the toilet)
    for i_particle=1:N_PARTICLES
        particles{i_particle}.Lm(7,:) = [15, 3.5];
        particles{i_particle}.LmP{7} = [[0.1, 0.1]; [0.1, 0.1]];
    end

    % Initialize Lm[8] (beacon in the water tap)
    for i_particle=1:N_PARTICLES
        particles{i_particle}.Lm(8,:) = [13.8, 3.5];
        particles{i_particle}.LmP{8} = [[0.1, 0.1]; [0.1, 0.1]];
    end
    
    % Initialize Lm[9] (beacon in the pitcher)
    for i_particle=1:N_PARTICLES
        particles{i_particle}.Lm(9,:) = [21, 9.5];
        particles{i_particle}.LmP{9}= [[0.1, 0.1]; [0.1, 0.1]];
    end

    % Initialize Lm[10] (beacon in the broom)
    for i_particle=1:N_PARTICLES
        particles{i_particle}.Lm(10,:) = [16.3, 6.3];
        particles{i_particle}.LmP{10} = [[0.1, 0.1]; [0.1, 0.1]];
    end

end