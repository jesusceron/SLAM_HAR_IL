function [trajectory_parts] = GetTrajectoryParts(beac_motion,step_events)
% Get strides where the beacons are moved. That indicates that is the
% begining of a new trajectory

% gap = 508; % movement readings have a delay of about 2 seconds (204.8Hz = 508 samples)
% pitcher_moving = [beac_motion(gap:end,1); zeros(gap-1,1)];
% toilet_moving = [beac_motion(gap:end,2); zeros(gap-1,1)];
% watertap_moving = [beac_motion(gap:end,3); zeros(gap-1,1)];
% broom_moving = [beac_motion(gap:end,4); zeros(gap-1,1)];

pitcher_moving = beac_motion(:,1);
toilet_moving = beac_motion(:,2);
watertap_moving = beac_motion(:,3);
broom_moving = beac_motion(:,4);


for i_beacon=1:4
    
    switch i_beacon
%         case 5 % door
%             door_movement = find(door_moving);
%             sample_first_door_movement = door_movement(1); % getting in to the room. p2:-520, p9:3800, p10:3200
%             % sample_second_door_movement = 59560; %p9:60300 P8:59560
%             for i_door_movement=1:size(door_movement,1)
%                 if door_movement(i_door_movement) > sample_first_door_movement + 20000
%                     sample_second_door_movement = door_movement(i_door_movement);% P2:-180 P5:-1400 p10:-300 P11:-900getting out to the room
%                 end
%             end
%             
%             for i_step=1:length(step_events)
%                 if step_events(i_step) > sample_first_door_movement
%                     first_step_in_door = i_step;
%                     door_position = [3.3 5];
%                     break
%                 end
%             end
%             
%             for i_step=1:length(step_events)
%                 if step_events(i_step) > sample_second_door_movement
%                     second_step_in_door = i_step;
%                     door_position = [3.3 5];
%                     break
%                 end
%             end
            
        case 2 % Toilet
            toilet_movement = find(toilet_moving);
            sample_toilet_movement = toilet_movement(1); % p8:12775;%
            for i_step=1:length(step_events)
                if step_events(i_step) > sample_toilet_movement
                    step_in_toilet = i_step;
                    toilet_position = [15 3.5];
                    break
                end
            end
            
        case 4 % Broom
            broom_movement = find(broom_moving);
            sample_broom_movement = broom_movement(1); % p8:44225;%
            for i_step=1:length(step_events)
                if step_events(i_step) > sample_broom_movement
                    step_in_broom = i_step ;
                    broom_position = [16.3 6.3];
                    break
                end
            end
            
        case 1 % pitcher
            pitcher_movement = find(pitcher_moving);
            sample_first_pitcher_movement = pitcher_movement(1); % getting in to the room. p2:-520, p9:3800, p10:3200
            % sample_second_pitcher_movement = 59560; %p9:60300 P8:59560
            for i_pitcher_movement=1:size(pitcher_movement,1)
                if pitcher_movement(i_pitcher_movement) > sample_first_pitcher_movement + 20000
                    sample_second_pitcher_movement = pitcher_movement(i_pitcher_movement);% P2:-180 P5:-1400 p10:-300 P11:-900getting out to the room
                end
            end
            
            for i_step=1:length(step_events)
                if step_events(i_step) > sample_first_pitcher_movement
                    first_step_in_pitcher = i_step;
                    first_pitcher_position = [21 9.5];
                    break
                end
            end
            
            for i_step=1:length(step_events)
                if step_events(i_step) > sample_second_pitcher_movement
                    second_step_in_pitcher = i_step;
                    second_pitcher_position = [19.64 4.3];
                    break
                end
            end
            
%         case 2 % TV controller
%             tvcontroller_movement = find(TVcontroller_moving);
%             sample_tvcontroller_movement = tvcontroller_movement(1); % p8:42271;%
%             for i_step=1:length(step_events)
%                 if step_events(i_step) > sample_tvcontroller_movement
%                     step_in_tvcontroller = i_step;
%                     tvcontroller_position = [8.7 6.8];
%                     break
%                 end
%             end
            
         case 3 % Water tap
            watertap_movement = find(watertap_moving);
            sample_watertap_movement = watertap_movement(1); % p8:42271;%
            for i_step=1:length(step_events)
                if step_events(i_step) > sample_watertap_movement
                    step_in_watertap = i_step;
                    watertap_position = [13.8 3.5];
                    break
                end
            end
    end
end

a = [step_in_toilet; step_in_broom; first_step_in_pitcher; second_step_in_pitcher; step_in_watertap];
pos = [toilet_position; broom_position; first_pitcher_position; second_pitcher_position; watertap_position];
[in, b] = sort(a);
p = [in pos(b,:)];

init_step = [1; p(:,1)];
final_step = [p(:,1);length(step_events)-1];
init_pos = [[3.2, 7]; p(:,2:3)];
final_pos = [p(:,2:3);[3, 7.5]];

trajectory_parts = [init_step, final_step, init_pos, final_pos];