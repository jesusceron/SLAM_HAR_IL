%%% participants id's
%%% december 18: @25282489, @65745876, @76325018, @9335772
%%% december 19: @4610326, @1002819840, @25269362, @34555639, @76306709, @25399881, @25271446, @10290567, @34325893, @34327818, 76323691
%%% december 22: @10534499, @19056006, @41534025, @41331063, @31239627, @1004301404, @14873733, @10526517, @34539568

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all
participant_id = int2str(25399881);
file_shimmer_data = strcat(participant_id,'_shimmer_foot.csv');
file_smartphone_data = strcat(participant_id,'_acc_gyr_beac_sync.csv');
[shimmer_data, smartphone_data] = LoadData(file_shimmer_data, file_smartphone_data);
PlotRawData(shimmer_data, smartphone_data);

%%%%% (1) Align shimmer and smartphone data in the first jump %%%%%
%%% shimmer_data(1:peak_difference,:)=[];  OR
%%% smartphone_data(1:peak_difference,:)=[];
%%% comment line 13

%%%%% (2) Drop samples until the first jump %%%%%
%%% smartphone_data(1:first_peak_sample,:) = [];
%%% shimmer_data(1:first_peak_sample,:) = [];

%%%%% (3) cut samples manually in shimmer_data and smartphone data for them to begin exactly in the first sample of a second 
%%%%% (by using matlab Workspace) %%%%%

% [edges_sh, edges_sm] = GetSecondEdges(shimmer_data, smartphone_data);
% cd = SyncData(edges_sh,edges_sm,shimmer_data,smartphone_data);
% PlotDataProcessed(cd);

%%%%% (4) drop samples after third final jump manually %%%%%
%%% cd(last_third_jump:end,:) = [];

% ExportFile(cd, participant_id)
% PlotFinalData(participant_id)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [shimmer_data, smartphone_data] = LoadData(file_shimmer_data, file_smartphone_data)
    %% Read data
    shimmer_data= readtable(strcat('E:\PycharmProjects\DataPreProcess\dataset\', file_shimmer_data));
    smartphone_data= readtable(strcat('E:\PycharmProjects\DataPreProcess\dataset\', file_smartphone_data));

    %% Set datetimes for synchronization
    time_shimmer = datetime(shimmer_data.Shimmer_93B6_Timestamp_Shimmer_CAL,'ConvertFrom','epochtime','TimeZone','America/Bogota',...
                 'TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
    time_shimmer.Hour = time_shimmer.Hour-5; % Time in Colombia: UTC-5
    shimmer_data.Shimmer_93B6_Timestamp_Shimmer_CAL = time_shimmer;

    time_smartphone = datetime(smartphone_data.Timestamps_ref,'ConvertFrom','epochtime','TimeZone','America/Bogota',...
                 'TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
    time_smartphone.Hour = time_smartphone.Hour-5; % Time in Colombia: UTC-5
    smartphone_data.Timestamps_ref = time_smartphone;
end

%% Plot data
function PlotRawData(shimmer_data, smartphone_data)
    figure(1)
    hold on;
    %plot(shimmer_data.Shimmer_93B6_Accel_WR_X_CAL,'DisplayName','Shimmer_Accel_X');
    %plot(shimmer_data.Shimmer_93B6_Accel_WR_Y_CAL,'DisplayName','Shimmer_Accel_Y');
    %plot(shimmer_data.Shimmer_93B6_Accel_WR_Z_CAL,'DisplayName','Shimmer_Accel_Z');

    %plot(shimmer_data.Shimmer_93B6_Gyro_X_CAL,'DisplayName','Shimmer_Gyro_X');
    %plot(shimmer_data.Shimmer_93B6_Gyro_Y_CAL,'DisplayName','Shimmer_Gyro_Y');
    plot(shimmer_data.Shimmer_93B6_Gyro_Z_CAL,'DisplayName','Shimmer_Gyro_Z');

    plot(smartphone_data.accX,'DisplayName','smartphone_accX');
    plot(smartphone_data.accY,'DisplayName','smartphone_accY');
    plot(smartphone_data.accZ,'DisplayName','smartphone_accZ');
    plot(smartphone_data.gyrX,'DisplayName','smartphone_gyrX');
    plot(smartphone_data.gyrY,'DisplayName','smartphone_gyrY');
    plot(smartphone_data.gyrZ,'DisplayName','smartphone_gyrZ');

    hold off;
    legend
end

function [edges_sh, edges_sm] = GetSecondEdges(shimmer_data, smartphone_data)
    
%     % Make shimmer_data and smartphone_data of the same size
%     if size(shimmer_data,1) > size(smartphone_data,1)
%         shimmer_data(size(smartphone_data,1):end,:)=[];
%     end

    %% Shimmer second edges
    i_initial_second_sample_sh=1;
    i_final_second_sample_sh=1;    
    edges_sh = [1];

    while i_final_second_sample_sh < size(shimmer_data,1)
        current_second_sh = floor(second(shimmer_data.Shimmer_93B6_Timestamp_Shimmer_CAL(i_initial_second_sample_sh)));
        i_second_sh = current_second_sh;
        i_final_second_sample_sh = i_initial_second_sample_sh;
        while i_second_sh == current_second_sh && i_final_second_sample_sh < size(shimmer_data,1)
            i_final_second_sample_sh = i_final_second_sample_sh + 1;
            i_second_sh = floor(second(shimmer_data.Shimmer_93B6_Timestamp_Shimmer_CAL(i_final_second_sample_sh)));
        end        
        i_initial_second_sample_sh = i_final_second_sample_sh;        
        edges_sh = [edges_sh, i_final_second_sample_sh];    % save edges on array

    end
    
    %% Smartphone second edges
    i_initial_second_sample_sm=1;
    i_final_second_sample_sm=1;
    edges_sm = [1];
    
    while i_final_second_sample_sm < size(smartphone_data,1)
        current_second_sm = floor(second(smartphone_data.Timestamps_ref(i_initial_second_sample_sm)));
        i_second_sm = current_second_sm;
        i_final_second_sample_sm = i_initial_second_sample_sm;
        while i_second_sm == current_second_sm && i_final_second_sample_sm < size(smartphone_data,1)
            i_final_second_sample_sm = i_final_second_sample_sm + 1;
            i_second_sm = floor(second(smartphone_data.Timestamps_ref(i_final_second_sample_sm)));
        end
        %disp([i_initial_second_sample_sm,i_final_second_sample_sm])
        i_initial_second_sample_sm = i_final_second_sample_sm;
        edges_sm = [edges_sm, i_final_second_sample_sm];    % save edges on array
    end
end

function cd = SyncData(edges_sh,edges_sm,shimmer_data,smartphone_data)
    
    % create a table with synchronized data of shimmer and smartphone data
    c = shimmer_data(1:1,:);
    d = smartphone_data(1:1,:);
    cd = [c d];
    cd(:,:) = [];
    
    if size(edges_sh,2)>size(edges_sm,2)
        final_edge = size(edges_sm,2);
    else
        final_edge = size(edges_sh,2);
    end
    
    for i=2:final_edge-1
        a = shimmer_data(edges_sh(i-1):edges_sh(i)-1,:);
        b = smartphone_data(edges_sm(i-1):edges_sm(i)-1,:);
        
        % do shimmer or smartphone data have more samples in this second?
        if size(a, 1) > size(b, 1)
            % if there are more samples of shimmer, create the remaining
            % samples in smartphone data to be equal to the size of shimmer
            % data and fill them with zeros
            b(size(b,1):size(a,1),2:end) = {0};
        else
            % if there are more samples of smarphone, cut smartphone data
            % to get the same size of shimmer data
            b(size(a,1)+1:end,:) = [];
        end
        
        ab = [a b];
        cd = [cd; ab];
    end
    
end

%% Plot data processed
function PlotDataProcessed(cd)
    figure(2)
    hold on;
    %plot(shimmer_data.Shimmer_93B6_Accel_WR_X_CAL,'DisplayName','Shimmer_Accel_X');
    %plot(shimmer_data.Shimmer_93B6_Accel_WR_Y_CAL,'DisplayName','Shimmer_Accel_Y');
    %plot(shimmer_data.Shimmer_93B6_Accel_WR_Z_CAL,'DisplayName','Shimmer_Accel_Z');

    %plot(shimmer_data.Shimmer_93B6_Gyro_X_CAL,'DisplayName','Shimmer_Gyro_X');
    %plot(shimmer_data.Shimmer_93B6_Gyro_Y_CAL,'DisplayName','Shimmer_Gyro_Y');
    plot(cd.Shimmer_93B6_Gyro_Z_CAL,'DisplayName','Shimmer_Gyro_Z');

    plot(cd.accX,'DisplayName','smartphone_accX');
    plot(cd.accY,'DisplayName','smartphone_accY');
    plot(cd.accZ,'DisplayName','smartphone_accZ');
    plot(cd.gyrX,'DisplayName','smartphone_gyrX');
    plot(cd.gyrY,'DisplayName','smartphone_gyrY');
    plot(cd.gyrZ,'DisplayName','smartphone_gyrZ');

    hold off;
    legend
end

%% Plot data processed
function ExportFile(cd, participant_id)

    % Remove unnecessary variables
    cd(:,([8:14,25:34])) = [];
    
    % Export file
    writetable(cd,strcat('dataset_synchronized\',participant_id,'.csv'),'WriteRowNames',true)
end

%% Plot final data
function PlotFinalData(participant_id)
    final_data= readtable(strcat('E:\Google Drive\matlab_slam3\dataset_synchronized\', participant_id));
    plot(final_data.Shimmer_93B6_Gyro_X_CAL,'DisplayName','shimmer_data.Shimmer_93B6_Gyro_X_CAL');hold on;plot(final_data.Shimmer_93B6_Gyro_Y_CAL,'DisplayName','shimmer_data.Shimmer_93B6_Gyro_Y_CAL');plot(final_data.Shimmer_93B6_Gyro_Z_CAL,'DisplayName','shimmer_data.Shimmer_93B6_Gyro_Z_CAL');hold off;
end
