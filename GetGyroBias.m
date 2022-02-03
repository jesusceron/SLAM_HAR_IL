%%% participants id's
%%% december 18: @25282489, @65745876, @76325018, @9335772
%%% december 19: @4610326, @1002819840, @25269362, @34555639, @76306709, @25399881, @25271446, @10290567, @34325893, @34327818, 76323691
%%% december 22: @10534499, @19056006, @41534025, @41331063, @31239627, @1004301404, @14873733, @10526517, @34539568

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Execute this code snippet first

% close all;
% clear all
% participant_id = int2str(31239627);
% file_shimmer_data = strcat(participant_id,'_shimmer_foot.csv');
% shimmer_data= readtable(strcat('E:\PycharmProjects\DataPreProcess\dataset\', file_shimmer_data));
% 
% acc_complete = table2array(shimmer_data(:,2:4));
% gyro_complete = table2array(shimmer_data(:,5:7));
% 
% %%% look for a standing still moment just before participant started to move
% figure(1)
% hold on;
% plot(gyro_complete(:,1),'DisplayName','Shimmer_Gyro_X');
% plot(gyro_complete(:,2),'DisplayName','Shimmer_Gyro_Y');
% plot(gyro_complete(:,3),'DisplayName','Shimmer_Gyro_Z');
% hold off;
% legend

%% After obtaining the range, execute this lines
w = 103200 : 105880; %change this range according to the plot.

acc = zeros(length(acc_complete),3);
acc(:,1) = acc_complete(:,2);
acc(:,2) = acc_complete(:,3); % y = z
acc(:,3) = acc_complete(:,1);

gyr = zeros(length(gyro_complete),3);
gyr(:,1) = gyro_complete(:,2);
gyr(:,2) = gyro_complete(:,3);
gyr(:,3) = gyro_complete(:,1);
gyr = deg2rad(gyr);

gravity=mean(sqrt(acc(w,1).^2 + acc(w,2).^2 + acc(w,3).^2));  % m/s^2
acc_mean=[mean(acc(w,1)),mean(acc(w,2)),mean(acc(w,3))];
bias_gyr=[mean(gyr(w, 1)), mean(gyr(w, 2)), mean(gyr(w, 3))];

result = [gravity,0,0;acc_mean;bias_gyr];
clearvars -except result