% Map shown in the "Trajectory ZUPT-KF" figure

% Map shown in the main UI
% figure
grid on
%title('Trajectory (ZUPT KF)');
axis equal
xlim([1, 22])
ylim([1, 13])
hold on
xlabel('x(m)')
ylabel('y(m)')


%rectangle('Position',[2 3.2 1.03 1.8])
%rectangle('Position',[2 2.2 1.03 1.01])
% rectangle('Position',[3.03 2.2 0.3755 1.01])
% rectangle('Position',[3.4050 2.2 0.3755 1.01])
% rectangle('Position',[3.7810 2.2 0.3755 1.01])
% rectangle('Position',[4.1565 2.2 0.3755 1.01])
% rectangle('Position',[4.5320 2.2 0.3755 1.01])
% rectangle('Position',[4.9075 2.2 0.3755 1.01])
% rectangle('Position',[5.2830 2.2 0.3755 1.01])
% rectangle('Position',[5.6585 2.2 0.3755 1.01])
% rectangle('Position',[6.0340 2.2 0.3755 1.01])
% rectangle('Position',[6.4095 2.2 0.3755 1.01])
% rectangle('Position',[6.7850 2.2 0.3755 1.01])
%rectangle('Position',[7.1605 2.2 0.99 1.01])



%% second floor
x = [2 3.8];
y = [8 8];
pl = line(x,y);
pl.Color = 'black';

x = [2 2];
y = [5 8];
pl = line(x,y);
pl.Color = 'black';

x = [3.8 3.8];
y = [5 8];
pl = line(x,y);
pl.Color = 'black';

x = [2 3];
y = [5 5];
pl = line(x,y);
pl.Color = 'black';

x = [2 2];
y = [5 2.2];
pl = line(x,y);
pl.Color = 'black';

x = [2 3.03];
y = [2.2 2.2];
pl = line(x,y);
pl.Color = 'black';

%% stairs
rectangle('Position',[3.03 2.2 0.318182 1.01])
rectangle('Position',[3.3482 2.2 0.318182 1.01])
rectangle('Position',[3.6664 2.2 0.318182 1.01])
rectangle('Position',[3.9845 2.2 0.318182 1.01])
rectangle('Position',[4.3027 2.2 0.318182 1.01])
rectangle('Position',[4.6209 2.2 0.318182 1.01])
rectangle('Position',[4.9391 2.2 0.318182 1.01])
rectangle('Position',[5.2573 2.2 0.318182 1.01])
rectangle('Position',[5.5755 2.2 0.318182 1.01])
rectangle('Position',[5.8939 2.2 0.318182 1.01])
rectangle('Position',[6.2118 2.2 0.318182 1.01])

x = [6.5300 7.53];
y = [2.2 2.2];
pl = line(x,y);
pl.Color = 'black';

x = [7.53 7.53];
y = [2.2 4];
pl = line(x,y);
pl.Color = 'black';


%% hallway
x = [3.8 7.53];
y = [5 5];
pl = line(x,y);
pl.Color = 'black';

x = [7.53 12.31]; %12.31-7.53 = 4.78 hallway long
y = [5 5];
pl = line(x,y);
pl.Color = 'black';

x = [7.53 14.13]; 
y = [4 4];
pl = line(x,y);
pl.Color = 'black';

%% Toilet
x = [13.58 13.58];
y = [4 3];
pl = line(x,y);
pl.Color = 'black';

x = [13.58 15.58];
y = [3 3];
pl = line(x,y);
pl.Color = 'black';

x = [15.58 15.58];
y = [3 4];
pl = line(x,y);
pl.Color = 'black';

x = [15.58 14.93];
y = [4 4];
pl = line(x,y);
pl.Color = 'black';

% Main door
x = [15.58 16.64];
y = [4 4];
pl = line(x,y);
pl.Color = 'black';

x = [16.64 16.64];
y = [4 4.8];
pl = line(x,y);
pl.Color = 'black';

x = [16.64 16.64];
y = [5.9 6.82];
pl = line(x,y);
pl.Color = 'black';

% Living room - kitchen
x = [16.64 21.67];
y = [6.82 8.9273];
pl = line(x,y);
pl.Color = 'black';

x = [12.31 10.722];
y = [5 8.781];
pl = line(x,y);
pl.Color = 'black';

x = [21.67 20.08];
y = [8.9273 12.713];
pl = line(x,y);
pl.Color = 'black';

x = [10.722 20.08];
y = [8.781 12.713];
pl = line(x,y);
pl.Color = 'black';

% Plot 'activity' Beacons
plot(3.3, 5,'s','markerfacecolor','r','color','r','MarkerSize',7,'DisplayName','Door')  % Door
plot(21, 9.5,'s','markerfacecolor','r','color','r','MarkerSize',7,'DisplayName','Toilet')  % Pitcher
plot(13.8, 3.5,'s','markerfacecolor','r','color','r','MarkerSize',7,'DisplayName','Broom')  % Water tap
plot(19.64,4.3,'s','markerfacecolor','r','color','r','MarkerSize',7,'DisplayName','Pitcher')  % Pitcher in garden
plot(16.3, 6.3,'s','markerfacecolor','r','color','r','MarkerSize',7,'DisplayName','Hair brush')  % Broom
plot(15, 3.5,'s','markerfacecolor','r','color','r','MarkerSize',7,'DisplayName','Toilet')  % Toilet
   
% % Plot fixed Beacons
plot(10.03, 4,'^','markerfacecolor','b','color','b','MarkerSize',4,'DisplayName','Door')  % Corridor
%viscircles([13, 4.4], 4,'Color','y','LineStyle','--');
plot(20.12, 11.1,'^','markerfacecolor','b','color','b','MarkerSize',4,'DisplayName','Kitchen')  % Kitchen
%viscircles([8, 6.2], 4,'Color','k','LineStyle','--');
plot(15, 8.3,'^','markerfacecolor','b','color','b','MarkerSize',4,'DisplayName','Chair')  % Chair
%viscircles([4.5, 6.8], 4,'Color',[0.4940, 0.1840, 0.5560],'LineStyle','--');
plot(17.3, 9.8,'^','markerfacecolor','b','color','b','MarkerSize',4,'DisplayName','Dining Table')  % Dining table
%viscircles([8.4, 1.3], 4,'Color',[0.3010, 0.7450, 0.9330],'LineStyle','--');
plot(13, 7.5585,'^','markerfacecolor','b','color','b','MarkerSize',4,'DisplayName','Living room')  % Living room
%viscircles([12.3, 1], 4,'Color',[0.4660, 0.6740, 0.1880],'LineStyle','--');

%%% linea perpendicular para calculo de posición de beacon puesto en la
%%% mesa
% x = [18.35 16];
% y = [7.5384 11.6098];
% pl = line(x,y);
% pl.Color = 'red';

%%% linea perpendicular para calculo de posición de beacon puesto en la
%%% sala
% x = [11.5 20];
% y = [6.9285 10.4985];
% pl = line(x,y);
% pl.Color = 'red';

% %% BEACONS MAP
% % initial/end point 
% plot(3, 7.5,'p','markerfacecolor','g','MarkerSize',8,'DisplayName','Start / End point')
% 
% % Plot fixed Beacons
% plot(16.6, 6.4,'^','markerfacecolor',[0.8500, 0.3250, 0.0980],'color',[0.8500, 0.3250, 0.0980],'MarkerSize',4,'DisplayName','Door')  % Room
% %viscircles([13, 4.4], 4,'Color','y','LineStyle','--');
% plot(11.7, 8,'^','markerfacecolor',[0.9290, 0.6940, 0.1250],'color',[0.9290, 0.6940, 0.1250],'MarkerSize',4,'DisplayName','Kitchen')  % Kitchen
% %viscircles([8, 6.2], 4,'Color','k','LineStyle','--');
% plot(8.1, 8.6,'^','markerfacecolor','b','color','b','MarkerSize',4,'DisplayName','Bathroom')  % Bathroom
% %viscircles([4.5, 6.8], 4,'Color',[0.4940, 0.1840, 0.5560],'LineStyle','--');
% plot(12.2, 3,'^','markerfacecolor',[0.75, 0.75, 0],'color',	[0.75, 0.75, 0],'MarkerSize',4,'DisplayName','Dining Table')  % Dining table
% %viscircles([8.4, 1.3], 4,'Color',[0.3010, 0.7450, 0.9330],'LineStyle','--');
% plot(15.8, 2.5,'^','markerfacecolor',[0, 0.75, 0.75],'color',[0, 0.75, 0.75],'MarkerSize',4,'DisplayName','Living room')  % Living room
% %viscircles([12.3, 1], 4,'Color',[0.4660, 0.6740, 0.1880],'LineStyle','--');
% 
% 
% % r = plot(fixed_beacons_results(:,1),fixed_beacons_results(:,2),'.','MarkerSize',10,'color','b'); % Bathroom
% % r.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % k = plot(fixed_beacons_results(:,3),fixed_beacons_results(:,4),'.','MarkerSize',10,'color',[0.9290, 0.6940, 0.1250]); % Kitchen
% % k.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % b=plot(fixed_beacons_results(:,5),fixed_beacons_results(:,6),'.','MarkerSize',10,'color',[0.8500, 0.3250, 0.0980]); % Room
% % b.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % d=plot(fixed_beacons_results(:,7),fixed_beacons_results(:,8),'.','MarkerSize',10,'color',[0, 0.75, 0.75]); % Living room
% % d.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % l=plot(fixed_beacons_results(:,9),fixed_beacons_results(:,10),'.','MarkerSize',10,'color',[0.75, 0.75, 0]); % Dining room
% % l.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% % Plot 'activity' Beacons
% plot(5.5, 5.5,'s','markerfacecolor',[0.4940, 0.1840, 0.5560],'color',[0.4940, 0.1840, 0.5560],'MarkerSize',5,'DisplayName','Door')  % Door
% plot(9.1, 8.3,'s','markerfacecolor','m','color','m','MarkerSize',5,'DisplayName','Toilet')  % Toiler lid
% plot(7.3, 8,'s','markerfacecolor',[0.4660, 0.6740, 0.1880],'color',[0.4660, 0.6740, 0.1880],'MarkerSize',5,'DisplayName','Broom')  % Broom
% plot(13.5, 7.4,'s','markerfacecolor','r','color','r','MarkerSize',5,'DisplayName','Pitcher')  % Pitcher
% plot(8.7, 6.3,'s','markerfacecolor','y','color','y','MarkerSize',5,'DisplayName','Hair brush')  % Hair brush
% 
% % d=plot(active_beacons_results(:,1),active_beacons_results(:,2),'.','MarkerSize',10,'markerfacecolor',[0.4940, 0.1840, 0.5560],'color',[0.4940, 0.1840, 0.5560]); % Door
% % d.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % br=plot(active_beacons_results(:,3),active_beacons_results(:,4),'.','MarkerSize',10,'markerfacecolor','y','color','y'); % Brush
% % br.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % bro=plot(active_beacons_results(:,5),active_beacons_results(:,6),'.','MarkerSize',10,'markerfacecolor',[0.4660, 0.6740, 0.1880],'color',[0.4660, 0.6740, 0.1880]); % Broom
% % bro.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % t=plot(active_beacons_results(:,7),active_beacons_results(:,8),'.','MarkerSize',10,'markerfacecolor','m','color','m'); % Toilet
% % t.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % p=plot(active_beacons_results(:,9),active_beacons_results(:,10),'.','MarkerSize',10,'markerfacecolor','r','color','r'); % Pitcher
% % p.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % legend
