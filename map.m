% Map shown in the "Trajectory ZUPT-KF" figure

% Map shown in the main UI
% figure
grid on

xlim([0,20])
ylim([1, 9])
hold on
xlabel('x(meters)')
ylabel('y(meters)')

% Whole the seminar room
rectangle('Position',[5.5 2 12.6 6.6])

% %stairs
rectangle('Position',[2.2 2.2500 1.5 0.3125])
rectangle('Position',[2.2 2.5625 1.5 0.3125])
rectangle('Position',[2.2 2.8750 1.5 0.3125])
rectangle('Position',[2.2 3.1875 1.5 0.3125])
rectangle('Position',[2.2 3.5000 1.5 0.3125])
rectangle('Position',[2.2 3.8125 1.5 0.3125])
rectangle('Position',[2.2 4.1250 1.5 0.3125])
rectangle('Position',[2.2 4.4375 1.5 0.3125])
rectangle('Position',[2.2 4.7500 1.5 0.3125])
rectangle('Position',[2.2 5.0625 1.5 0.3125])
rectangle('Position',[2.2 5.3750 1.5 0.3125])
rectangle('Position',[2.2 5.6875 1.5 0.3125])
rectangle('Position',[2.2 6 1.5 0.3125])
rectangle('Position',[2.2 6.3125 1.5 0.3125])
rectangle('Position',[2.2 6.6250 1.5 0.3125])

%zones without use:
rectangle('Position',[5.5 2 4.3 3], 'FaceColor','k')
rectangle('Position',[5.5 6.4 1.5 2.2], 'FaceColor','k')

% Room and living room
rectangle('Position',[11.8 2 0.8 1.8], 'FaceColor','k') % dining table
rectangle('Position',[13.9 2 0.6 1.4], 'FaceColor','k') % wall btw dining and living room
rectangle('Position',[15.1 2 1.4 0.6], 'FaceColor','k') % Living table
rectangle('Position',[13.9 4.6 4.2 0.4], 'FaceColor','k') % wall btw room and living room
rectangle('Position',[16.5 6 1.6 1.6], 'FaceColor','k') % bed

% Kitchen
rectangle('Position',[8.1 6 3.1 0.6], 'FaceColor','k') % Wall btw bathroom-kitchen and main corridor
rectangle('Position',[12.5 6 1.4 0.6], 'FaceColor','k') % Wall btw kitchen and corridor
rectangle('Position',[9.5 6.6 0.8 2], 'FaceColor','k') % Wall btw bathroom and kitchen
rectangle('Position',[13.1 6.6 0.8 2], 'FaceColor','k') % Wall btw kitchen and room
rectangle('Position',[10.3 8 2.8 0.6], 'FaceColor','k') % Tables of the kitchen

%% BEACONS MAP
% initial/end point 
plot(3, 7.5,'p','markerfacecolor','g','MarkerSize',8,'DisplayName','Start / End point')

% Plot fixed Beacons
plot(16.6, 6.4,'^','markerfacecolor',[0.8500, 0.3250, 0.0980],'color',[0.8500, 0.3250, 0.0980],'MarkerSize',4,'DisplayName','Door')  % Room
%viscircles([13, 4.4], 4,'Color','y','LineStyle','--');
plot(11.7, 8,'^','markerfacecolor',[0.9290, 0.6940, 0.1250],'color',[0.9290, 0.6940, 0.1250],'MarkerSize',4,'DisplayName','Kitchen')  % Kitchen
%viscircles([8, 6.2], 4,'Color','k','LineStyle','--');
plot(8.1, 8.6,'^','markerfacecolor','b','color','b','MarkerSize',4,'DisplayName','Bathroom')  % Bathroom
%viscircles([4.5, 6.8], 4,'Color',[0.4940, 0.1840, 0.5560],'LineStyle','--');
plot(12.2, 3,'^','markerfacecolor',[0.75, 0.75, 0],'color',	[0.75, 0.75, 0],'MarkerSize',4,'DisplayName','Dining Table')  % Dining table
%viscircles([8.4, 1.3], 4,'Color',[0.3010, 0.7450, 0.9330],'LineStyle','--');
plot(15.8, 2.5,'^','markerfacecolor',[0, 0.75, 0.75],'color',[0, 0.75, 0.75],'MarkerSize',4,'DisplayName','Living room')  % Living room
%viscircles([12.3, 1], 4,'Color',[0.4660, 0.6740, 0.1880],'LineStyle','--');


% r = plot(fixed_beacons_results(:,1),fixed_beacons_results(:,2),'.','MarkerSize',10,'color','b'); % Bathroom
% r.Annotation.LegendInformation.IconDisplayStyle = 'off';
% k = plot(fixed_beacons_results(:,3),fixed_beacons_results(:,4),'.','MarkerSize',10,'color',[0.9290, 0.6940, 0.1250]); % Kitchen
% k.Annotation.LegendInformation.IconDisplayStyle = 'off';
% b=plot(fixed_beacons_results(:,5),fixed_beacons_results(:,6),'.','MarkerSize',10,'color',[0.8500, 0.3250, 0.0980]); % Room
% b.Annotation.LegendInformation.IconDisplayStyle = 'off';
% d=plot(fixed_beacons_results(:,7),fixed_beacons_results(:,8),'.','MarkerSize',10,'color',[0, 0.75, 0.75]); % Living room
% d.Annotation.LegendInformation.IconDisplayStyle = 'off';
% l=plot(fixed_beacons_results(:,9),fixed_beacons_results(:,10),'.','MarkerSize',10,'color',[0.75, 0.75, 0]); % Dining room
% l.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Plot 'activity' Beacons
plot(5.5, 5.5,'s','markerfacecolor',[0.4940, 0.1840, 0.5560],'color',[0.4940, 0.1840, 0.5560],'MarkerSize',5,'DisplayName','Door')  % Door
plot(9.1, 8.3,'s','markerfacecolor','m','color','m','MarkerSize',5,'DisplayName','Toilet')  % Toiler lid
plot(7.3, 8,'s','markerfacecolor',[0.4660, 0.6740, 0.1880],'color',[0.4660, 0.6740, 0.1880],'MarkerSize',5,'DisplayName','Broom')  % Broom
plot(13.5, 7.4,'s','markerfacecolor','r','color','r','MarkerSize',5,'DisplayName','Pitcher')  % Pitcher
plot(8.7, 6.3,'s','markerfacecolor','y','color','y','MarkerSize',5,'DisplayName','Hair brush')  % Hair brush

% d=plot(active_beacons_results(:,1),active_beacons_results(:,2),'.','MarkerSize',10,'markerfacecolor',[0.4940, 0.1840, 0.5560],'color',[0.4940, 0.1840, 0.5560]); % Door
% d.Annotation.LegendInformation.IconDisplayStyle = 'off';
% br=plot(active_beacons_results(:,3),active_beacons_results(:,4),'.','MarkerSize',10,'markerfacecolor','y','color','y'); % Brush
% br.Annotation.LegendInformation.IconDisplayStyle = 'off';
% bro=plot(active_beacons_results(:,5),active_beacons_results(:,6),'.','MarkerSize',10,'markerfacecolor',[0.4660, 0.6740, 0.1880],'color',[0.4660, 0.6740, 0.1880]); % Broom
% bro.Annotation.LegendInformation.IconDisplayStyle = 'off';
% t=plot(active_beacons_results(:,7),active_beacons_results(:,8),'.','MarkerSize',10,'markerfacecolor','m','color','m'); % Toilet
% t.Annotation.LegendInformation.IconDisplayStyle = 'off';
% p=plot(active_beacons_results(:,9),active_beacons_results(:,10),'.','MarkerSize',10,'markerfacecolor','r','color','r'); % Pitcher
% p.Annotation.LegendInformation.IconDisplayStyle = 'off';
% legend
