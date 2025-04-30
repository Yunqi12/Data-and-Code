% calculate the seismic moment (energy) level per month after mainshock
clc;
clear all;close all;
file_path = 'WP_50k_MS_20240807.csv';
disp(['CSV file name: ', file_path]);
dataTable = readtable(file_path);
write_matrix = table2array(dataTable);
currentDateTime = datetime(write_matrix(:, [3, 4, 5, 8, 9, 10]));
mag = write_matrix(:, 6);
Energy = zeros(2,1);

for i = 1:length(mag)
    Mw(i) = Ml2Mw(mag(i));
end
Mw(1) = 5.9;
Energy = power(10,3/2.*(Mw+10.7));

startdate = datetime(2021,9,1);
enddate = startdate+calmonths(1);
count1 = 1;
while startdate <= currentDateTime(end)
    index = find(currentDateTime >= startdate & currentDateTime < enddate);
    Event_month_afterMS(count1,1) = year(startdate); Event_month_afterMS(count1,2) = month(startdate);
    Event_month_afterMS(count1,3) = sum(Energy(index));
    startdate = enddate;
    enddate = startdate+calmonths(1);
    count1 = count1 +1;
end

%% Normal Distribution -- background
subplot(2,1,1);
load('Energy_month_beforeMS.mat'); % has been calculated by using Energy_monthly_beforeMS code
index = find(Event_month(:,3) ~= 0);
x = log10(sort(Event_month(index,3)));
hold on;
pdf_values = normpdf(x, mean(x), std(x));
plot(power(10,x),pdf_values, 'LineWidth', 2);
xscale log;
xlabel('Cumulative Seismic Moment / Month','fontsize',14);
ylabel('Density','fontsize',14);
title('Pre-Mainshock','fontsize',14);grid on;


%% compare background with current seismicity energy released level
subplot(2,1,2);
dates = datetime(Event_month_afterMS(:,1), Event_month_afterMS(:,2), 1);
time_lower_bound = min(currentDateTime)-days(40); % the axis limit for time
time_upper_bound = max(currentDateTime);
% 2 standard deviation of pre-mainshock background seismicity level
up = power(10,mean(x)+2*std(x));
bottom = power(10,mean(x)-2*std(x));
hold on;
x_fill = [time_lower_bound, time_upper_bound, time_upper_bound, time_lower_bound];
y_fill = [bottom, bottom, up, up];
% Fill the color between the lines
color = [0.8 0.8 0.8];
fill(x_fill, y_fill, color, 'FaceAlpha', 0.5,'EdgeColor','none'); % Blue color with 50% transparency
% 1 standard deviation of pre-mainshock background seismicity level
hold on;
up = power(10,mean(x)+std(x));
bottom = power(10,mean(x)-std(x));
x_fill = [time_lower_bound, time_upper_bound, time_upper_bound, time_lower_bound];
y_fill = [bottom, bottom, up, up];
color = [0.6 0.6 0.6];
fill(x_fill, y_fill, color, 'FaceAlpha', 0.5,'EdgeColor','none'); % Blue color with 50% transparency
hold on;
xlim([time_lower_bound time_upper_bound]);
stem(dates,Event_month_afterMS(:,3),'LineStyle','-.',...
     'MarkerFaceColor','black',...
     'MarkerEdgeColor','black'); 
yscale log;
xlabel('Time','FontSize',14);
ylabel('Cumulative Seismic Moment / Month','FontSize',14);
title('Post-Mainshock','FontSize',14)
grid on;
% legend('Energy Released','Cumulative Events above Mc','Cumulative Events');
legend('2 Std','1 Std','Energy Released','Location','northeast');
