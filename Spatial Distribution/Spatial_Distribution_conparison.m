clear all;clc; close all;
file_path = 'WP_50k_MS_20240807.csv';

Mc = 0.8;
dataTable = readtable(file_path);
write_matrix = table2array(dataTable);
validindex = find(write_matrix(:,6)>= Mc);
write_matrix = write_matrix(validindex,:);
lon = write_matrix(:,1);
lat = write_matrix(:,2);
% Year = write_matrix(:,3);
% Month = write_matrix(:,4);
% Day = write_matrix(:,5);
mag = write_matrix(:,6);
depth = write_matrix(:,7);
% Hour = write_matrix(:,8);
% Minute = write_matrix(:,9);
% Second = write_matrix(:,10);


% calculate the elapsed day
currentDateTime = datetime(write_matrix(:, [3, 4, 5, 8, 9, 10]));
Time_sinceMS(:,1) = days(currentDateTime(:)-currentDateTime(1));

%% Spatial Distribution
% calculate the distance gap
proj = projcrs(7855);
proj.GeographicCRS.Name;
[Lon_meter,Lat_meter] = projfwd(proj,lat,lon);
Dep = depth*1000;
Distance_sinceMS(:,1) = sqrt(power(Lon_meter(:)-Lon_meter(1),2)+power(Lat_meter(:)-Lat_meter(1),2)+power(Dep(:)-Dep(1),2)); 

%% Time slice
%% pdf
figure;
% initial 30 days
index_0_30 = find(Time_sinceMS >= 0 & Time_sinceMS < 30);
x = sort(Distance_sinceMS(index_0_30))/1000;
hold on;
pdf_values = normpdf(x, mean(x), std(x));
% pdf_values = gampdf(x, mean(x));
plot(x,pdf_values, 'LineWidth', 2);
disp(sprintf('initial 30 days, distance from MS hypocenter: %.2f +/- %.2f', ...
    mean(x), std(x)))


% 30 days to 90 days
index_30_90 = find(Time_sinceMS >= 30 & Time_sinceMS < 90);
x = sort(Distance_sinceMS(index_30_90))/1000;
hold on;
pdf_values = normpdf(x, mean(x), std(x));
% pdf_values = gampdf(x, mean(x));
plot(x,pdf_values, 'LineWidth', 2);
disp(sprintf('30 days to 90 days, distance from MS hypocenter: %.2f +/- %.2f', ...
    mean(x), std(x)))


% 90 days to 500 days
index_90_500 = find(Time_sinceMS >= 90 & Time_sinceMS < 500);
x = sort(Distance_sinceMS(index_90_500))/1000;
hold on;
pdf_values = normpdf(x, mean(x), std(x));
% pdf_values = gampdf(x, mean(x));
plot(x,pdf_values, 'LineWidth', 2);
disp(sprintf('90 days to 500 days, distance from MS hypocenter: %.2f +/- %.2f', ...
    mean(x), std(x)))


% more than 500 days
index_morethan_500 = find(Time_sinceMS >= 500);
x = sort(Distance_sinceMS(index_morethan_500))/1000;
hold on;
pdf_values = normpdf(x, mean(x), std(x));
% pdf_values = gampdf(x, mean(x));
plot(x,pdf_values, 'LineWidth', 2);
disp(sprintf('more than 500 days, distance from MS hypocenter: %.2f +/- %.2f', ...
    mean(x), std(x)))

hold on;

hold on;
load('pre_ms_pdf_withoutMc.mat');
plot(pre_ms_pdf_withoutMc(:,1), pre_ms_pdf_withoutMc(:,2),':', 'LineWidth', 2);
x = sort(pre_ms_pdf_withoutMc(:,1));
disp(sprintf('before MS, distance from MS hypocenter: %.2f +/- %.2f', ...
    mean(x), std(x)))


legend({'Initial 30 Days (elapsed)','30 to 90 days (elapsed)','90 to 500 days (elapsed)','Since 500 days (elapsed)','Background Seismicity (synthetic)'},'fontsize',12);
ylabel('Density','FontSize',14);xlabel('Distance from Mainshock (KM)','FontSize',14);grid on;

%% Function
function [x_max_curvature,y_max_curvature] = max_curvature(x,y)
    % Calculate the first derivatives
    dx = gradient(x);
    dy = gradient(y);
    
    % Calculate the second derivatives
    d2x = gradient(dx);
    d2y = gradient(dy);
    
    % Calculate curvature
    curvature = abs(d2y .* dx - d2x .* dy) ./ (dx.^2 + dy.^2).^(3/2);
    
    % Find the maximum curvature and its index
    [max_curvature_value, idx] = max(curvature);
    
    % Display the maximum curvature and the corresponding point
    x_max_curvature = x(idx)
    y_max_curvature = y(idx)
end



