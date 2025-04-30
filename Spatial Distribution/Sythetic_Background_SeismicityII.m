clear all;clc; close all;
file_path = 'WP_50k_MS_20240807.csv';
dataTable = readtable(file_path);
write_matrix = table2array(dataTable);
MS_lon = write_matrix(1,1);
MS_lat = write_matrix(1,2);
MS_dep = write_matrix(1,7)*1000;
MS_time = datetime(write_matrix(1, [3, 4, 5, 8, 9, 10]));

mag = write_matrix(:,6);

index = find(mag >= 3);
file_path = 'WP_50k_from2000_beforeMS.csv';
% Mc = 0; % no Mc cutoff

dataTable = [];
write_matrix = [];
dataTable = readtable(file_path);
write_matrix = table2array(dataTable);
% validindex = find(write_matrix(:,6)>= Mc);
% write_matrix = write_matrix(validindex,:);

depth = write_matrix(:,7);
validindex = find(depth>=0);
write_matrix = write_matrix(validindex,:);

validindex = find(datetime(write_matrix(:, [3, 4, 5, 8, 9, 10])) >= datetime(2000,01,01));
write_matrix = write_matrix(validindex,:);

lon = write_matrix(:,1);
lat = write_matrix(:,2);
mag = write_matrix(:,6);
depth = write_matrix(:,7);

proj = projcrs(7855);
proj.GeographicCRS.Name;


figure(1); % the distribution of pre-mainshock seismicity in depth
subplot(2,1,1);
a = histogram(depth);
xlabel('Depth (KM)','FontSize',14);
ylabel('Frequency','FontSize',14);
title('Histogram (Catalog)','FontSize',14);
grid on;


subplot(2,1,2); % the normal distribution of pre-mainshock seismicity in depth
x = sort(depth);
pdf_values = normpdf(x, mean(x), std(x));
plot(x,pdf_values,'linewidth',2);
grid on;

%% 10,000 points (synthetic background seismicity II contains 10,000 events)
% Parameters
n = 10000; % Number of random points to generate
radius = 50; % Radius of the circle --50km

% Generate random angles uniformly distributed between 0 and 2*pi
theta = 2 * pi * rand(1, n);

% Generate random radii with a uniform distribution in the area
r = radius * sqrt(rand(1, n));
surface_distance = r';

depth_distance = normrnd(mean(x), std(x), [1, n]);
index_depth = find(depth_distance >= 0);
while length(index_depth) < n
   depth_distance = [depth_distance, normrnd(mean(x), std(x), [1, 20])];
   index_depth = find(depth_distance >= 0);
end
depth_distance = depth_distance(index_depth(1:n));

pdf_values = normpdf(sort(depth_distance'), mean(depth_distance'), std(depth_distance'));
hold on;
plot(sort(depth_distance'),pdf_values,'r:','LineWidth',2); % the normal distribution of Synthetic Background Seismcity II in depth

xlabel('Depth (KM)','FontSize',14);
ylabel('Density','FontSize',14);

legend('Catalog','Synthetic Catalog');
title('Probability Density Distribution','FontSize',14);

% calculate the time gap
currentDateTime = datetime(write_matrix(:, [3, 4, 5, 8, 9, 10]));
Time_beforeMS(:,1) = days(MS_time-currentDateTime(:));


% calculate the distance gap
[MS_Lon_meter,MS_Lat_meter] = projfwd(proj,MS_lat,MS_lon);


[Lon_meter,Lat_meter] = projfwd(proj,lat,lon);
Dep = depth*1000;

Distance_beforeMS(:,1) = sqrt(power(Lon_meter(:)-MS_Lon_meter,2)+power(Lat_meter(:)-MS_Lat_meter,2)+power(Dep(:)-MS_dep,2)); 

figure(2);
subplot(2,2,1); % observed seismicity spatial distribution
scatter(currentDateTime,Distance_beforeMS/1000,'filled');
xlabel('Year','FontSize',14);
ylabel('Distance(KM)','FontSize',14);
grid on;title('Observed Catalog','FontSize',14)
% xscale log

subplot(2,2,2); % compariosn of spatial distribution between observed seismicity and the sythetic catalogue
x = sort(Distance_beforeMS)/1000;
pdf_values = normpdf(x, mean(x), std(x));
pre_ms_pdf_withoutMc = [x,pdf_values];

plot(pre_ms_pdf_withoutMc(:,1),pre_ms_pdf_withoutMc(:,2),'linewidth',2);


distance = sqrt(power(depth_distance' - MS_dep/1000,2)+power(surface_distance,2));
hold on;
x = sort(distance);
pdf_values = normpdf(x, mean(x), std(x));
plot(x,pdf_values,'linewidth',2);
ylabel('Density','FontSize',14);
xlabel('KM','FontSize',14)

legend('Observed','Synthetic');
grid on;
title('Probability Density Function','fontsize',14);

subplot(2,2,3);
histogram(pre_ms_pdf_withoutMc(:,1),'BinWidth',1);
title('Observed Catalog (Distance)','FontSize',14);
grid on;
xlabel('Distance (KM)','FontSize',14);
ylabel('Frequency','FontSize',14);

subplot(2,2,4);
histogram(distance,'BinWidth',1);
title('Synthetic Catalog (Distance)','FontSize',14);
grid on;
xlabel('Distance (KM)','FontSize',14);
ylabel('Frequency','FontSize',14);