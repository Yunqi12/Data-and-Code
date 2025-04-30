clc;close all;
clear all;

%% before the MS -- ML
file_path = 'WP_50k_from2000_beforeMS.csv';
disp(['CSV file name: ', file_path]);
dataTable = readtable(file_path);
write_matrix = table2array(dataTable);

startdate=datetime('2000-01-01 00:00:00');%year,month,day,hour,minute,second
enddate=datetime('2031-01-01 00:00:00');%year,month,day,hour,minute,second
currentDateTime = datetime(write_matrix(:, [3, 4, 5, 8, 9, 10]));
index = find(currentDateTime >= startdate & currentDateTime < enddate);
write_matrix = write_matrix(index,:);
Ml = write_matrix(:,6);

%% GR law calculation
magnitudeVector=write_matrix(:,6);
%divide in bins
mag_range=min(magnitudeVector):0.1:7;
[number,bin]= hist(magnitudeVector,mag_range'); % calculate the number per bin
mag_count_pre=[bin,number'];
for i=1:length(mag_count_pre)
    mag_count_pre(i,3)=sum(mag_count_pre(i:end,2));
end
%% Mc MAXC (Max Curvature) method
mCatalog=write_matrix;
nMethod =1; % MAX CURVE OR nMethod=5; % BEST COMBINATION OR nMethod=6; % EMR
fBinning = 0.1; % size of our magnitude bins
fMcCorrection = 0; % Used to add a correction because maximum curvature often underestimates Mc
fMc = 2.03;

subplot(1,2,1)
vSel = mCatalog(:,6) >= fMc;
mCat = mCatalog(vSel,:);
fBinning=0.1;
[fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
a1_pre=fAvalue;b1_pre=fBvalue;
disp(sprintf('Maximum Likelihood Estimation: b=%.2f+/-%.2f,a=%.2f %', b1_pre,fStdDev,a1_pre));

semilogy(mag_count_pre(:,1),mag_count_pre(:,3),'k^','Markersize',6) % cumulative events%
hold on;
Mc=fMc;
Max_GR=5.8;
x=fMc:0.1:Max_GR;
semilogy(x,power(10,a1_pre-b1_pre*x),'k','Linewidth',2); % GR law
grid on;

xlabel('Magnitude','fontsize',14);
ylabel('Frequency','fontsize',14);

%% before the MS -- MW
file_path = 'WP_50k_from2000_beforeMS.csv';
dataTable = readtable(file_path);
write_matrix = table2array(dataTable);
startdate=datetime('2000-01-01 00:00:00');%year,month,day,hour,minute,second
enddate=datetime('2031-01-01 00:00:00');%year,month,day,hour,minute,second
currentDateTime = datetime(write_matrix(:, [3, 4, 5, 8, 9, 10]));
index = find(currentDateTime >= startdate & currentDateTime < enddate);
write_matrix = write_matrix(index,:);

%% Convert Ml to Mw
Ml = write_matrix(:,6);

Mw = zeros(2,1);
for i = 1:length(Ml)
    Mw(i) = Ml2Mw(Ml(i));
end
Mw = round(Mw,1);
write_matrix(:,6) = Mw;

%% GR law calculation
magnitudeVector=write_matrix(:,6);
%divide in bins
mag_range=min(magnitudeVector):0.1:7;
[number,bin]= hist(magnitudeVector,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end

mCatalog=write_matrix;
nMethod =1; % MAX CURVE OR nMethod=5; % BEST COMBINATION OR nMethod=6; % EMR
fBinning = 0.1; % size of our magnitude bins
fMcCorrection = 0; % Used to add a correction because maximum curvature often underestimates Mc
fMc = 2.3;

hold on;
vSel = mCatalog(:,6) >= fMc;
mCat = mCatalog(vSel,:);
fBinning=0.1;
[fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
a1=fAvalue;b1=fBvalue;
disp(sprintf('Maximum Likelihood Estimation: b=%.2f+/-%.2f,a=%.2f %', b1,fStdDev,a1));
semilogy(mag_count(:,1),mag_count(:,3),'k.','Markersize',15) % cumulative events%

hold on;
Mc=fMc;
Max_GR=5.8;
x=fMc:0.1:Max_GR;
semilogy(x,power(10,a1-b1*x),'k--','Linewidth',2); % GR law
grid on;
title('Background Seismicity');
set(gca,'fontsize',12);

legend({'M_L Catalogue','M_L Catalogue (Regression Line)','M_W Catalogue','M_W Catalogue (Regression Line)'},'FontSize',12,'Location','southwest');

subplot(1,2,2)
%% WPAS ML
file_path = 'WP_50k_MS_20240807.csv';
% file_path = 'WP_50k_afterMS.csv';
dataTable = readtable(file_path);
write_matrix = table2array(dataTable);
startdate=datetime('2021-09-21 23:15:53');%year,month,day,hour,minute,second
enddate=datetime('2031-01-01 00:00:00');%year,month,day,hour,minute,second
currentDateTime = datetime(write_matrix(:, [3, 4, 5, 8, 9, 10]));
index = find(currentDateTime >= startdate & currentDateTime < enddate);
write_matrix = write_matrix(index,:);
Ml = write_matrix(:,6);

%% GR law calculation
magnitudeVector=write_matrix(:,6);
%divide in bins
mag_range=min(magnitudeVector):0.1:7;
[number,bin]= hist(magnitudeVector,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end

%% Mc MAXC (Max Curvature) method
mCatalog=write_matrix;
nMethod =1; % MAX CURVE OR nMethod=5; % BEST COMBINATION OR nMethod=6; % EMR
fBinning = 0.1; % size of our magnitude bins
fMcCorrection = 0; % Used to add a correction because maximum curvature often underestimates Mc
fMc = 0.8;


vSel = mCatalog(:,6) >= fMc;
mCat = mCatalog(vSel,:);
fBinning=0.1;
[fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
a1=fAvalue;b1=fBvalue;
disp(sprintf('Maximum Likelihood Estimation: b=%.2f+/-%.2f,a=%.2f %', b1,fStdDev,a1));

semilogy(mag_count(:,1),mag_count(:,3),'k^','Markersize',6) % cumulative events%
hold on;
Mc=fMc;
Max_GR=5.8;
x=fMc:0.1:Max_GR;
semilogy(x,power(10,a1-b1*x),'k','Linewidth',2); % GR law
grid on;
xlabel('Magnitude','fontsize',14);
ylabel('Frequency','fontsize',14);

%% before the MS -- MW
file_path = 'WP_50k_MS_20240807.csv';

disp(['CSV file name: ', file_path]);
dataTable = readtable(file_path);
write_matrix = table2array(dataTable);
startdate=datetime('2021-09-21 23:15:53');%year,month,day,hour,minute,second
enddate=datetime('2031-01-01 00:00:00');%year,month,day,hour,minute,second
currentDateTime = datetime(write_matrix(:, [3, 4, 5, 8, 9, 10]));
index = find(currentDateTime >= startdate & currentDateTime < enddate);
write_matrix = write_matrix(index,:);
%% Convert Ml to Mw
Ml = write_matrix(:,6);
Mw = zeros(2,1);
for i = 1:length(Ml)
    Mw(i) = Ml2Mw(Ml(i));
end
Mw = round(Mw,1);
write_matrix(:,6) = Mw;

%% GR law calculation
magnitudeVector=write_matrix(:,6);
%divide in bins
mag_range=min(magnitudeVector):0.1:7;
[number,bin]= hist(magnitudeVector,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end

mCatalog=write_matrix;
nMethod =1; % MAX CURVE OR nMethod=5; % BEST COMBINATION OR nMethod=6; % EMR
fBinning = 0.1; % size of our magnitude bins
fMcCorrection = 0; % Used to add a correction because maximum curvature often underestimates Mc
fMc = 1.68;
hold on;
vSel = mCatalog(:,6) >= fMc;
mCat = mCatalog(vSel,:);
fBinning=0.1;
[fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
a1=fAvalue;b1=fBvalue;
disp(sprintf('Maximum Likelihood Estimation: b=%.2f+/-%.2f,a=%.2f %', b1,fStdDev,a1));
semilogy(mag_count(:,1),mag_count(:,3),'k.','Markersize',15) % cumulative events%
hold on;
Mc=fMc;
Max_GR=5.8;
x=fMc:0.1:Max_GR;
semilogy(x,power(10,a1-b1*x),'k--','Linewidth',2); % GR law
grid on;
ylim([power(10,-1),power(10,4)]);
set(gca,'fontsize',12);
title('WPAS');