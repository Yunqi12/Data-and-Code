clear all; clc; close all;
%% the Mc was calculated in the Mc folder
%% before MS -- Mw catalogue
file_path = 'WP_50k_from2000_beforeMS.csv';
disp(['CSV file name: ', file_path]);
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
% Mw(1) = 5.9;
write_matrix(:,6) = Mw;
magnitude=write_matrix(:,6);
%divide in bins
mag_range=min(magnitude):0.1:7;
[number,bin]= hist(magnitude,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end

%% calculate b-value (b-value changed with the Mc value)
fMc = 0:0.1:4;
b_value_MLE = zeros(2,1);
mCatalog = write_matrix;
fBinning=0.1;
for i = 1:length(fMc)    
    vSel = []; mCat = [];
    vSel = mCatalog(:,6) >= fMc(i);
    mCat = mCatalog(vSel,:);    
    [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
    b_value_MLE(i,1) = fBvalue;
    b_value_MLE(i,2) = fStdDev;
end

%% figrue
figure;
magnitudeVector=Mw;
%divide in bins
mag_range=min(magnitudeVector):0.1:7;
[number,bin]= hist(magnitudeVector,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end

subplot(2,2,2);
errorbar(fMc,b_value_MLE(:,1),b_value_MLE(:,2),'k','LineWidth',1.5);
hold on;
ylimit = [0 4];
plot([1.8 1.8],ylimit,'k--','LineWidth',1.5);
plot([1.71 1.71],ylimit,'g', 'LineWidth', 1.5);
fill([1.71-0.18 1.71+0.18 1.71+0.18 1.71-0.18], [min(ylimit),min(ylimit),max(ylimit),max(ylimit)],'g','FaceAlpha',0.1,'EdgeColor', 'none');

plot([2.9 2.9],ylimit,'k-.','LineWidth',1.5);
plot([2.89 2.89],ylimit,'b', 'LineWidth', 1.5);
fill([2.89-0.18 2.89+0.18 2.89+0.18 2.89-0.18], [min(ylimit),min(ylimit),max(ylimit),max(ylimit)],'b','FaceAlpha',0.1,'EdgeColor', 'none');

plot([2.3 2.3],ylimit,'r','LineWidth',1.5);
fill([2.3-0.18 2.3+0.18 2.3+0.18 2.3-0.18], [min(ylimit),min(ylimit),max(ylimit),max(ylimit)],'r','FaceAlpha',0.1,'EdgeColor', 'none');

% legend({'b-value','MAXC','MAXC+Boostrapping (Ave)','MAXC+Boostrapping (Std)','MBS','MBS+Boostrapping (Ave)','MBS+Boostrapping (Std)','Preferred Mc','Preferred Mc (Range)'}, ...
%    'Location','northwest','fontsize',12);
xlabel('Magnitude','fontsize',14);
ylabel('b-value','fontsize',14);
title('Background Seismicity (M_W)','fontsize',14);
set(gca,'fontsize',12);
grid on;


%% Before MS ---ML
write_matrix(:,6) = Ml;
magnitude = [];mag_count = [];
magnitude=write_matrix(:,6);
%divide in bins
mag_range=min(magnitude):0.1:7;
[number,bin]= hist(magnitude,mag_range'); % calculate the number per bin

mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end
%% b-value calculation
fMc = 0:0.1:4;
b_value_MLE = zeros(2,1);
mCatalog = write_matrix;
fBinning=0.1;
for i = 1:length(fMc)    
    vSel = []; mCat = [];
    vSel = mCatalog(:,6) >= fMc(i);
    mCat = mCatalog(vSel,:);    
    [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
    b_value_MLE(i,1) = fBvalue;
    b_value_MLE(i,2) = fStdDev;
    vMag = mCat(:,6);
end

magnitudeVector=Ml;
%divide in bins
mag_range=min(magnitudeVector):0.1:7;
[number,bin]= hist(magnitudeVector,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end

subplot(2,2,1);
errorbar(fMc,b_value_MLE(:,1),b_value_MLE(:,2),'k','LineWidth',1.5);
hold on;
ylimit = [0 4];
plot([1 1],ylimit,'k--','LineWidth',1.5);
plot([0.9 0.9],ylimit,'g', 'LineWidth', 1.5);
fill([0.9-0.24 0.9+0.24 0.9+0.247 0.9-0.24], [min(ylimit),min(ylimit),max(ylimit),max(ylimit)],'g','FaceAlpha',0.1,'EdgeColor', 'none');

plot([3.3 3.3],ylimit,'k-.','LineWidth',1.5);
plot([3.16 3.16],ylimit,'b', 'LineWidth', 1.5);
fill([3.16-0.22 3.16+0.22 3.16+0.22 3.16-0.22], [min(ylimit),min(ylimit),max(ylimit),max(ylimit)],'b','FaceAlpha',0.1,'EdgeColor', 'none');

plot([2.03 2.03],ylimit,'r','LineWidth',1.5);
fill([2.03-0.23 2.03+0.23 2.03+0.23 2.03-0.23], [min(ylimit),min(ylimit),max(ylimit),max(ylimit)],'r','FaceAlpha',0.1,'EdgeColor', 'none');

xlabel('Magnitude','fontsize',14);
ylabel('b-value','fontsize',14);
title('Background Seismicity (M_L)','fontsize',14);
set(gca,'fontsize',12);
grid on;


%% after MS ---Mw
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
% Mw(1) = 5.9;
write_matrix(:,6) = Mw;

%% GR law calculation
magnitude=write_matrix(:,6);
%divide in bins
mag_range=min(magnitude):0.1:7;
[number,bin]= hist(magnitude,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end

fMc = 0:0.1:3;
b_value_MLE = zeros(2,1);
mCatalog = write_matrix;
fBinning=0.1;
for i = 1:length(fMc)
    vSel = []; mCat = [];
    vSel = mCatalog(:,6) >= fMc(i);
    mCat = mCatalog(vSel,:);    
    [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
    b_value_MLE(i,1) = fBvalue;
    b_value_MLE(i,2) = fStdDev;
end


subplot(2,2,4);
errorbar(fMc,b_value_MLE(:,1),b_value_MLE(:,2),'k','LineWidth',1.5);
hold on;
ylimit = [0 4];
plot([1.5 1.5],ylimit,'k--','LineWidth',1.5);
plot([1.5 1.5],ylimit,'g', 'LineWidth', 1.5);
fill([1.5-0.02 1.5+0.02 1.5+0.02 1.5-0.02], [min(ylimit),min(ylimit),max(ylimit),max(ylimit)],'g','FaceAlpha',0.1,'EdgeColor', 'none');

plot([1.9 1.9],ylimit,'k-.','LineWidth',1.5);
plot([1.85 1.85],ylimit,'b', 'LineWidth', 1.5);
fill([1.85-0.08 1.85+0.08 1.85+0.08 1.85-0.08], [min(ylimit),min(ylimit),max(ylimit),max(ylimit)],'b','FaceAlpha',0.1,'EdgeColor', 'none');

plot([1.68 1.68],ylimit,'r','LineWidth',1.5);
fill([1.68-0.06 1.68+0.06 1.68+0.06 1.68-0.06], [min(ylimit),min(ylimit),max(ylimit),max(ylimit)],'r','FaceAlpha',0.1,'EdgeColor', 'none');

xlabel('Magnitude','fontsize',14);
ylabel('b-value','fontsize',14);
title('WPAS (M_W)','fontsize',14);
set(gca,'fontsize',12);
grid on;


%% AfterMS ---Ml
write_matrix(:,6) = Ml;

%% GR law calculation
magnitude=write_matrix(:,6);
%divide in bins
mag_range=min(magnitude):0.1:7;
[number,bin]= hist(magnitude,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end

fMc = 0:0.1:3;
b_value_MLE = zeros(2,1);
mCatalog = write_matrix;
fBinning=0.1;
for i = 1:length(fMc)    
    vSel = []; mCat = [];
    vSel = mCatalog(:,6) >= fMc(i);
    mCat = mCatalog(vSel,:);    
    [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
    b_value_MLE(i,1) = fBvalue;
    b_value_MLE(i,2) = fStdDev;
    vMag = mCat(:,6);
end

magnitudeVector=Ml;
%divide in bins
mag_range=min(magnitudeVector):0.1:7;
[number,bin]= hist(magnitudeVector,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end

subplot(2,2,3);
errorbar(fMc,b_value_MLE(:,1),b_value_MLE(:,2),'k','LineWidth',1.5);
hold on;
ylimit = [0 4];
plot([0.6 0.6],ylimit,'k--','LineWidth',1.5);
plot([0.61 0.61],ylimit,'g', 'LineWidth', 1.5);
fill([0.61-0.06 0.61+0.06 0.61+0.06 0.61-0.06], [min(ylimit),min(ylimit),max(ylimit),max(ylimit)],'g','FaceAlpha',0.1,'EdgeColor', 'none');

plot([0.9 0.9],ylimit,'k-.','LineWidth',1.5);
plot([1 1],ylimit,'b', 'LineWidth', 1.5);
fill([1-0.13 1+0.13 1+0.13 1-0.13], [min(ylimit),min(ylimit),max(ylimit),max(ylimit)],'b','FaceAlpha',0.1,'EdgeColor', 'none');

plot([0.8 0.8],ylimit,'r','LineWidth',1.5);
fill([0.8-0.1 0.8+0.1 0.8+0.1 0.8-0.1], [min(ylimit),min(ylimit),max(ylimit),max(ylimit)],'r','FaceAlpha',0.1,'EdgeColor', 'none');

legend({'b-value','MAXC','MAXC+Boostrapping (Ave)','MAXC+Boostrapping (Std)','MBS','MBS+Boostrapping (Ave)','MBS+Boostrapping (Std)','Preferred Mc','Preferred Mc (Range)'}, ...
    'Location','northwest','FontSize',12);


xlabel('Magnitude','fontsize',14);
ylabel('b-value','fontsize',14);
title('WPAS (M_L)','fontsize',14);
set(gca,'fontsize',12);
grid on;

