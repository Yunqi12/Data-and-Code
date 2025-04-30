%calculate the Mc value based on 4 methods
clc;close all;clear all;
%% Pre-mainshock ML catalogue
file_path = 'WP_50k_from2000_beforeMS.csv';
disp('====================================================');
disp(['Processing File (pre-MS, ML): ', file_path]);
disp('====================================================');

dataTable = readtable(file_path);
write_matrix = table2array(dataTable);
startdate=datetime('2000-01-01 00:00:00');%year,month,day,hour,minute,second
enddate=datetime('2031-01-01 00:00:00');%year,month,day,hour,minute,second
currentDateTime = datetime(write_matrix(:, [3, 4, 5, 8, 9, 10]));
index = find(currentDateTime >= startdate & currentDateTime < enddate);
write_matrix = write_matrix(index,:);
Ml = write_matrix(:,6);
magnitudeVector=write_matrix(:,6);
%divide in bins
mag_range=min(magnitudeVector):0.1:7;
[number,bin]= hist(magnitudeVector,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end

%% MaxC method -- used as the minimum magnitude cutoff (underestimate the Mc)
mCatalog=write_matrix;
nMethod =1; % MAX CURVE OR nMethod=5; % BEST COMBINATION OR nMethod=6; % EMR
fBinning = 0.1; % size of our magnitude bins
fMcCorrection = 0; % Used to add a correction because maximum curvature often underestimates Mc
[fMc_maxc] = calc_Mc(mCatalog, nMethod, fBinning, fMcCorrection);
disp(['<strong>MAXC: </strong>' num2str(fMc_maxc)]);
%% Boostrapping+MAXC
magnitudeVector_sort = sort(magnitudeVector);
out=bootrsp(magnitudeVector_sort,10000);
boostrap_Maxc_Mc = zeros(10000,1)+1;
for k = 1:1:10000
    mCatalog = zeros(length(out(:,k)),6);
    mCatalog(:,6) = out(:,k);
    [fMc_maxc] = calc_Mc(mCatalog, nMethod, fBinning, fMcCorrection);
    boostrap_Maxc_Mc(k) = fMc_maxc;
end
index = ~isnan(boostrap_Maxc_Mc);
valid_boostrap_Mc = find(index == 1);
fMc_boot_Mc = mean(boostrap_Maxc_Mc(valid_boostrap_Mc));
fMc_boot_Mc_std = std(boostrap_Maxc_Mc(valid_boostrap_Mc));
disp(['<strong>MAXC + Boostrapping: </strong>' num2str(fMc_boot_Mc) ' ± ' num2str(fMc_boot_Mc_std)]);
disp(['The number of Valid Bootstrapping: ', num2str(length(valid_boostrap_Mc))]);

%% Boostrapping+MBS
boostrap_Mc = zeros(10000,1)+1;
for k = 1:1:10000
    minM = fMc_maxc;   
    fBinning=0.1;
    b_MLE = zeros(2,1);
    b_unc = zeros(2,1);
    b_ave = zeros(2,1);
    num = 1;
    for i = minM:0.1:4
        fMc = i;
        vSel = out(:,k) >= fMc;
        mCat = out(vSel,k);    
        [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
        b_MLE(num) = fBvalue;
        b_unc(num) = fStdDev;
        b_ave_value = 0;
        b_range = 0.5;
        for j = fMc:0.1:fMc+b_range
            vSel = out(:,k) >= j;
            mCat = out(vSel,k);    
            [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
            b_ave_value = b_ave_value + fBvalue;
        end
        b_ave_global = 1;
        num = num + 1;
    end
    index = find(abs(b_ave_global-b_MLE) <= b_unc);
    mag_Mc = minM:0.1:4;
    if isempty(index)
        boostrap_Mc(k) = nan;
    else
        boostrap_Mc(k) = mag_Mc(index(1));
    end
end

index = ~isnan(boostrap_Mc);
valid_boostrap_Mc = find(index == 1);
fMc_boot = mean(boostrap_Mc(valid_boostrap_Mc));
fMc_boot_std = std(boostrap_Mc(valid_boostrap_Mc));
disp(['<strong>MBS + Boostrapping: </strong>' num2str(fMc_boot) ' ± ' num2str(fMc_boot_std)]);
disp(['The number of Valid Bootstrapping: ', num2str(length(valid_boostrap_Mc))]);
%% MBS method (Mc by b-value stability)
minM = fMc_maxc;
% minM = 0;
fBinning=0.1;
b_MLE = zeros(2,1);
b_unc = zeros(2,1);
b_ave = zeros(2,1);
num = 1;
for i = minM:0.1:4
    fMc = i;
    vSel = write_matrix(:,6) >= fMc;
    mCat = write_matrix(vSel,:);    
    [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
    b_MLE(num) = fBvalue;
    b_unc(num) = fStdDev;
    b_ave_value = 0;
    b_range = 0.5;
    for j = fMc:0.1:fMc+b_range
        vSel = write_matrix(:,6) >= j;
        mCat = write_matrix(vSel,:);    
        [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
        b_ave_value = b_ave_value + fBvalue;
    end
    b_ave(num) = b_ave_value*fBinning/0.6;
    num = num + 1;
end
b_ave_global = ones(length(b_ave),1);
index = find(abs(b_ave_global-b_MLE) <= b_unc);
mag_Mc = minM:0.1:4;
fMc_MBS = mag_Mc(index(1));
disp(['<strong>MBS: </strong>' num2str(fMc_MBS)]);

%% figure
figure;
subplot(2,2,1);
semilogy(mag_count(:,1),mag_count(:,2),'k*','Markersize',8); % events per bin
hold on;
semilogy(mag_count(:,1),mag_count(:,3),'r.','Markersize',15) % cumulative events%
hold on;
y_limits = ylim;

plot([fMc_maxc fMc_maxc], y_limits, 'LineWidth', 1.2, 'Color', '#EDB120');
plot([fMc_boot_Mc fMc_boot_Mc], y_limits, 'g--', 'LineWidth', 1.2); 
fill([fMc_boot_Mc-fMc_boot_Mc_std fMc_boot_Mc+fMc_boot_Mc_std fMc_boot_Mc+fMc_boot_Mc_std fMc_boot_Mc-fMc_boot_Mc_std],[y_limits(1) y_limits(1) y_limits(2) y_limits(2)],'g', 'FaceAlpha', 0.1,'EdgeColor','none');

plot([fMc_MBS fMc_MBS], y_limits, 'LineWidth', 1.2, 'Color', '#7E2F8E');
plot([fMc_boot fMc_boot], y_limits, 'b--', 'LineWidth', 1.2); 

fill([fMc_boot-fMc_boot_std fMc_boot+fMc_boot_std fMc_boot+fMc_boot_std fMc_boot-fMc_boot_std],[y_limits(1) y_limits(1) y_limits(2) y_limits(2)],'b', 'FaceAlpha', 0.1,'EdgeColor','none');

grid on; legend('Non Cum. FMD','Cum. FMD','MAXC','MAXC+Boostrapping (Ave.)','MAXC+Boostrapping (Std)','MBS','MBS+Boostrapping (Ave.)','MBS+Boostrapping (Std)');
xlabel('Magnitude','FontSize', 14);
ylabel('Frequency','FontSize', 14);
title('Background Seismicity (Ml)','FontSize', 14);

%% Pre-mainshock MW catalogue
file_path = 'WP_50k_from2000_beforeMS.csv';
disp('====================================================');
disp(['Processing File (pre-MS, MW): ', file_path]);
disp('====================================================');

% read the dataset
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

%% count number of events per bin
magnitudeVector=write_matrix(:,6);
%divide in bins
mag_range=min(magnitudeVector):0.1:7;
[number,bin]= hist(magnitudeVector,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end


%% MaxC method -- used as the minimum magnitude cutoff (underestimate the Mc)
mCatalog=write_matrix;
nMethod =1; % MAX CURVE OR nMethod=5; % BEST COMBINATION OR nMethod=6; % EMR
fBinning = 0.1; % size of our magnitude bins
fMcCorrection = 0; % Used to add a correction because maximum curvature often underestimates Mc
[fMc_maxc] = calc_Mc(mCatalog, nMethod, fBinning, fMcCorrection);
disp(['<strong>MAXC: </strong>' num2str(fMc_maxc)]);
%% Boostrapping+MAXC
magnitudeVector_sort = sort(magnitudeVector);

out=bootrsp(magnitudeVector_sort,10000);
boostrap_Maxc_Mc = zeros(10000,1)+1;
for k = 1:1:10000
    mCatalog = zeros(length(out(:,k)),6);
    mCatalog(:,6) = out(:,k);
    [fMc_maxc] = calc_Mc(mCatalog, nMethod, fBinning, fMcCorrection);
    boostrap_Maxc_Mc(k) = fMc_maxc;
end


index = ~isnan(boostrap_Maxc_Mc);
valid_boostrap_Mc = find(index == 1); 
fMc_boot_Mc = mean(boostrap_Maxc_Mc(valid_boostrap_Mc));
fMc_boot_Mc_std = std(boostrap_Maxc_Mc(valid_boostrap_Mc));
disp(['<strong>MAXC + Boostrapping: </strong>' num2str(fMc_boot_Mc) ' ± ' num2str(fMc_boot_Mc_std)]);
disp(['The number of Valid Bootstrapping: ', num2str(length(valid_boostrap_Mc))]);
%% Boostrapping+MBS
boostrap_Mc = zeros(10000,1)+1;
for k = 1:1:10000
    minM = fMc_maxc;   
    fBinning=0.1;
    b_MLE = zeros(2,1);
    b_unc = zeros(2,1);
    b_ave = zeros(2,1);
    num = 1;
    for i = minM:0.1:4
        fMc = i;
        vSel = out(:,k) >= fMc;
        mCat = out(vSel,k);    
        [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
        b_MLE(num) = fBvalue;
        b_unc(num) = fStdDev;
        b_ave_value = 0;
        b_range = 0.5;
        for j = fMc:0.1:fMc+b_range
            vSel = out(:,k) >= j;
            mCat = out(vSel,k);    
            [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
            b_ave_value = b_ave_value + fBvalue;
        end
        b_ave_global = 1;
        num = num + 1;
    end
    index = find(abs(b_ave_global-b_MLE) <= b_unc);
    mag_Mc = minM:0.1:4;
    if isempty(index)
        boostrap_Mc(k) = nan;
    else
        boostrap_Mc(k) = mag_Mc(index(1));
    end
end


index = ~isnan(boostrap_Mc);
valid_boostrap_Mc = find(index == 1);
fMc_boot = mean(boostrap_Mc(valid_boostrap_Mc));
fMc_boot_std = std(boostrap_Mc(valid_boostrap_Mc));
disp(['<strong>MBS + Boostrapping: </strong>' num2str(fMc_boot) ' ± ' num2str(fMc_boot_std)]);
disp(['The number of Valid Bootstrapping: ', num2str(length(valid_boostrap_Mc))]);

%% MBS method (Mc by b-value stability)
minM = fMc_maxc;
fBinning=0.1;
b_MLE = zeros(2,1);
b_unc = zeros(2,1);
b_ave = zeros(2,1);
num = 1;
for i = minM:0.1:4
    fMc = i;
    vSel = write_matrix(:,6) >= fMc;
    mCat = write_matrix(vSel,:);    
    [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
    b_MLE(num) = fBvalue;
    b_unc(num) = fStdDev;
    b_ave_value = 0;
    b_range = 0.5;
    for j = fMc:0.1:fMc+b_range
        vSel = write_matrix(:,6) >= j;
        mCat = write_matrix(vSel,:);    
        [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
        b_ave_value = b_ave_value + fBvalue;
    end
    b_ave(num) = b_ave_value*fBinning/0.6;
    num = num + 1;
end
b_ave_global = ones(length(b_ave),1);
index = find(abs(b_ave_global-b_MLE) <= b_unc);
mag_Mc = minM:0.1:4;
fMc_MBS = mag_Mc(index(1));

disp(['<strong>MBS: </strong>' num2str(fMc_MBS)]);

%% figure
subplot(2,2,2);
semilogy(mag_count(:,1),mag_count(:,2),'k*','Markersize',8); % events per bin]
hold on;
semilogy(mag_count(:,1),mag_count(:,3),'r.','Markersize',15) % cumulative events%
hold on;
y_limits = ylim;

plot([fMc_maxc fMc_maxc], y_limits, 'LineWidth', 1.2, 'Color', '#EDB120');
plot([fMc_boot_Mc fMc_boot_Mc], y_limits, 'g--', 'LineWidth', 1.2);  
fill([fMc_boot_Mc-fMc_boot_Mc_std fMc_boot_Mc+fMc_boot_Mc_std fMc_boot_Mc+fMc_boot_Mc_std fMc_boot_Mc-fMc_boot_Mc_std],[y_limits(1) y_limits(1) y_limits(2) y_limits(2)],'g', 'FaceAlpha', 0.1,'EdgeColor','none');

plot([fMc_MBS fMc_MBS], y_limits, 'LineWidth', 1.2, 'Color', '#7E2F8E');
plot([fMc_boot fMc_boot], y_limits, 'b--', 'LineWidth', 1.2);

fill([fMc_boot-fMc_boot_std fMc_boot+fMc_boot_std fMc_boot+fMc_boot_std fMc_boot-fMc_boot_std],[y_limits(1) y_limits(1) y_limits(2) y_limits(2)],'b', 'FaceAlpha', 0.1,'EdgeColor','none');

grid on; legend('Non Cum. FMD','Cum. FMD','MAXC','MAXC+Boostrapping (Ave.)','MAXC+Boostrapping (Std)','MBS','MBS+Boostrapping (Ave.)','MBS+Boostrapping (Std)');
xlabel('Magnitude','FontSize', 14);
ylabel('Frequency','FontSize', 14);
title('Background Seismicity (Mw)','FontSize', 14);

%% WPAS ML catalogue
file_path = 'WP_50k_MS_20240807.csv';
disp('====================================================');
disp(['Processing File (Post-MS, ML): ', file_path]);
disp('====================================================');

dataTable = readtable(file_path);
write_matrix = table2array(dataTable);

startdate= datetime('2021-09-21 23:15:53');%year,month,day,hour,minute,second
enddate=datetime('2031-01-01 00:00:00');%year,month,day,hour,minute,second
currentDateTime = datetime(write_matrix(:, [3, 4, 5, 8, 9, 10]));
index = find(currentDateTime >= startdate & currentDateTime < enddate);

write_matrix = write_matrix(index,:);
Ml = write_matrix(:,6);
magnitudeVector=write_matrix(:,6);
%divide in bins
mag_range=min(magnitudeVector):0.1:7;
[number,bin]= hist(magnitudeVector,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end


%% MaxC method -- used as the minimum magnitude cutoff (underestimate the Mc)
mCatalog=write_matrix;
nMethod =1; % MAX CURVE OR nMethod=5; % BEST COMBINATION OR nMethod=6; % EMR
fBinning = 0.1; % size of our magnitude bins
fMcCorrection = 0; % Used to add a correction because maximum curvature often underestimates Mc
[fMc_maxc] = calc_Mc(mCatalog, nMethod, fBinning, fMcCorrection);
disp(['<strong>MAXC: </strong>' num2str(fMc_maxc)]);
%% Boostrapping+MAXC
magnitudeVector_sort = sort(magnitudeVector);

out=bootrsp(magnitudeVector_sort,10000);
boostrap_Maxc_Mc = zeros(10000,1)+1;
for k = 1:1:10000
    mCatalog = zeros(length(out(:,k)),6);
    mCatalog(:,6) = out(:,k);
    [fMc_maxc] = calc_Mc(mCatalog, nMethod, fBinning, fMcCorrection);
    boostrap_Maxc_Mc(k) = fMc_maxc;
end

index = ~isnan(boostrap_Maxc_Mc);
valid_boostrap_Mc = find(index == 1);
fMc_boot_Mc = mean(boostrap_Maxc_Mc(valid_boostrap_Mc));
fMc_boot_Mc_std = std(boostrap_Maxc_Mc(valid_boostrap_Mc));
disp(['<strong>MAXC + Boostrapping: </strong>' num2str(fMc_boot_Mc) ' ± ' num2str(fMc_boot_Mc_std)]);
disp(['The number of Valid Bootstrapping: ', num2str(length(valid_boostrap_Mc))]);
%% Boostrapping+MBS
boostrap_Mc = zeros(10000,1)+1;


for k = 1:1:10000
    minM = fMc_maxc;   
    fBinning=0.1;
    b_MLE = zeros(2,1);
    b_unc = zeros(2,1);
    b_ave = zeros(2,1);
    num = 1;
    for i = minM:0.1:4
        fMc = i;
        vSel = out(:,k) >= fMc;
        mCat = out(vSel,k);    
        [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
        b_MLE(num) = fBvalue;
        b_unc(num) = fStdDev;
        b_ave_value = 0;
        b_range = 0.5;
        for j = fMc:0.1:fMc+b_range
            vSel = out(:,k) >= j;
            mCat = out(vSel,k);    
            [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
            b_ave_value = b_ave_value + fBvalue;
        end
        b_ave(num) = b_ave_value*fBinning/0.6;
        %b_ave_global = 1;
        num = num + 1;
    end
    %index = find(abs(b_ave_global-b_MLE) <= b_unc);
    index = find(abs(b_ave-b_MLE) <= b_unc);
    mag_Mc = minM:0.1:4;
    if isempty(index)
        boostrap_Mc(k) = nan;
    else
        boostrap_Mc(k) = mag_Mc(index(1));
    end
end

index = ~isnan(boostrap_Mc);
valid_boostrap_Mc = find(index == 1); 
fMc_boot = mean(boostrap_Mc(valid_boostrap_Mc));
fMc_boot_std = std(boostrap_Mc(valid_boostrap_Mc));
disp(['<strong>MBS + Boostrapping: </strong>' num2str(fMc_boot) ' ± ' num2str(fMc_boot_std)]);
disp(['The number of Valid Bootstrapping: ', num2str(length(valid_boostrap_Mc))]);


%% MBS method (Mc by b-value stability)
minM = fMc_maxc;
fBinning=0.1;
b_MLE = zeros(2,1);
b_unc = zeros(2,1);
b_ave = zeros(2,1);
num = 1;
for i = minM:0.1:4
    fMc = i;
    vSel = write_matrix(:,6) >= fMc;
    mCat = write_matrix(vSel,:);    
    [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
    b_MLE(num) = fBvalue;
    b_unc(num) = fStdDev;
    b_ave_value = 0;
    b_range = 0.5;
    for j = fMc:0.1:fMc+b_range
        vSel = write_matrix(:,6) >= j;
        mCat = write_matrix(vSel,:);    
        [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
        b_ave_value = b_ave_value + fBvalue;
    end
    b_ave(num) = b_ave_value*fBinning/0.6;
    num = num + 1;
end
b_ave_global = ones(length(b_ave),1);
% index = find(abs(b_ave_global-b_MLE) <= b_unc);
index = find(abs(b_ave-b_MLE) <= b_unc);
mag_Mc = minM:0.1:4;
fMc_MBS = mag_Mc(index(1));

disp(['<strong>MBS: </strong>' num2str(fMc_MBS)]);

%% figure
subplot(2,2,3);
semilogy(mag_count(:,1),mag_count(:,2),'k*','Markersize',8); % events per bin]
hold on;
semilogy(mag_count(:,1),mag_count(:,3),'r.','Markersize',15) % cumulative events%
hold on;
ylim([power(10,0),power(10,4)]);
y_limits = ylim;


plot([fMc_maxc fMc_maxc], y_limits, 'LineWidth', 1.2, 'Color', '#EDB120');
plot([fMc_boot_Mc fMc_boot_Mc], y_limits, 'g--', 'LineWidth', 1.2);  
fill([fMc_boot_Mc-fMc_boot_Mc_std fMc_boot_Mc+fMc_boot_Mc_std fMc_boot_Mc+fMc_boot_Mc_std fMc_boot_Mc-fMc_boot_Mc_std],[y_limits(1) y_limits(1) y_limits(2) y_limits(2)],'g', 'FaceAlpha', 0.1,'EdgeColor','none');

plot([fMc_MBS fMc_MBS], y_limits, 'LineWidth', 1.2, 'Color', '#7E2F8E');
plot([fMc_boot fMc_boot], y_limits, 'b--', 'LineWidth', 1.2); 

fill([fMc_boot-fMc_boot_std fMc_boot+fMc_boot_std fMc_boot+fMc_boot_std fMc_boot-fMc_boot_std],[y_limits(1) y_limits(1) y_limits(2) y_limits(2)],'b', 'FaceAlpha', 0.1,'EdgeColor','none');
grid on; legend('Non Cum. FMD','Cum. FMD','MAXC','MAXC+Boostrapping (Ave.)','MAXC+Boostrapping (Std)','MBS','MBS+Boostrapping (Ave.)','MBS+Boostrapping (Std)');
xlabel('Magnitude','FontSize', 14);
ylabel('Frequency','FontSize', 14);
title('WPAS (Ml)','FontSize', 14);

%% WPAS MW catalogue
file_path = 'WP_50k_MS_20240807.csv';
disp('====================================================');
disp(['Processing File (Post-MS, MW): ', file_path]);
disp('====================================================');

dataTable = readtable(file_path);
write_matrix = table2array(dataTable);

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

%% Count events per bin
magnitudeVector=write_matrix(:,6);
%divide in bins
mag_range=min(magnitudeVector):0.1:7;
[number,bin]= hist(magnitudeVector,mag_range'); % calculate the number per bin
mag_count=[bin,number'];
for i=1:length(mag_count)
    mag_count(i,3)=sum(mag_count(i:end,2));
end


%% MaxC method -- used as the minimum magnitude cutoff (underestimate the Mc)
mCatalog=write_matrix;
nMethod =1; % MAX CURVE OR nMethod=5; % BEST COMBINATION OR nMethod=6; % EMR
fBinning = 0.1; % size of our magnitude bins
fMcCorrection = 0; % Used to add a correction because maximum curvature often underestimates Mc
[fMc_maxc] = calc_Mc(mCatalog, nMethod, fBinning, fMcCorrection);
disp(['<strong>MAXC: </strong>' num2str(fMc_maxc)]);
%% Boostrapping+MAXC
magnitudeVector_sort = sort(magnitudeVector);

out=bootrsp(magnitudeVector_sort,10000);
boostrap_Maxc_Mc = zeros(10000,1)+1;
for k = 1:1:10000
    mCatalog = zeros(length(out(:,k)),6);
    mCatalog(:,6) = out(:,k);
    [fMc_maxc] = calc_Mc(mCatalog, nMethod, fBinning, fMcCorrection);
    boostrap_Maxc_Mc(k) = fMc_maxc;
end

index = ~isnan(boostrap_Maxc_Mc);
valid_boostrap_Mc = find(index == 1);
fMc_boot_Mc = mean(boostrap_Maxc_Mc(valid_boostrap_Mc));
fMc_boot_Mc_std = std(boostrap_Maxc_Mc(valid_boostrap_Mc));
disp(['<strong>MAXC + Boostrapping: </strong>' num2str(fMc_boot_Mc) ' ± ' num2str(fMc_boot_Mc_std)]);
disp(['The number of Valid Bootstrapping: ', num2str(length(valid_boostrap_Mc))]);
%% Boostrapping+MBS
boostrap_Mc = zeros(10000,1)+1;
for k = 1:1:10000
    minM = fMc_maxc;   
    fBinning=0.1;
    b_MLE = zeros(2,1);
    b_unc = zeros(2,1);
    b_ave = zeros(2,1);
    num = 1;
    for i = minM:0.1:4
        fMc = i;
        vSel = out(:,k) >= fMc;
        mCat = out(vSel,k);    
        [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
        b_MLE(num) = fBvalue;
        b_unc(num) = fStdDev;
        b_ave_value = 0;
        b_range = 0.5;
        for j = fMc:0.1:fMc+b_range
            vSel = out(:,k) >= j;
            mCat = out(vSel,k);    
            [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
            b_ave_value = b_ave_value + fBvalue;
        end
        b_ave(num) = b_ave_value*fBinning/0.6;
        % b_ave_global = 1;
        num = num + 1;
    end
    % index = find(abs(b_ave_global-b_MLE) <= b_unc);
    index = find(abs(b_ave-b_MLE) <= b_unc);
    mag_Mc = minM:0.1:4;
    if isempty(index)
        boostrap_Mc(k) = nan;
    else
        boostrap_Mc(k) = mag_Mc(index(1));
    end
end

index = ~isnan(boostrap_Mc);
valid_boostrap_Mc = find(index == 1);
fMc_boot = mean(boostrap_Mc(valid_boostrap_Mc));
fMc_boot_std = std(boostrap_Mc(valid_boostrap_Mc));
disp(['<strong>MBS + Boostrapping: </strong>' num2str(fMc_boot) ' ± ' num2str(fMc_boot_std)]);
disp(['The number of Valid Bootstrapping: ', num2str(length(valid_boostrap_Mc))]);

%% MBS method (Mc by b-value stability)
minM = fMc_maxc;
% minM = 0;
fBinning=0.1;
b_MLE = zeros(2,1);
b_unc = zeros(2,1);
b_ave = zeros(2,1);
num = 1;
for i = minM:0.1:4
    fMc = i;
    vSel = write_matrix(:,6) >= fMc;
    mCat = write_matrix(vSel,:);    
    [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
    b_MLE(num) = fBvalue;
    b_unc(num) = fStdDev;
    b_ave_value = 0;
    b_range = 0.5;
    for j = fMc:0.1:fMc+b_range
        vSel = write_matrix(:,6) >= j;
        mCat = write_matrix(vSel,:);    
        [fMeanMag, fBvalue, fStdDev, fAvalue] =  calc_bmemag(mCat, fBinning);
        b_ave_value = b_ave_value + fBvalue;
    end
    b_ave(num) = b_ave_value*fBinning/0.6;
    num = num + 1;
end
b_ave_global = ones(length(b_ave),1);
% index = find(abs(b_ave_global-b_MLE) <= b_unc);
index = find(abs(b_ave-b_MLE) <= b_unc);
mag_Mc = minM:0.1:4;
fMc_MBS = mag_Mc(index(1));


disp(['<strong>MBS: </strong>' num2str(fMc_MBS)]);
%% Figure
subplot(2,2,4);
semilogy(mag_count(:,1),mag_count(:,2),'k*','Markersize',8); % events per bin]
hold on;
semilogy(mag_count(:,1),mag_count(:,3),'r.','Markersize',15) % cumulative events%
hold on;
y_limits = ylim;

plot([fMc_maxc fMc_maxc], y_limits, 'LineWidth', 1.2, 'Color', '#EDB120');
plot([fMc_boot_Mc fMc_boot_Mc], y_limits, 'g--', 'LineWidth', 1.2); 
fill([fMc_boot_Mc-fMc_boot_Mc_std fMc_boot_Mc+fMc_boot_Mc_std fMc_boot_Mc+fMc_boot_Mc_std fMc_boot_Mc-fMc_boot_Mc_std],[y_limits(1) y_limits(1) y_limits(2) y_limits(2)],'g', 'FaceAlpha', 0.1,'EdgeColor','none');

plot([fMc_MBS fMc_MBS], y_limits, 'LineWidth', 1.2, 'Color', '#7E2F8E');
plot([fMc_boot fMc_boot], y_limits, 'b--', 'LineWidth', 1.2); 

fill([fMc_boot-fMc_boot_std fMc_boot+fMc_boot_std fMc_boot+fMc_boot_std fMc_boot-fMc_boot_std],[y_limits(1) y_limits(1) y_limits(2) y_limits(2)],'b', 'FaceAlpha', 0.1,'EdgeColor','none');

grid on; legend('Non Cum. FMD','Cum. FMD','MAXC','MAXC+Boostrapping (Ave.)','MAXC+Boostrapping (Std)','MBS','MBS+Boostrapping (Ave.)','MBS+Boostrapping (Std)');
xlabel('Magnitude','FontSize', 14);
ylabel('Frequency','FontSize', 14);
title('WPAS (Mw)','FontSize', 14);


