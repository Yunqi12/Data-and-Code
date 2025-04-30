% calculate the seismic moment (energy) release per month before the mainshock
clc;
clear all;
file_path = 'WP_50k_from2000_beforeMS.csv';
disp(['CSV file name: ', file_path]);
dataTable = readtable(file_path);
write_matrix = table2array(dataTable);
currentDateTime = datetime(write_matrix(:, [3, 4, 5, 8, 9, 10]));
mag = write_matrix(:, 6);
Energy = zeros(2,1);

for i = 1:length(mag)
    Mw(i) = Ml2Mw(mag(i));
end
Energy = power(10,3/2.*(Mw+10.7)); % emprical equation for Mo conversion


index = find(currentDateTime >= datetime(2000,01,01));
write_matrix = write_matrix(index,:);
currentDateTime = datetime(write_matrix(:, [3, 4, 5, 8, 9, 10]));


%% calculate month seismic moment released rate
startdate = datetime(2000,1,1);
enddate = startdate+calmonths(1);
cum_event = 0;
count1 = 1;
while startdate <= currentDateTime(end)
    index = find(currentDateTime >= startdate & currentDateTime < enddate);
    Event_month(count1,1) = year(startdate); Event_month(count1,2) = month(startdate);
    Event_month(count1,3) = sum(Energy(index));
    startdate = enddate;
    enddate = startdate+calmonths(1);
    count1 = count1 +1;
end

save('Energy_month_beforeMS.mat','Event_month');
