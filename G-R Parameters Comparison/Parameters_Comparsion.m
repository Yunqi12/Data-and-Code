%% this is the plot in the result's part
%% Ebel (2009)
clear all; clc; close all;
dataset = readtable('other_SCRs_study.csv');  % For numeric data
dataset = sortrows(dataset, dataset.Properties.VariableNames{1});  % Sort by the first column name

MS_Mag = table2array(dataset(:,1));
Mag_type = dataset(:,2);
Mc = table2array(dataset(:,3));
b_value = table2array(dataset(:,4));
p_value = table2array(dataset(:,5));
a_value = table2array(dataset(:,6));


% separate different magnitude type
index_Mw = find(strcmp(Mag_type.Mag_Type, 'MW'));
index_Ml = find(strcmp(Mag_type.Mag_Type, 'ML'));

% Plot the figure
figure;

% Use tiledlayout to control spacing
t = tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact'); 

% First subplot (b-value)
nexttile;
[x_Ml,y_Ml] = normal_distribution(b_value(index_Ml));
plot(x_Ml, y_Ml, '-', 'LineWidth', 2,'Color',[0.8500 0.3250 0.0980]);  
hold on;
[x_Mw,y_Mw] = normal_distribution(b_value(index_Mw));
plot(x_Mw, y_Mw, '-', 'LineWidth', 2,'color',[0 0.4470 0.7410]);  
y_limits = [min(ylim)-0.2,max(ylim)+0.2]; 
plot([0.76 0.76],y_limits,'k-','LineWidth',2); % WPAS ML
fill([0.76-0.02,0.76-0.02,0.76+0.02,0.76+0.02],[max(y_limits),min(y_limits),min(y_limits),max(y_limits)], ...
    [0.8500 0.3250 0.0980], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % WPAS ML (1std)
plot([1.07 1.07],y_limits,'k:','LineWidth',2); % WPAS MW
fill([1.07-0.03,1.07-0.03,1.07+0.03,1.07+0.03],[max(y_limits),min(y_limits),min(y_limits),max(y_limits)], ...
    [0 0.4470 0.7410], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % WPAS MW (1std)
plot([0.91 0.91],y_limits,'k--','LineWidth',2); % Generic Model
grid on;
legend({'M_L (Ebel)','M_W (Ebel)','WPAS M_L','WPAS M_L (1std)','WPAS M_W','WPAS M_W (1std)','Generic Model'},'FontSize',12);
ylim(y_limits);
xlabel('b value','FontSize',12);
ylabel('Probability Density','FontSize',12);
set(gca, 'FontSize', 14);

% Second subplot (p-value)
nexttile;
[x_Ml,y_Ml] = normal_distribution(p_value(index_Ml));
plot(x_Ml, y_Ml, '-', 'LineWidth', 2,'Color',[0.8500 0.3250 0.0980]);  
hold on;
[x_Mw,y_Mw] = normal_distribution(p_value(index_Mw));
plot(x_Mw, y_Mw, '-', 'LineWidth', 2,'color',[0 0.4470 0.7410]);  
[x,y] = normal_distribution(p_value);
plot(x, y, '-', 'LineWidth', 2,'color',[0.9290 0.6940 0.1250]);  
y_limits = [min(ylim)-0.2,max(ylim)+0.2]; 
plot([0.84 0.84],y_limits,'k-','LineWidth',2); % WPAS ML
plot([0.83 0.83],y_limits,'k:','LineWidth',2); % WPAS MW
plot([1.08 1.08],y_limits,'k--','LineWidth',2); % Generic Model
grid on;
legend({'M_L (Ebel)','M_W (Ebel)','M_L and M_W (Ebel)','WPAS M_L','WPAS M_W','Generic Model'},'FontSize',12);
ylim(y_limits);
xlabel('p value','FontSize',12);
ylabel('Probability Density','FontSize',12);
set(gca, 'FontSize', 14);

% Third subplot (a-value)
nexttile;
[x_Ml,y_Ml] = normal_distribution(a_value(index_Ml));
plot(x_Ml, y_Ml, '-', 'LineWidth', 2,'Color',[0.8500 0.3250 0.0980]);  
hold on;
[x_Mw,y_Mw] = normal_distribution(a_value(index_Mw));
plot(x_Mw, y_Mw, '-', 'LineWidth', 2,'color',[0 0.4470 0.7410]);  
y_limits = [min(ylim)-0.2,max(ylim)+0.2]; 
plot([-1.98 -1.98],y_limits,'k-','LineWidth',2); % WPAS ML
plot([-2.77 -2.77],y_limits,'k:','LineWidth',2); % WPAS MW
plot([-1.67 -1.67],y_limits,'k--','LineWidth',2); % Generic Model
grid on;
legend({'M_L (Ebel)','M_W (Ebel)','WPAS M_L','WPAS M_W','Generic Model'},'FontSize',12);
ylim(y_limits);
xlabel('a value','FontSize',12);
ylabel('Probability Density','FontSize',12);
set(gca, 'FontSize', 14);


function [x_values,y_values] = normal_distribution(data)
    mu = mean(data);        
    sigma = std(data);      
    % Generate x-values for the normal distribution curve
    x_values = linspace(min(data), max(data), 100);  
    y_values = normpdf(x_values, mu, sigma);  
end
