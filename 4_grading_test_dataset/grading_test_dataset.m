clc
close all
clear all

%% read in data

training_data = readmatrix('SOFT_STIFF_DATA.xlsx');

soft = training_data(:,15) < 1;
stiff = training_data(:,15) > 0;

test_data = readmatrix('AZA_DATA');

%% grading test dataset with logistic regression models

manual_grade = test_data(:,13);

cell_area = test_data(:,1);

DAPI_area = test_data(:,7);

cell_ar = test_data(:,5);

nuc_ar = test_data(:,11);

DOA = test_data(:,14);

% threshold values and coefficients determined from
% 'logistic_regression_model_k_fold_cross_validation'
grade_cell_area = cell_area > 1516.30;

grade_DAPI_area = DAPI_area > 140.96;

grade_cell_ar = cell_ar > 6.17;

grade_nuc_ar = nuc_ar > 1.86;

grade_DOA = DOA > 0.72;

temp_MLR = 1./(1+exp(-(-9.24+4.04e-4*cell_area+0.01*DAPI_area+-0.26*cell_ar+0.54*nuc_ar+10.21*DOA)));
grade_MLR = temp_MLR > 0.5;

%% percent corrrect for test dataset

% percent_correct = ((total # of cells) - (# of cells graded incorrectly)) / (total # of
% cells) * 100
cell_area_percent_correct = (length(test_data)-sum(abs(manual_grade - grade_cell_area)))/length(test_data)*100;
nuc_area_percent_correct = (length(test_data)-sum(abs(manual_grade - grade_DAPI_area)))/length(test_data)*100;
cell_ar_percent_correct = (length(test_data)-sum(abs(manual_grade - grade_cell_ar)))/length(test_data)*100;
nuc_ar_percent_correct = (length(test_data)-sum(abs(manual_grade - grade_nuc_ar)))/length(test_data)*100;
DOA_percent_correct = (length(test_data)-sum(abs(manual_grade - grade_DOA)))/length(test_data)*100;
MLR_percent_correct = (length(test_data)-sum(abs(manual_grade - grade_MLR)))/length(test_data)*100;

data_to_plot = [cell_area_percent_correct,nuc_area_percent_correct,cell_ar_percent_correct,nuc_ar_percent_correct,...
    DOA_percent_correct,MLR_percent_correct];

figure
h = bar(data_to_plot); hold on;
h.FaceColor = [0.6 0.6 0.6];
h.EdgeColor = 'k';
xlim([0 7])
ylim([40 100])
ylabel('% correct','FontWeight','bold')
title("performance on test dataset")
xticklabels({'cell area','nuc area','cell AR','nuc AR','DOA','MLR'})
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(get(gca, 'Title'), 'FontWeight', 'bold');
set(gca,'FontSize',25);
set(h,'LineWidth', 4);
pbaspect([2 1 1])
% exportgraphics(gcf,'test_percent_correct.png','Resolution',300)

%% percent quiscient/activated in soft, stiff, and AZA condition compared between manual grading and MLR

% grading training dataset with logistic regression models

training_data_manual_grade = training_data(:,13);

% threshold values and coefficients determined from
% 'logistic_regression_model_k_fold_cross_validation'
training_data_cell_area = training_data(:,1);
training_data_DAPI_area = training_data(:,7);
training_data_cell_ar = training_data(:,5);
training_data_nuc_ar = training_data(:,11);
training_data_DOA = training_data(:,14);

training_cell_area_grade = training_data_cell_area > 1516.30;
training_nuc_area_grade = training_data_DAPI_area > 140.96;
training_cell_ar_grade = training_data_cell_ar > 6.17;
training_nuc_ar_grade = training_data_nuc_ar > 1.86;
training_DOA_grade = training_data_DOA > 0.72;
training_data_MLR = 1./(1+exp(-(-9.24+4.04e-4*training_data_cell_area+0.01*training_data_DAPI_area+-0.26*training_data_cell_ar+0.54*training_data_nuc_ar+10.21*training_data_DOA)));
training_MLR_grade = training_data_MLR > 0.50;

% Computations: 
% percent_q = ((total # of cells in condition) - (total # of activated
% cells in condition)) / (total # of cells in condition) * 100

% percent_a = (total # of activated cells in condition) / (total # of cells
% in condition) * 100

% percent quiscient/activated in soft condition (manual grading)
soft_manual_q = (sum(soft)-sum(training_data_manual_grade(soft)))/sum(soft)*100;
soft_manual_a = sum(training_data_manual_grade(soft))/sum(soft)*100;

% percent quiscient/activated in stiff condition (manual grading)
stiff_manual_q = (sum(stiff)-sum(training_data_manual_grade(stiff)))/sum(stiff)*100;
stiff_manual_a = sum(training_data_manual_grade(stiff))/sum(stiff)*100;

% percent quiscient/activated in AZA condition (manual grading)
test_manual_q = ((length(test_data)-sum(manual_grade))/length(test_data))*100;
test_manual_a = (sum(manual_grade)/length(test_data))*100;

% percent quiscient/activated in soft condition (MLR grading)
soft_MLR_q = (sum(soft)-sum(training_MLR_grade(soft)))/sum(soft)*100;
soft_MLR_a = sum(training_MLR_grade(soft))/sum(soft)*100;

% percent quiscient/activated in stiff condition (MLR grading)
stiff_MLR_q = (sum(stiff)-sum(training_MLR_grade(stiff)))/sum(stiff)*100;
stiff_MLR_a = sum(training_MLR_grade(stiff))/sum(stiff)*100;

% percent quiscient/activated in AZA condition (MLR grading)
test_MLR_q = ((length(test_data)-sum(grade_MLR))/length(test_data))*100;
test_MLR_a = (sum(grade_MLR)/length(test_data))*100;

%% making pie charts

% Labels for the pie chart
labels = {'quiescent', 'activated'};

% Manually set RGB colors for each slice
colors = [0.80 0 0.80; 0 0.80 0];

% Soft manual grading
values = [soft_manual_q,soft_manual_a];
figure
p = pie(values,labels); hold on;
colormap(colors);
set(p(1:2:end), 'LineWidth', 16);
title('soft manual grading')
% exportgraphics(gcf,'soft_manual.png','Resolution',300)

% Stiff manual grading
values = [stiff_manual_q,stiff_manual_a];
figure
p = pie(values,labels);
colormap(colors);
set(p(1:2:end), 'LineWidth', 16);
title('stiff manual grading')
% exportgraphics(gcf,'stiff_manual.png','Resolution',300)

% Az manual grading
values = [test_manual_q,test_manual_a];
figure
p = pie(values,labels);
colormap(colors);
set(p(1:2:end), 'LineWidth', 16);
title('AZA manual grading')
% exportgraphics(gcf,'az_manual.png','Resolution',300)

% Soft MLR grading
values = [soft_MLR_q,soft_MLR_a];
figure
p = pie(values,labels);
colormap(colors);
set(p(1:2:end), 'LineWidth', 16);
title('soft MLR grading')
% exportgraphics(gcf,'soft_MLR_ar.png','Resolution',300)

% Stiff MLR grading
values = [stiff_MLR_q,stiff_MLR_a];
figure
p = pie(values,labels);
colormap(colors);
set(p(1:2:end), 'LineWidth', 16);
title('stiff MLR grading')
% exportgraphics(gcf,'stiff_MLR_ar.png','Resolution',300)

% Az MLR grading
values = [test_MLR_q,test_MLR_a];
figure
p = pie(values,labels);
colormap(colors);
set(p(1:2:end), 'LineWidth', 16);
title('AZA MLR grading')
% exportgraphics(gcf,'az_MLR_ar.png','Resolution',300)