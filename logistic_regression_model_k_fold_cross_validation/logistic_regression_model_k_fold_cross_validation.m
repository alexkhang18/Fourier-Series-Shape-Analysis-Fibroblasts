%% logistic_regression_model_k_fold_cross_validation

clc
clear all
close all

%% import data

DATA = readmatrix('DATA.xlsx');

% concatenates activation data and converts unactivated cells to have a
% label of 2 and activated cells to have a label of 1
activation = DATA(:,13);
sp = logical(activation);
sp = ~sp + 1;

%% logistic regression models

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cell area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_area = DATA(:,1);

% fits the logistic regresison model
B_cell_area = mnrfit(cell_area,sp);

% computes model parameters mu and s
B_mu_cell_area = -B_cell_area(1)/B_cell_area(2);
B_s_cell_area = 1/B_cell_area(2);

% computes x value at which cell has 99% probability of being activated
x_99 = (-log(1/0.99 - 1) - B_cell_area(1)) / B_cell_area(2);

% generates x and y values of logistic regression model for plotting
x_LR = linspace(0,max([x_99;cell_area]),1000);
y_LR = 1./(1+exp(-(B_cell_area(1)+B_cell_area(2)*x_LR)));

% finds x,y location where activaton probability is 50%
star = [B_mu_cell_area,1./(1+exp(-(B_cell_area(1)+B_cell_area(2)*B_mu_cell_area)))];

% computes negative log-likelihood 
activation_probability = 1./(1+exp(-(B_cell_area(1)+B_cell_area(2)*cell_area))); 
nll_cell_area = sum(-activation.*log(activation_probability)-(1-activation).*log(1-activation_probability));

% plots logistic regression model
figure
plot(x_LR,y_LR,'LineWidth',16); hold on;
scatter(star(1),star(2),1000,'pentagram','r','filled'); hold on;
xlabel('cell area (\mum^2)','FontWeight','bold')
ylabel({'activation';'probability'},'FontWeight','bold')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
title('cell area')
% exportgraphics(gcf,'cell_area_logistic_regression.png','Resolution',300)

% computes accuracy of logistic regression model to perform binary grading
% percent correct = (total cells - incorrectly graded cells) / total cells * 100 
predicted_activation = activation_probability>0.5; % binary grading
percent_correct_cell_area = ((length(activation)-sum(nonzeros(abs(predicted_activation-activation))))/length(activation))*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nuclear area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nuclear_area = DATA(:,7);

% fits the logistic regresison model
B_nuclear_area = mnrfit(nuclear_area,sp);

% computes model parameters mu and s
B_mu_nuclear_area = -B_nuclear_area(1)/B_nuclear_area(2);
B_s_nuclear_area = 1/B_nuclear_area(2);

% computes x value at which cell has 99% probability of being activated
x_99 = (-log(1/0.99 - 1) - B_nuclear_area(1)) / B_nuclear_area(2);

% generates x and y values of logistic regression model for plotting
x_LR = linspace(0,max([x_99;nuclear_area]),1000);
y_LR = 1./(1+exp(-(B_nuclear_area(1)+B_nuclear_area(2)*x_LR)));

% finds x,y location where activaton probability is 50%
star = [B_mu_nuclear_area,1./(1+exp(-(B_nuclear_area(1)+B_nuclear_area(2)*B_mu_nuclear_area)))];

% computes negative log-likelihood 
activation_probability = 1./(1+exp(-(B_nuclear_area(1)+B_nuclear_area(2)*nuclear_area))); 
nll_nuclear_area = sum(-activation.*log(activation_probability)-(1-activation).*log(1-activation_probability));

% plots logistic regression model
figure
plot(x_LR,y_LR,'LineWidth',16); hold on;
scatter(star(1),star(2),1000,'pentagram','r','filled'); hold on;
xlabel('nuclear area (\mum^2)','FontWeight','bold')
ylabel({'activation';'probability'},'FontWeight','bold')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
title('nuclear area')
% exportgraphics(gcf,'nuclear_area_logistic_regression.png','Resolution',300)

% computes accuracy of logistic regression model to perform binary grading
% percent correct = (total cells - incorrectly graded cells) / total cells * 100 
predicted_activation = activation_probability>0.5; % binary grading
percent_correct_nuclear_area = ((length(activation)-sum(nonzeros(abs(predicted_activation-activation))))/length(activation))*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cell aspect ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_aspect_ratio = DATA(:,5);

% fits the logistic regresison model
B_cell_aspect_ratio = mnrfit(cell_aspect_ratio,sp);

% computes model parameters mu and s
B_mu_cell_aspect_ratio = -B_cell_aspect_ratio(1)/B_cell_aspect_ratio(2);
B_s_cell_aspect_ratio = 1/B_cell_aspect_ratio(2);

% computes x value at which cell has 99% probability of being activated
x_99 = (-log(1/0.99 - 1) - B_cell_aspect_ratio(1)) / B_cell_aspect_ratio(2);

% generates x and y values of logistic regression model for plotting
x_LR = linspace(0,max([x_99;cell_aspect_ratio]),1000);
y_LR = 1./(1+exp(-(B_cell_aspect_ratio(1)+B_cell_aspect_ratio(2)*x_LR)));

% finds x,y location where activaton probability is 50%
star = [B_mu_cell_aspect_ratio,1./(1+exp(-(B_cell_aspect_ratio(1)+B_cell_aspect_ratio(2)*B_mu_cell_aspect_ratio)))];

% computes negative log-likelihood 
activation_probability = 1./(1+exp(-(B_cell_aspect_ratio(1)+B_cell_aspect_ratio(2)*cell_aspect_ratio))); 
nll_cell_aspect_ratio = sum(-activation.*log(activation_probability)-(1-activation).*log(1-activation_probability));

% plots logistic regression model
figure
plot(x_LR,y_LR,'LineWidth',16); hold on;
scatter(star(1),star(2),1000,'pentagram','r','filled'); hold on;
xlabel('cell aspect ratio (\mum^2)','FontWeight','bold')
ylabel({'activation';'probability'},'FontWeight','bold')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
title('cell aspect ratio')
% exportgraphics(gcf,'cell_aspect_ratio_logistic_regression.png','Resolution',300)

% computes accuracy of logistic regression model to perform binary grading
% percent correct = (total cells - incorrectly graded cells) / total cells * 100 
predicted_activation = activation_probability>0.5; % binary grading
percent_correct_cell_aspect_ratio = ((length(activation)-sum(nonzeros(abs(predicted_activation-activation))))/length(activation))*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nuclear aspect ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nuclear_aspect_ratio = DATA(:,11);

% fits the logistic regresison model
B_nuclear_aspect_ratio = mnrfit(nuclear_aspect_ratio,sp);

% computes model parameters mu and s
B_mu_nuclear_aspect_ratio = -B_nuclear_aspect_ratio(1)/B_nuclear_aspect_ratio(2);
B_s_nuclear_aspect_ratio = 1/B_nuclear_aspect_ratio(2);

% computes x value at which cell has 99% probability of being activated
x_99 = (-log(1/0.99 - 1) - B_nuclear_aspect_ratio(1)) / B_nuclear_aspect_ratio(2);

% generates x and y values of logistic regression model for plotting
x_LR = linspace(0,max([x_99;nuclear_aspect_ratio]),1000);
y_LR = 1./(1+exp(-(B_nuclear_aspect_ratio(1)+B_nuclear_aspect_ratio(2)*x_LR)));

% finds x,y location where activaton probability is 50%
star = [B_mu_nuclear_aspect_ratio,1./(1+exp(-(B_nuclear_aspect_ratio(1)+B_nuclear_aspect_ratio(2)*B_mu_nuclear_aspect_ratio)))];

% computes negative log-likelihood 
activation_probability = 1./(1+exp(-(B_nuclear_aspect_ratio(1)+B_nuclear_aspect_ratio(2)*nuclear_aspect_ratio))); 
nll_nuclear_aspect_ratio = sum(-activation.*log(activation_probability)-(1-activation).*log(1-activation_probability));

% plots logistic regression model
figure
plot(x_LR,y_LR,'LineWidth',16); hold on;
scatter(star(1),star(2),1000,'pentagram','r','filled'); hold on;
xlabel('nuclear aspect ratio (\mum^2)','FontWeight','bold')
ylabel({'activation';'probability'},'FontWeight','bold')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
title('nuclear aspect ratio')
% exportgraphics(gcf,'nuclear_aspect_ratio_logistic_regression.png','Resolution',300)

% computes accuracy of logistic regression model to perform binary grading
% percent correct = (total cells - incorrectly graded cells) / total cells * 100 
predicted_activation = activation_probability>0.5; % binary grading
percent_correct_nuclear_aspect_ratio = ((length(activation)-sum(nonzeros(abs(predicted_activation-activation))))/length(activation))*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Degree of Anisotropy (DOA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DOA = DATA(:,14);

% fits the logistic regresison model
B_DOA = mnrfit(DOA,sp);

% computes model parameters mu and s
B_mu_DOA = -B_DOA(1)/B_DOA(2);
B_s_DOA = 1/B_DOA(2);

% computes x value at which cell has 99% probability of being activated
x_99 = (-log(1/0.99 - 1) - B_DOA(1)) / B_DOA(2);

% generates x and y values of logistic regression model for plotting
x_LR = linspace(0,1,1000);
y_LR = 1./(1+exp(-(B_DOA(1)+B_DOA(2)*x_LR)));

% finds x,y location where activaton probability is 50%
star = [B_mu_DOA,1./(1+exp(-(B_DOA(1)+B_DOA(2)*B_mu_DOA)))];

% computes negative log-likelihood 
activation_probability = 1./(1+exp(-(B_DOA(1)+B_DOA(2)*DOA))); 
nll_DOA = sum(-activation.*log(activation_probability)-(1-activation).*log(1-activation_probability));

% plots logistic regression model
figure
plot(x_LR,y_LR,'LineWidth',16); hold on;
scatter(star(1),star(2),1000,'pentagram','r','filled'); hold on;
xlabel('DOA','FontWeight','bold')
ylabel({'activation';'probability'},'FontWeight','bold')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
title('Degree of Anisotropy')
% exportgraphics(gcf,'DOA_logistic_regression.png','Resolution',300)

% computes accuracy of logistic regression model to perform binary grading
% percent correct = (total cells - incorrectly graded cells) / total cells * 100 
predicted_activation = activation_probability>0.5; % binary grading
percent_correct_DOA = ((length(activation)-sum(nonzeros(abs(predicted_activation-activation))))/length(activation))*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Multi-variable logistic regression (MLR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create color scale
mg_r = linspace(1.0,0,1000)';
mg_g = linspace(0,1.0,1000)';
mg_b = linspace(1.0,0,1000)';
colors = [mg_r, mg_g, mg_b];

% concatenates all measures for MLR
meas = [DOA,cell_area,nuclear_area,cell_aspect_ratio,nuclear_aspect_ratio];

% fits MLR
B_mlr = mnrfit(meas,sp);

% computes mean of all measures
mean_meas = mean(meas);

% generates XYZ values of MLR model for plotting
x = linspace(0,max(cell_area));
y = linspace(0,1);
[X,Y] = meshgrid(x,y);
Z = 1./(1+exp(-(B_mlr(1)+B_mlr(2)*Y+B_mlr(3)*X+B_mlr(4)*mean_meas(3)+B_mlr(5)*mean_meas(4)+B_mlr(6)*mean_meas(5))));

% computes negative log-likelihood 
MLR_activation_probability = 1./(1+exp(-(B_mlr(1)+B_mlr(2)*meas(:,1)+B_mlr(3)*meas(:,2)+B_mlr(4)*meas(:,3)+B_mlr(5)*meas(:,4)+B_mlr(6)*meas(:,5))));
nll_MLR = sum(-activation.*log(MLR_activation_probability)-(1-activation).*log(1-MLR_activation_probability ));

% plots MLR model as countour plot
fig = figure;
contourf(X,Y,Z,1000,'LineColor','none')
xticks([0 4500 9000])
xlabel('cell area (\mum^2)')
ylabel('DOA')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
colormap(colors)
cb = colorbar(); 
cb.LineWidth = 4;
cb.Limits =[0 1];
ylabel(cb,'activation probability','FontSize',30,'Rotation',270,'FontWeight','bold')
% exportgraphics(gcf,'MLR_logistic_regression.png','Resolution',300)

% computes accuracy of logistic regression model to perform binary grading
% percent correct = (total cells - incorrectly graded cells) / total cells * 100 
predicted_activation = MLR_activation_probability>0.5;
percent_correct_MLR = ((length(activation)-sum(nonzeros(abs(predicted_activation-activation))))/length(activation))*100;

%% receiver operating characteristic curves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cell area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[TPR_cell_area,FPR_cell_area] = ROC_curve(cell_area,activation); %generates ROC curve
AUC_cell_area = abs(trapz(FPR_cell_area,TPR_cell_area)); %area under the curve (AUC)

% plot ROC curve
figure
plot([0 1],[0 1], 'LineWidth',8,'LineStyle',':','Color',[0.5,0.5,0.5]); hold on;
plot(FPR_cell_area,TPR_cell_area,'LineWidth',16,'Color','k'); hold off;
xlabel('FPR','FontWeight','bold')
ylabel('TPR','FontWeight','bold')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
title('cell area')
% exportgraphics(gcf,'cell_area_ROC.png','Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nuclear area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[TPR_nuclear_area,FPR_nuclear_area] = ROC_curve(nuclear_area,activation); %generates ROC curve
AUC_nuclear_area = abs(trapz(FPR_nuclear_area,TPR_nuclear_area)); %area under the curve (AUC)

% plot ROC curve
figure
plot([0 1],[0 1], 'LineWidth',8,'LineStyle',':','Color',[0.5,0.5,0.5]); hold on;
plot(FPR_nuclear_area,TPR_nuclear_area,'LineWidth',16,'Color','k'); hold off;
xlabel('FPR','FontWeight','bold')
ylabel('TPR','FontWeight','bold')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
title('nuclear area')
% exportgraphics(gcf,'nuc_area_ROC.png','Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cell aspect ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[TPR_cell_aspect_ratio,FPR_cell_aspect_ratio] = ROC_curve(cell_aspect_ratio,activation); %generates ROC curve
AUC_cell_aspect_ratio = abs(trapz(FPR_cell_aspect_ratio,TPR_cell_aspect_ratio)); %area under the curve (AUC)

% plot ROC curve
figure
plot([0 1],[0 1], 'LineWidth',8,'LineStyle',':','Color',[0.5,0.5,0.5]); hold on;
plot(FPR_cell_aspect_ratio,TPR_cell_aspect_ratio,'LineWidth',16,'Color','k'); hold off;
xlabel('FPR','FontWeight','bold')
ylabel('TPR','FontWeight','bold')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
title('cell aspect ratio')
% exportgraphics(gcf,'cell_aspect_ratio_ROC.png','Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nuclear aspect ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[TPR_nuclear_aspect_ratio,FPR_nuclear_aspect_ratio] = ROC_curve(nuclear_aspect_ratio,activation); %generates ROC curve
AUC_nuclear_aspect_ratio = abs(trapz(FPR_nuclear_aspect_ratio,TPR_nuclear_aspect_ratio)); %area under the curve (AUC)

% plot ROC curve
figure
plot([0 1],[0 1], 'LineWidth',8,'LineStyle',':','Color',[0.5,0.5,0.5]); hold on;
plot(FPR_nuclear_aspect_ratio,TPR_nuclear_aspect_ratio,'LineWidth',16,'Color','k'); hold off;
xlabel('FPR','FontWeight','bold')
ylabel('TPR','FontWeight','bold')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
title('nuclear aspect ratio')
% exportgraphics(gcf,'nuclear_aspect_ratio_ROC.png','Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DOA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[TPR_DOA,FPR_DOA] = ROC_curve(DOA,activation); %generates ROC curve
AUC_DOA = abs(trapz(FPR_DOA,TPR_DOA)); %area under the curve (AUC)

% plot ROC curve
figure
plot([0 1],[0 1], 'LineWidth',8,'LineStyle',':','Color',[0.5,0.5,0.5]); hold on;
plot(FPR_DOA,TPR_DOA,'LineWidth',16,'Color','k'); hold off;
xlabel('FPR','FontWeight','bold')
ylabel('TPR','FontWeight','bold')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
title('DOA')
% exportgraphics(gcf,'DOA_ROC.png','Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MLR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[TPR_MLR,FPR_MLR] = ROC_curve(MLR_activation_probability,activation); %generates ROC curve
AUC_MLR = abs(trapz(FPR_MLR,TPR_MLR)); %area under the curve (AUC)

% plot ROC curve
figure
plot([0 1],[0 1], 'LineWidth',8,'LineStyle',':','Color',[0.5,0.5,0.5]); hold on;
plot(FPR_MLR,TPR_MLR,'LineWidth',16,'Color','k'); hold off;
xlabel('FPR','FontWeight','bold')
ylabel('TPR','FontWeight','bold')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
title('MLR')
% exportgraphics(gcf,'MLR_ROC.png','Resolution',300)

%% k-fold cross validation

% sets seed for random generator for consistent results
rng(7)

% concatenates all measures
meas = [DOA,cell_area,nuclear_area,cell_aspect_ratio,nuclear_aspect_ratio];

% generateds random indices to shuffle meas
ran_idx = randperm(size(DATA,1),size(DATA,1));

% meas2 is a shuffled version of meas for random selection of
% training/validation datasets
meas2 = meas(ran_idx,:);

% activation2 is a shuffled version of activation for consistency with
% meas2
activation2 = activation(ran_idx,:);

% converts unactivated cells to have a label of 2 and activated cells to have a label of 1
sp2 = sp(ran_idx);

% k-fold cross validation with k = 10
for k = 1:10

% decides which indices are in the validation dataset
if k < 10
    val_idx = round(size(DATA,1)/10)*(k-1)+1:round(size(DATA,1)/10)*k;
else
    val_idx = round(size(DATA,1)/10)*(k-1)+1:size(DATA,1);
end

% concatenates validation dataset
val_data = meas2(val_idx,:);

% concatenates training dataset
train_data = meas2;
train_data(val_idx,:) = []; % deletes validation dataset from training dataset

% converts unactivated cells to have a label of 2 and activated cells to have a label of 1
sp_temp = sp2;
sp_temp(val_idx) = []; % deletes validation dataset labels

% DOA 
DOA_train = train_data(:,1); % training data
DOA_test = val_data(:,1); % validation data
[B,dev,stats] = mnrfit(DOA_train,sp_temp); % fits the logistic regresison model
activation_probability = 1./(1+exp(-(B(1)+B(2)*DOA_test))); % predicted activation probability
predicted_activation2 = activation_probability>0.5; % binary prediction
percent_correct_test_DOA(k) = ((length(activation2(val_idx))-sum(nonzeros(abs(predicted_activation2-activation2(val_idx)))))/length(activation2(val_idx)))*100;

% cell area
cell_area_train = train_data(:,2);% training data
cell_area_test = val_data(:,2);% validation data
[B,dev,stats] = mnrfit(cell_area_train,sp_temp);% fits the logistic regresison model
activation_probability = 1./(1+exp(-(B(1)+B(2)*cell_area_test)));% predicted activation probability
predicted_activation2 = activation_probability>0.5; % binary prediction
percent_correct_test_cell_area(k) = ((length(activation2(val_idx))-sum(nonzeros(abs(predicted_activation2-activation2(val_idx)))))/length(activation2(val_idx)))*100;

% nuc area
nuc_area_train = train_data(:,3);% training data
nuc_area_test = val_data(:,3);% validation data
[B,dev,stats] = mnrfit(nuc_area_train,sp_temp);% fits the logistic regresison model
activation_probability = 1./(1+exp(-(B(1)+B(2)*nuc_area_test)));% predicted activation probability
predicted_activation2 = activation_probability>0.5; % binary prediction
percent_correct_test_nuc_area(k) = ((length(activation2(val_idx))-sum(nonzeros(abs(predicted_activation2-activation2(val_idx)))))/length(activation2(val_idx)))*100;

% cell aspect ratio
cell_aspect_ratio_train = train_data(:,4);% training data
cell_aspect_ratio_test = val_data(:,4);% validation data
[B,dev,stats] = mnrfit(cell_aspect_ratio_train,sp_temp);% fits the logistic regresison model
activation_probability = 1./(1+exp(-(B(1)+B(2)*cell_aspect_ratio_test)));% predicted activation probability
predicted_activation2 = activation_probability>0.5; % binary prediction
percent_correct_test_cell_aspect_ratio(k) = ((length(activation2(val_idx))-sum(nonzeros(abs(predicted_activation2-activation2(val_idx)))))/length(activation2(val_idx)))*100;

% nuc aspect ratio
nuc_aspect_ratio_train = train_data(:,5);% training data
nuc_aspect_ratio_test = val_data(:,5);% validation data
[B,dev,stats] = mnrfit(nuc_aspect_ratio_train,sp_temp);% fits the logistic regresison model
activation_probability = 1./(1+exp(-(B(1)+B(2)*nuc_aspect_ratio_test)));% predicted activation probability
predicted_activation2 = activation_probability>0.5; % binary prediction
percent_correct_test_nuc_aspect_ratio(k) = ((length(activation2(val_idx))-sum(nonzeros(abs(predicted_activation2-activation2(val_idx)))))/length(activation2(val_idx)))*100;

% mlr
mlr_train = train_data;% training data
mlr_test = val_data;% validation data
[B,dev,stats] = mnrfit(mlr_train,sp_temp);% fits the logistic regresison model
activation_probability = 1./(1+exp(-(B(1)+B(2)*mlr_test(:,1)+B(3)*mlr_test(:,2)+B(4)*mlr_test(:,3)+B(5)*mlr_test(:,4)+B(6)*mlr_test(:,5))));% predicted activation probability
predicted_activation2 = activation_probability>0.5; % binary prediction
percent_correct_test_MLR(k) = ((length(activation2(val_idx))-sum(nonzeros(abs(predicted_activation2-activation2(val_idx)))))/length(activation2(val_idx)))*100;

end

% compiles k-fold cross validation results to be plotted
data_to_plot = [percent_correct_test_cell_area',percent_correct_test_nuc_area',percent_correct_test_cell_aspect_ratio',percent_correct_test_nuc_aspect_ratio',percent_correct_test_DOA',percent_correct_test_MLR'];

% generates random x-values to superimpose data points on top of bar plot
x_scatter_data = [ones(length(percent_correct_test_cell_area),1).*(1+(rand(length(percent_correct_test_cell_area),1)-0.5)/2);...
    1+ones(length(percent_correct_test_nuc_area),1).*(1+(rand(length(percent_correct_test_nuc_area),1)-0.5)/2);...
    2+ones(length(percent_correct_test_cell_aspect_ratio),1).*(1+(rand(length(percent_correct_test_cell_aspect_ratio),1)-0.5)/2);...
    3+ones(length(percent_correct_test_nuc_aspect_ratio),1).*(1+(rand(length(percent_correct_test_nuc_aspect_ratio),1)-0.5)/2);...
    4+ones(length(percent_correct_test_DOA),1).*(1+(rand(length(percent_correct_test_DOA),1)-0.5)/2);...
    5+ones(length(percent_correct_test_MLR),1).*(1+(rand(length(percent_correct_test_MLR),1)-0.5)/2)];

% one way anova on accuracy of logistic regression models
[~,~,stats] = anova1(data_to_plot);
[c,~,~,gnames] = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A") = gnames(tbl.("Group A"));
tbl.("Group B") = gnames(tbl.("Group B"))

% generates x, y, and std for bar plot
x = 1:size(data_to_plot,2);
data = mean(data_to_plot);
err = std(data_to_plot);

% Create bar plot
figure;
h = bar(x, data);
h.FaceColor = [0.6 0.6 0.6];
h.EdgeColor = 'k';
xlim([0 7])
ylim([40 100])
hold on;
% Add error bars
er = errorbar(x, data, err, err);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 4;
er.CapSize = 24;
% Scatter the underlying data points
scatter(x_scatter_data,data_to_plot(:),200,'r.'); hold off;
ylabel('% correct','FontWeight','bold')
title("k-fold cross-validation")
xticklabels({'cell area','nuc area','cell AR','nuc AR','DOA','MLR'})
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(get(gca, 'Title'), 'FontWeight', 'bold');
set(gca,'FontSize',25);
set(h,'LineWidth', 4);
pbaspect([2 1 1])
xlim([0.25 6.75])
% exportgraphics(gcf,'k_fold_cross_validation_bars.png','Resolution',300)
