%% MWE Cells PCA and eigenshape generation

clc
clear all
close all

% plot settings
set(gca,'DefaultLineLineWidth',8)
set(groot,'defaultAxesFontSize',40)
set(groot,'defaultAxesLineWidth',4)
xlims = [-75 75];
ylims = xlims;

% reads data
f_coeffs = readmatrix('cell_data.xlsx');

% removes NaNs (from Excel sheet headings and subheadings)
f_coeffs = f_coeffs(3:end,:);

% identifiers for cells seeded on soft and stiff hydrogels 
stiff_idx = f_coeffs(:,end)>0;
soft_idx = f_coeffs(:,end)<1;

% removes soft/stiff identifer
f_coeffs = f_coeffs(:,1:end-1);

% standarization of variables
mean_f_coeffs = mean(f_coeffs);
std_f_coeffs = std(f_coeffs);
f_coeffs_z = (f_coeffs - mean_f_coeffs)./std_f_coeffs;

% PCA
A = cov(f_coeffs_z); % co-variance matrix 
[V,D] = eig(A); % Get Eigenvalues and Eigenvectors 
Eig = diag(D); % concatenates Eigenvalues
[val,idx] = sort(Eig,'descend'); % sorts Eigenvalues
PV = Eig(idx); % sorts Eigenvalues
PC = V(:,idx); % sorts Eigenvectors

% PC scores
score = f_coeffs_z*PC;

% Create color bar
mg_r = linspace(1.0,0,1000)';
mg_g = linspace(0,1.0,1000)';
mg_b = linspace(1.0,0,1000)';
colors = [mg_r, mg_g, mg_b];

% query points for Fourier shape evaluation
thetas = linspace(0,2*pi,10000);

% Plotting of Eigenshapes from -2*STD to 2*STD on PC axis 1-3
for k = 1:3 % iterates for PC axes 1-3

    % starts at -2*STD
    score2 = mean(score);
    score2_std = std(score);
    score2(k) = score2(k) - 2*score2_std(k);
    
    for i = 1:5 % Iterates from -2*STD to 2*STD
        
        % Adds STD at each iteration 
        if i > 1
            score2(k) = score2(k) + score2_std(k);
        end
    
        % Computes recovered coefficients 
        PC_inv = inv(PC);
        rec_coeff = score2*PC_inv;
        rec_coeff = (rec_coeff.*std_f_coeffs)+mean_f_coeffs; % converts from standarized values back to absolute values
        rec_DOA = rec_coeff(:,end); % concatenates DOA value
        rec_coeff(:,end) = []; % removes DOA value

        % sets up coloring of shapes based on DOA
        color_indx2 = round(rec_DOA*1000);
        color_indx2(color_indx2 > 1000) = 1000;
        color_indx2(color_indx2 < 1) = 1;

        % Concatenates coefficients 
        a_CELL = rec_coeff(1:31); % a0 - a30
        b_CELL = rec_coeff(32:62); % b0 - b30
        a_DAPI = rec_coeff(63:73); % a0 - a10
        b_DAPI = rec_coeff(74:84); % b0 - b10
    
        % Evaluate Fourier Series with recovered coefficients (cell body) 
        [x_CELL] = fourier_series_evaluate(a_CELL,thetas);
        [y_CELL] = fourier_series_evaluate(b_CELL,thetas);
        pgon_CELL = polyshape(x_CELL,y_CELL);
    
        % Evaluate Fourier Series with recovered coefficients (nuclei) 
        [x_DAPI] = fourier_series_evaluate(a_DAPI,thetas);
        [y_DAPI] = fourier_series_evaluate(b_DAPI,thetas);
        pgon_DAPI = polyshape(x_DAPI,y_DAPI);
    
        % Plots shapes at PC Axis 
        figure (1)
        plot(pgon_CELL,'FaceColor',colors(color_indx2,:),'EdgeColor','k','LineStyle','-','FaceAlpha',1,'LineWidth',8); hold on;
        plot(pgon_DAPI,'FaceColor','b','EdgeColor','k','LineStyle','-','FaceAlpha',1,'LineWidth',8); hold off;
        xlim(xlims)
        ylim(ylims)
        set(gca,'XTick',[],'YTick',[])
        set(gca,'xcolor','k','ycolor','k')
        title(strcat('Axis',num2str(k),',Shape',num2str(i)))
        box on
        pbaspect([1 1 1])
        exportgraphics(gcf,strcat('CELL_and_DAPI_eigenshape_PC_',num2str(k),'_SD_',num2str(i),'.png'),'Resolution',300)
    
    end
    
    % pre-allocate x-values for kernel density estimation
    pts = linspace(min(score(:,k)),max(score(:,k)),100);
    
    % kernel density estimation of PC scores between soft and stiff
    figure (500+k)
    [f_soft,xi_soft] = ksdensity(score(soft_idx,k),pts); % KDE for soft
    plot(xi_soft,f_soft,'k','LineWidth',16); hold on; % Plot for soft
    [f_stiff,xi_stiff] = ksdensity(score(stiff_idx,k),pts); % KDE for stiff
    plot(xi_stiff,f_stiff,'r','LineWidth',16); hold off; % Plot for stiff
    xlabel(strcat("PC",num2str(k)),'FontWeight','bold')
    ylabel("PDE",'FontWeight','bold')
    set(get(gca, 'XAxis'), 'FontWeight', 'bold');
    set(get(gca, 'YAxis'), 'FontWeight', 'bold');
    pbaspect([1 1 1])
    exportgraphics(gcf,strcat('SOFT_and_STIFF_KDE_PCA_',num2str(k),'.png'),'Resolution',300)

end

