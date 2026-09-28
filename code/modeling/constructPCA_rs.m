%% Run PCA for separation of random samples between WTs and 2E2s
% 24 Jan 2025
% Cyriel Huijer
% Last updated = 22 Oct 2025 (correct run was added in path, this matches
% the PCA in the presentation 20250127_updateMeeting_validation.pptx

clc
clear

% Load in rs solutions:
load("C:/Users/chuijer/surfdrive/PhD/Projects/WesternOntario_collab/6_results/metabolic_modeling/20250120_rsWithRosemary/run4/20250122_WT_solutions_2K_rerun.mat");
load("C:/Users/chuijer/surfdrive/PhD/Projects/WesternOntario_collab/6_results/metabolic_modeling/20250120_rsWithRosemary/run4/20250122_KO_solutions_2K_rerun.mat");
load("iCHOv1.mat");

% Load in res. file: (from run 4)
res = readtable("20250122_rs_results_after_filter.csv");


% PCA of random samples:
data_WT = full(solutions_WT);  % WT data (6663 x 2000)
data_KO = full(solutions_KO);  % KO data (6663 x 2000)
rxn_idx = findRxnIDs(iCHOv1,res.rxn);
% select significant rows from table
data_WT = data_WT(rxn_idx,:);
data_KO = data_KO(rxn_idx,:);


combined_data = [data_WT, data_KO]';  % Transpose to make samples as rows (4000 x 121)
% Create labels for conditions
labels = [repmat({'WT'}, 2000, 1); repmat({'KO'}, 2000, 1)];

% Perform PCA
[coeff, score, ~, ~, explained] = pca(combined_data);

% Plot PCA (PC1 vs PC2)
figure;
colors = [hex2rgb('#FF7043'); hex2rgb('#0070C0')]; % KO = blue, WT = orange
gscatter(score(:, 1), score(:, 2), labels, colors, 'xo', 10);
xlabel(sprintf('PC1 (%.2f%% Variance Explained)', explained(1)));
ylabel(sprintf('PC2 (%.2f%% Variance Explained)', explained(2)));
title('PCA of Random Samples (WT vs KO)');
legend('WT', '2E2');
grid on;

%exportgraphics(gcf, 'PCA_plot.pdf', 'ContentType', 'vector');




