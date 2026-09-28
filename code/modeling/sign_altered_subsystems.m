%% Get statistically sign reactions table + plot
% 24 Jan 2025
% Cyriel Huijer
% Last updated = 22 Oct 2025 (correct run was added in path, this matches
% the sign diff plot in the presentation 20250127_updateMeeting_validation.pptx

clear
clc

% set wd to code 
cd /code

% Load in rs solutions:
load("data/modeling/solutions/20250122_WT_solutions_2K_rerun.mat");
load("data/modeling/solutions/20250122_KO_solutions_2K_rerun.mat");
load("data/modeling/WT_model.mat")
load("data/modeling/KO_model.mat")

% set out dir root for function
out_dir_root = "";

WT_model.id = 'WT_model';
KO_model.id = 'KO_model';

res = evaluateSampling_updated(WT_model, KO_model, solutions_WT, solutions_KO, 0.01,0,out_dir_root);
filter_rows = abs(res.mean_sample_WT_model) > 1e-5 & abs(res.mean_sample_KO_model) > 1e-5;
res = res(filter_rows, :);
res = res(abs(res.raw_diff) > 0.3, :);
alt_subs = countAltSubsystems_CH_2(WT_model,KO_model,res,out_dir_root);
writetable(alt_subs, 'data/modeling/alt_subs_updated.csv');

