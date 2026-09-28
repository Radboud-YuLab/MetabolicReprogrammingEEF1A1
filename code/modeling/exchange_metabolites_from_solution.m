cd /code

load('modeling/models/iCHOv1.mat');
model = iCHOv1;

[exchangeRxns, exchIndex] = getExchangeRxns(model);

% load solutions:

load('data/modeling/solutions/20250122_WT_solutions_2K_rerun.mat');
load('data/modeling/solutions/20250122_KO_solutions_2K_rerun.mat');

solutions_WT_exchRxn = solutions_WT(exchIndex, :);
solutions_KO_exchRxn = solutions_KO(exchIndex, :);

WT_mean = mean(solutions_WT_exchRxn, 2);
KO_mean = mean(solutions_KO_exchRxn, 2);

overview = table(exchangeRxns(:), WT_mean, KO_mean,'VariableNames', {'ExchangeRxn', 'WT_mean', 'KO_mean'});

% filter reactions with flux smaller than 1e-5:
filter_criteria = ~(abs(overview.WT_mean) < 1e-7 & abs(overview.KO_mean) < 1e-7);

overview_filtered = overview(filter_criteria, :);

% Preallocate cell array to hold metabolites per reaction
metabolites_per_rxn = cell(length(exchangeRxns), 1);

for k = 1:length(exchangeRxns)
    rxnID = exchangeRxns{k};
    rxnIndex = find(strcmp(model.rxns, rxnID));
    
    % Find metabolites involved (non-zero stoichiometry)
    metIndices = find(model.S(:, rxnIndex) ~= 0);
    metsInReaction = model.mets(metIndices);
    
    % Store just the metabolite IDs (no stoich)
    metabolites_per_rxn{k} = metsInReaction;
end

% Optionally create a table
metaboliteTable = table(exchangeRxns, metabolites_per_rxn, ...
    'VariableNames', {'ExchangeRxn', 'met'});

metaboliteTable = metaboliteTable(filter_criteria, :);

final_table = table(metaboliteTable.met, overview_filtered.WT_mean,overview_filtered.KO_mean, ...
    'VariableNames', {'met', 'WT_mean','KO_mean'});

writetable(final_table, 'data/modeling/exchange_fluxes_solution.csv');