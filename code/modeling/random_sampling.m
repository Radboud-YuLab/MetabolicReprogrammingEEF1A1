% script that does random sampling for CHO 2E2 manuscript
% written by: Cyriel Huijer
% last updated: 22 June 2026

initCobraToolbox;
changeCobraSolver('gurobi');
setRavenSolver('gurobi');

%%
clear
clc

cd C:\Users\chuijer\Documents\nextcloud_backup\PhD\Projects\WesternOntario_collab\9_manuscript\code

% Load model
load("modeling/models/iCHOv1.mat")
random_sampling = "no";
% Change model name to model
model = iCHOv1;

% Set variable for media components
mediaComps = {'L-Alanine', 'L-Arginine', 'L-Asparagine', 'L-Aspartate', 'L-Cysteine', ...
    'L-Glutamine', 'L-Glutamate', 'Glycine', 'L-Histidine', 'L-Lysine', ...
    'L-Methionine', 'L-Phenylalanine', 'L-Proline', 'L-Serine', 'L-Threonine', ...
    'L-Tryptophan', 'L-Tyrosine', 'L-Valine', 'D-Glucose', ...
    'L-Leucine', 'L-Isoleucine', 'Alpha-Tocopherol', 'Aquacob(III)alamin', ...
    'Fe2+ mitochondria', 'Gamma-Tocopherol', 'H2O H2O', 'Hypoxanthine', ...
    'Myo-Inositol', 'Linoleic acid (all cis C18:2) n-6', 'Alpha-Linolenic acid, C18:3, n-3', ...
    'Lipoate', 'Nicotinamide', 'O2 O2', '(R)-Pantothenate', 'Phosphate', ...
    'Pyridoxine', 'Retinoate', 'Riboflavin C17H20N4O6', 'Sulfate', 'Thiamin', ...
    'Thymidine C10H14N2O5', 'Folate', 'Choline C5H14NO', 'Biotin'}';

% Only allow uptake + release for media components
model = setExchangeBounds(model,mediaComps,-1000,1000,true);
clear mediaComps
% without SK_Tyr_ggn_c the biomass production = 0
model = setParam(model,'lb','SK_Tyr_ggn_c',-0.100); % lb of -0.100 is standard value in model. 

% Load in measured exchange fluxes:
exch_ub = readtable('data/exchange_fluxes/for_modeling/exch_fluxes_ub_lm_sign_diff_metNames_AlaRemoved.csv');
exch_lb = readtable('data/exchange_fluxes/for_modeling/exch_fluxes_lb_lm_sign_diff_metNames_AlaRemoved.csv');

% Conversion to mmol/gDW/h by division with the cell dry weight (taken from
% https://www.researchgate.net/figure/Biomass-composition-of-various-CHO-cell-lines-under-defined-conditions-Error-bars_fig1_342579887)
exch_ub_mmol_gDW_h = exch_ub;
exch_lb_mmol_gDW_h = exch_lb;
exch_ub_mmol_gDW_h.KO = exch_ub.KO * (1 / (264 * 10^-12));
exch_lb_mmol_gDW_h.KO = exch_lb.KO * (1 / (264 * 10^-12));
exch_ub_mmol_gDW_h.WT = exch_ub.WT * (1 / (264 * 10^-12));
exch_lb_mmol_gDW_h.WT = exch_lb.WT * (1 / (264 * 10^-12));

clear exch_ub exch_lb

% Constrain model with exchange fluxes from both models
KO_model = setExchangeBounds(model, ...
                                  exch_lb_mmol_gDW_h.metabolite, ...
                                  exch_lb_mmol_gDW_h.KO, ...
                                  exch_ub_mmol_gDW_h.KO, ...
                                  false);
WT_model = setExchangeBounds(model, ...
                                  exch_lb_mmol_gDW_h.metabolite, ...
                                  exch_lb_mmol_gDW_h.WT, ...
                                  exch_ub_mmol_gDW_h.WT, ...
                                  false);
clear model iCHOv1

%printConstraints(KO_model, -999, 999)
sol = optimizeCbModel(KO_model); %0.8581
%disp(sol.f);
sol = optimizeCbModel(WT_model); %0.8581
%disp(sol.f);
clear sol


% Add growth constraint to models:
biomass_rxn = find(WT_model.c == 1);
WT_model = setParam(WT_model,'eq',biomass_rxn,0.0289);
KO_model = setParam(KO_model,'eq',biomass_rxn,0.0289);
clear biomass_rxn

% Add additional constraints:
% Add constraints:
AA_constraints = {'EX_gly_e','EX_thr__L_e','EX_arg__L_e','EX_cys__L_e', ...
    'EX_ile__L_e','EX_his__L_e','EX_lys__L_e','EX_leu__L_e'};
WT_model = setParam(WT_model,'lb',AA_constraints,-1);
WT_model = setParam(WT_model,'ub',AA_constraints,0);
KO_model = setParam(KO_model,'lb',AA_constraints,-1);
KO_model = setParam(KO_model,'ub',AA_constraints,0);

% Add additional constraints from flux solution:
extra_constraints = {'EX_fe2_e','EX_pydxn_e','EX_retn_e','EX_ncam_e',... 
    'EX_hxan_e','EX_lnlc_e','EX_fol_e','EX_lnlnca_e','EX_chol_e',...
    'EX_thymd_e'};
WT_model = setParam(WT_model,'lb',extra_constraints,-1e-5);
WT_model = setParam(WT_model,'ub',extra_constraints,0);
KO_model = setParam(KO_model,'lb',extra_constraints,-1e-5);
KO_model = setParam(KO_model,'ub',extra_constraints,0);
% Add inositol constraint:
WT_model = setParam(WT_model,'lb','EX_inost_e',-1e-5);
WT_model = setParam(WT_model,'ub','EX_inost_e',1e-5);
KO_model = setParam(KO_model,'lb','EX_inost_e',-1e-5);
KO_model = setParam(KO_model,'ub','EX_inost_e',1e-5);

% Set constraints for alanine
WT_model = setParam(WT_model,'lb','EX_ala__L_e',0);
WT_model = setParam(WT_model,'ub','EX_ala__L_e',1000);
KO_model = setParam(KO_model,'lb','EX_ala__L_e',0);
KO_model = setParam(KO_model,'ub','EX_ala__L_e',1000);

% Set constraints on release exchange fluxes:
release_constraints = {'EX_lac__D_e','EX_pyr_e','EX_acald_e','EX_3aib__D_e',...
    'EX_acac_e','EX_HC00250_e','EX_2hb_e','EX_3mop_e','EX_thf_e','EX_bhb_e',...
    'EX_HC00229_e','EX_acetone_e','EX_4abut_e','EX_3aib_e','EX_ocdca_e','EX_hdca_e',...
    'EX_ttdca_e','EX_glyc_e','EX_succ_e','EX_mal__L_e','EX_dag_cho_e','EX_sbt__D_e',...
    'EX_ac_e','EX_2mcit_e','EX_cit_e','EX_drib_e','EX_ptrc_e','EX_4mop_e','EX_akg_e',...
    'EX_Rtotal_e','EX_mag_cho_e','EX_tag_cho_e','EX_mercplaccys_e','EX_3mob_e',...
    'EX_34hpp_e','EX_fru_e','EX_5oxpro_e','EX_gluala_e','EX_ddca_e','EX_arach_e',...
    'EX_Rtotal2_e','EX_Rtotal3_e','EX_HC01444_e','EX_ala_B_e','EX_pglyc_cho_e',...
    'EX_crmp_cho_e','EX_ps_cho_e','EX_sph1p_e','EX_slfcys_e','EX_pheacgln_e',...
    'EX_4hphac_e','EX_orn_e','EX_tymsf_e','EX_34dhphe_e','EX_gthrd_e','EX_dopa_e',...
    'EX_mthgxl_e','EX_h2o2_e','EX_hco3_e','EX_nh4_e','EX_so4_e','EX_taur_e',...
    'EX_sl__L_e','EX_urea_e','EX_gua_e','EX_creat_e','EX_urate_e','EX_ade_e',...
    'EX_no_e','EX_orot_e','EX_ura_e','EX_cbasp_e','EX_hista_e','EX_gsn_e','EX_adn_e',...
    'EX_cgly_e','EX_citr__L_e','EX_dgsn_e','EX_35cgmp_e','EX_ahcys_e','EX_dad_2_e',...
    'EX_atp_e','EX_amp_e','EX_camp_e','EX_fad_e','DM_dgtp_m','EX_5mta_e','DM_dgtp_n',...
    'DM_datp_m','DM_datp_n','EX_aicar_e','EX_ins_e','EX_din_e','EX_thym_e','EX_ser__D_e',...
    'EX_cytd_e','EX_gthox_e','EX_dopasf_e','EX_nrpphrsf_e','EX_tdchola_e','EX_tchola_e',...
    'EX_digalsgalside_cho_e'};
WT_model = setParam(WT_model,'lb',release_constraints,0);
WT_model = setParam(WT_model,'ub',release_constraints,1e-5);
KO_model = setParam(KO_model,'lb',release_constraints,0);
KO_model = setParam(KO_model,'ub',release_constraints,1e-5);
clear AA_constraints extra_constraints release_constraints

% Change objective to NGAM reaction 
WT_model = changeObjective(WT_model,'ATPM');
KO_model = changeObjective(KO_model,'ATPM');
% Check max ATP production
sol_WT_ATP = optimizeCbModel(WT_model);
disp("WT ATP production before transport: ")
disp(sol_WT_ATP.f)
sol_KO_ATP = optimizeCbModel(KO_model);
disp("KO ATP production before transport: ")
disp(sol_KO_ATP.f)

% Get exchange rxns:
exch_rxns = getExchangeRxns(WT_model);
exch_idx = find(ismember(WT_model.rxns,exch_rxns));
WT_flux = sol_WT_ATP.x(exch_idx);
WT_lb = WT_model.lb(exch_idx);
WT_ub = WT_model.ub(exch_idx);
KO_flux = sol_KO_ATP.x(exch_idx);
KO_lb = KO_model.lb(exch_idx);
KO_ub = KO_model.ub(exch_idx);
% Create a table with solutions
exch_fluxes_solution = table(exch_rxns, WT_flux, WT_lb, WT_ub,  KO_flux,KO_lb,KO_ub);% Extract relevant data
exch_fluxes_solution = sortrows(exch_fluxes_solution,"WT_flux","descend");

rxns_from_sampling = {'PYRt2p','r2516','NADHtpu','NADtx','Htx','L_LACtm',...
    'r1464','r1290','ILEtmi','PROtm','DIDPtn','ALAtN1','PRODt2r','PRO_Dtde',...
    'DITPtn','SUCCtm','AKGMALtm','r0911','ASPGLUm','ACACt2m'};
additional_rxns = {'r0836','CITRt2m','FUMSO3tm','NACASPtm','r0835','ACt2m',...
    'CITRt2m','r0829','r0915','r0834'};
rxns_from_sampling = [rxns_from_sampling,additional_rxns];
WT_model = setParam(WT_model,'lb',rxns_from_sampling,-1e-5);
WT_model = setParam(WT_model,'ub',rxns_from_sampling,1e-5);
KO_model = setParam(KO_model,'lb',rxns_from_sampling,-1e-5);
KO_model = setParam(KO_model,'ub',rxns_from_sampling,1e-5);

sol_WT_ATP_transp = optimizeCbModel(WT_model);
disp(sol_WT_ATP_transp.f);
sol_KO_ATP_transp = optimizeCbModel(KO_model);
disp(sol_KO_ATP_transp.f);

% Get exchange rxns:
exch_rxns = getExchangeRxns(WT_model);
exch_idx = find(ismember(WT_model.rxns,exch_rxns));
WT_flux = sol_WT_ATP_transp.x(exch_idx);
WT_lb = WT_model.lb(exch_idx);
WT_ub = WT_model.ub(exch_idx);
KO_flux = sol_KO_ATP_transp.x(exch_idx);
KO_lb = KO_model.lb(exch_idx);
KO_ub = KO_model.ub(exch_idx);
% Create a table with solutions
exch_fluxes_solution = table(exch_rxns, WT_flux, WT_lb, WT_ub,  KO_flux,KO_lb,KO_ub);% Extract relevant data
exch_fluxes_solution = sortrows(exch_fluxes_solution,"WT_flux","descend");

if random_sampling == "yes"
    WT_model = setParam(WT_model,'eq','ATPM',sol_WT_ATP_transp.f,2);
    KO_model = setParam(KO_model,'eq','ATPM',sol_KO_ATP_transp.f,2);
    % Run random sampling with 2000 samples
    solutions_WT = randomSampling(WT_model, 2000, true, false, true, [], false);
    solutions_KO = randomSampling(KO_model, 2000, true, false, true, [], false);

    WT_model.id = 'WT_model';
    KO_model.id = 'KO_model';
    out_dir_root = "data/modeling/solution_eval/20250120_rsWithRosemary";

    res = evaluateSampling_updated(WT_model, KO_model, solutions_WT, solutions_KO, 0.01,0,out_dir_root);

    filter_rows = abs(res.mean_sample_WT_model) > 1e-5 & abs(res.mean_sample_WT_model) > 1e-5;
    res = res(filter_rows, :);
    res.mets_in = cellfun(@(x) strjoin(x, '; '), res.mets_in, 'UniformOutput', false);
    res.mets_out = cellfun(@(x) strjoin(x, '; '), res.mets_out, 'UniformOutput', false);
end

% res is saved udner /code/data/modeling/rxn_results_sampling.csv
% KO_model and WT_model saved under /code/data/modeling/KO_model.mat or
% WT_model.mat
% after filter criteria, res is stored as
% /code/data/modeling/20250122_rs_results_after_filter.csv (this is used
% later)

clear exch_rxns exch_idx WT_flux WT_lb WT_ub KO_flux KO_lb KO_ub  



 
