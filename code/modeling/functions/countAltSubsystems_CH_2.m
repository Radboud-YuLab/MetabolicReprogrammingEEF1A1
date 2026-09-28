function [res_subsystems] = countAltSubsystems_CH_2(model_ref,model,sampling_res, out_dir)

% Count unique reactions per model subsystem
subsystem_changed = sampling_res.subsystem;
all_subs = vertcat(model_ref.subSystems(:));
[cnt_unique,~,idx] = unique(all_subs);
n_counts = accumarray(idx(:),1);
temp = table(categorical(cnt_unique),n_counts);
temp.Properties.VariableNames = {'subsytem','counts_total'};
% Count unique subsystems per significant metabolic reactions

%Write subsystems
[C,~,ic] = unique(vertcat( subsystem_changed(:) ));
res_subsystems = table(categorical(C), accumarray(ic,1));
res_subsystems.Properties.VariableNames = {'subsytem','counts_altered'};
%res_subsystems = rmmissing(res_subsystems);
% join overlapping subsystem only
res_subsystems = innerjoin(res_subsystems,temp);
res_subsystems.pct_rxn_change = 100*res_subsystems.counts_altered./res_subsystems.counts_total;
res_subsystems = sortrows(res_subsystems,4,'descend');
out_name = strcat(model_ref.id,"_",model.id,"_SamplingSubsystemsAltered.csv"); % CHANGED BY CYRIEL
writetable(res_subsystems,fullfile(out_dir,out_name),'Delimiter',',');  

[m,~] = size(res_subsystems);
out_name_plot = strcat(model_ref.id,"_",model.id,"_SamplingSubsystemsAltered.jpg");
f = figure('visible','off','WindowState','maximized'); % only if visible is set to 'on' the figure is printed with correct size
set(gca, 'XTickLabelRotation', 40, 'FontSize',12,'FontName','Times');
set(gcf,'position',[0 0 2400 2000]);
grid on;
barh(flip(res_subsystems.pct_rxn_change),0.40);
yticks(1:m);
yticklabels(flip(res_subsystems.subsytem));
xticks([0 20 40 60 80 100])
ylabel('Metabolic subsystem')
xlabel('% Reactions with statistically different sample flux')
legend(model_ref.id, 'Location', 'east','Interpreter','none')
title(strcat('Flux sampling comparision per metabolic subystem (normalized by reaction count) between',{' '},model_ref.id,{' '},'and',{' '},model.id),'Interpreter','none');
print(fullfile(out_dir,out_name_plot),'-djpeg','-r300');
