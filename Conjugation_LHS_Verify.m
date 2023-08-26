% This script runs IN TANDEM WITH Conjugation_DataPrep. If you don't have
% the output variables of that script (i.e. the data), this won't work.


%% Fitting Gam Equations

% Notes copied from Mathematica script:
% This model includes formation- and growth of transconjugants, T. The
% formation of these is based on reaction-diffusion coincidence of bacteria
% B and yeast Y, with a constant transfer rate, \[Gamma], that is thus only
% dependent on the concentrations of B and Y. T grows with the same
% parameters as Y (i.e. there's assumed no fitness cost). Because of the
% inclusion of T, all terms that are Y-dependent are converted to (Y +
% T)-dependence, e.g. both strains produce limited amino acid and are
% summed in production term; the strains reach carrying capacity together.
% Parameter value for \[Gamma] of 4E-3 /hr based on Volkova et. al.
rows = 65:80;

for z = 1:numel(rows)
Idx = rows(z);

number_guesses = 1;


% For "fixed" parameters, set sampling range for solver (multipliers)
Perc_LW = Conjugation.ByExperiment.Perc_AA(Idx);
% Sample = Conjugation.ByExperiment.BactYeastPair(Idx);
Sample = Conjugation.ByExperiment.AltNom(Idx);

% Parameters time bounds: Rb, Ry, Kb, Ky, Cb, Cy, Gb (global AA needed by Bact),
% Gy, Ab, Ay, kb, ky, D, Gam. Middle row is actually fixed

% param_bounds_lower = AllFits.Saved_Params.ConjPlotParams(z,:);
param_bounds_lower = Conjugation.ByExperiment.Mean_Params_wDeath{Idx};

compile(z,:) = Conjugation.ByExperiment.Mean_Params_wDeath{Idx};
param_bounds_lower(:,7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100);
param_bounds_lower(:,8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100);

param_bounds_upper = param_bounds_lower;

t_nixOne = find(Conjugation.ByExperiment.Spark_Time{Idx}.FitDay~=1);


% Set init conditions of B, Y, and T [B0 Y0 T0] based on experimental setups
Binit = nanmean(Conjugation.ByExperiment.Reads_Spark{Idx}(t_nixOne(3:5),1,:),'all');
Yinit = nanmean(Conjugation.ByExperiment.Reads_Spark{Idx}(t_nixOne(3:5),2,:),'all');
Tinit = 0;
initial_conditions = [Binit Yinit Tinit]; 

% Set up hourly time vector over 6 days
t_spark = Conjugation.ByExperiment.Spark_Time{Idx}; % For batch solver, only give hours between dilutions
dilution_days = unique(t_spark.FitDay)'; % For batch solver, enter todal days diluted
num_replicates = size(Conjugation.ByExperiment.Reads_Spark{Idx},3);

dilution_days(1) = [];


% Sample latin hypercube: first line samples linearly
params_LHS = lhsu(param_bounds_lower,param_bounds_upper,number_guesses);

%Need to txform data to put reads only on 24-hour marks
fluorescence_measured = (Conjugation.ByExperiment.Reads_Spark{Idx,1}(t_nixOne,[1 2 4],:));

% Run model w/ each parameter set
fluorescence_model = zeros([size(fluorescence_measured,1) size(fluorescence_measured,2) num_replicates number_guesses]);


% Batch Solver for SPARK DATA
    guess = 1;
    params_temp = params_LHS(guess,:);
    y_batch = [];
    
    for d = dilution_days
        t_measured = t_spark.FitTime(t_spark.FitDay==d,1);
        if d==2
            New_conditions = initial_conditions;
        end
        
        [~,y_add] = ode45(@(t,y) Conjugation_ODEs_wDeath(t,y,params_temp),t_measured,New_conditions);
        y_batch = vertcat(y_batch,y_add);
        
        New_conditions = y_batch(end,:,:) ./ 10;
    end
    fluorescence_model(:,:,guess) = y_batch;


% Manually calculate least-squares error and rank parameter sets

% Reshape measurements and model predictions for easy matrix operations
fluorescence_measured_4D = repmat(fluorescence_measured,[1 1 1 number_guesses]);

% Calculate RSS error between measurements and model predictions
RSS_error = (fluorescence_measured_4D(:,1:2,:,:) - fluorescence_model(:,1:2,:,:)).^2; %Only looking at error between Bact & Yeast growth, NOT TKC
RSS_error = sum(RSS_error,1:3); 
RSS_error = squeeze(RSS_error);

% Rank parameters and model predictions by RSS error
[RSS_error, RSS_sort_idx] = sort(RSS_error,'ascend');
params_LHS_ranked = params_LHS(RSS_sort_idx,:);
fluorescence_model_ranked = fluorescence_model(:,:,RSS_sort_idx);

% Plot measuremetns and top model prediction
number_top_fits_plot = 1;

% Transform data for plot into what is actually plotted in other analyses
B_meas = squeeze(fluorescence_measured(:,1,:)); % ./ AllFits.Saved_Params.K.Bact_M2('Cells_200uL');
Y_meas = squeeze(fluorescence_measured(:,2,:)); % ./ AllFits.Saved_Params.K.Yeast('Cells_200uL');
T_meas = squeeze(fluorescence_measured(:,3,:)); 

B_mod = fluorescence_model_ranked(:,1,1); % ./ AllFits.Saved_Params.K.Bact_M2('Cells_200uL');
Y_mod = fluorescence_model_ranked(:,2,1); % ./ AllFits.Saved_Params.K.Yeast('Cells_200uL');
T_mod = fluorescence_model_ranked(:,3,1) .* AllFits.Saved_Params.Fluor_Conversion.Gamma(1);

t_tot = t_spark.CumTime(t_nixOne);

% i = Idx-min(rows)+1;

subplot(4,4,z)
hold on
plot(t_tot,B_meas,'r.'); % Plot bacteria measurements
plot(t_tot,Y_meas,'b.'); % Plot yeast measurements
plot(t_tot,B_mod,'r'); % Plot bacteria model
plot(t_tot,Y_mod,'b'); % Plot yeast model
ylim([0 6000])
xlabel('t (hr)');
ylabel('Normalized FU');
text = string(Sample) + ' at ' + num2str(Perc_LW) + '% LW';
% title(text,'Interpreter','None');
title(text);
hold off

% subplot(4,4,i)
% hold on
% plot(t_tot,T_meas,'g.'); % Plot yeast measurements
% plot(t_tot,T_mod,'g'); % Plot yeast model
% title('Predicted Total Transconjugants, ?=0.004');
% % ylim([0 150]);
% text = string(Sample) + ' at ' + num2str(Perc_LW) + '% LW';
% title(text,'interpreter','none')
% hold off

% sgtitle('Best fits for growth parameters with mannose added (no clumping)')
end 
    
% end