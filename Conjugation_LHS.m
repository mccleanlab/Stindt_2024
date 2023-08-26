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
for iter = 1
for z = 1:numel(rows)
Idx = rows(z);

tic
iter
Idx
% For "fixed" parameters, set sampling range for solver (multipliers)
WideRange = [0.001 1 100];
TightRange = [0.75 1 1.25];
Perc_LW = Conjugation.ByExperiment.Perc_AA(Idx);
Sample = Conjugation.ByExperiment.BactYeastPair(Idx);
Exp = categorical("230221_DynMannose");

BactFcn = string(Conjugation.ByExperiment.BactFcn(Idx));
YeastFcn = string(Conjugation.ByExperiment.YeastFcn(Idx));

% Set number guesses (simulations) for model to perform
number_guesses = 1E5;


% Idx = find(Conjugation.ByExperiment.Perc_AA==Perc_LW & Conjugation.ByExperiment.BactYeastPair==Sample &...
%     Conjugation.ByExperiment.Experiment==Exp);

% % Parameters time bounds: Rb, Ry, Kb, Ky, Cb, Cy, Gb (global AA needed by Bact),
% % Gy, Ab, Ay, kb, ky, D, Gam. Middle row is actually fixed
% params(:,1) = AllFits.Saved_Params.R.Bact_TR6('Fluor/hr') .* TightRange;
% params(:,2) = AllFits.Saved_Params.R.Yeast('Fluor/hr') .* TightRange;
% params(:,3) = AllFits.Saved_Params.K.Bact_TR6('Fluor') .* TightRange;
% params(:,4) = AllFits.Saved_Params.K.Yeast('Fluor') .* TightRange;
% params(:,5) = AllFits.Saved_Params.C.Cb .* [0.1 1 10];
% params(:,6) = AllFits.Saved_Params.C.Cy .* [0.1 1 10];
% params(:,7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100) .* [1 1 1];
% params(:,8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100) .* [1 1 1];
% % A subscripts define WHAT CONSUMES the AA, but are saved based on WHAT
% % SECRETES that AA, because there are differences b/w cross & wt rates
% params(:,9) = AllFits.Saved_Params.A.(YeastFcn)('Molar/Fluor*hr') .* WideRange .* 10;
% params(:,10) = AllFits.Saved_Params.A.(BactFcn)('Molar/Fluor*hr') .* WideRange .* 10;
% params(:,11) = AllFits.Saved_Params.k.(BactFcn)('Molar') .* TightRange;
% params(:,12) = AllFits.Saved_Params.k.(YeastFcn)('Molar') .* TightRange;
% params(:,13) = AllFits.Saved_Params.D.L .* TightRange;
% 
% % Set a wider range for Gamma since that's what we're checking out
% % GamRange = [0.001  1 10];
% params(:,14) = 4E-3 * AllFits.Saved_Params.Fluor_Conversion.Gamma(1);

% params(:,15) = [0.1 0.2 0.8];
% params(:,16) = [0.1 0.2 0.8];

% OR use some saved version for params
% clear params
params(2,:) = Conjugation.ByExperiment.Mean_Params_wDeath{Idx};% + [-.2 0 0 800 -0.9 0.15 0 0 0 0 0 0 0 0 0 0];
% Fits = table2array(BatchFits);
% params(2,1:6) = Fits(6,1:6);
% params(2,15:16) = Fits(6,7:8);

% params(2,1) = .7;
% params(2,5) = 1.6;
% params(2,6) = .9;
% params(2,3) = 4E3;
% params(2,4) = 2E3;

params(1,:) = params(2,:) .* 0.5;
params(3,:) = params(2,:) .* 1.5;

% params(:,5) = .5 .* [.5 1 2.2];
% params(:,6) = .5 .* [.5 1 2.2];

params(3,1:6) = [.85 .85 4000 3500 1.5 1.5];
params(3,15:16) = .85;

params(1,1:6) = [.5 .4 3000 2500 .4 .4];
params(1,15:16) = .1;

params(:,7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100) .* [1 1 1];
params(:,8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100) .* [1 1 1];
% params(:,9) = AllFits.Saved_Params.A.(YeastFcn)('Molar/Fluor*hr') .* WideRange .* 10;
% params(:,10) = AllFits.Saved_Params.A.(BactFcn)('Molar/Fluor*hr') .* WideRange .* 10;

if BactFcn=="wtB"
    params(2,11) = 0;
% elseif BactFcn=="crossB"
%     params(3,9) = 1E-5;
%     params(1,9) = 1E-10;
end
if YeastFcn=="wtY"
    params(2,12) = 0;
% elseif YeastFcn=="crossY"
%     params(3,12) = 1E-4;
%     params(3,12) = 1E-10;
end


% params(:,11) = AllFits.Saved_Params.k.(BactFcn)('Molar') .* TightRange;
% params(:,12) = AllFits.Saved_Params.k.(YeastFcn)('Molar') .* TightRange;
% params(:,13) = AllFits.Saved_Params.D.L .* TightRange;
params(:,14) = .004 * AllFits.Saved_Params.Fluor_Conversion.Gamma(1) .* [.01 1 2.5];

% OPTIONAL: add death terms for bact & yeast in positions 15 & 16 (must
% change solver below to be "...ODEs_wDeath"
% params(:,15) = [0.1 0.2 0.8];
% params(:,16) = [0.1 0.2 0.8];


param_bounds_lower = params(1,:);
param_bounds_upper = params(3,:);

t_nixOne = find(Conjugation.ByExperiment.Spark_Time{Idx}.FitDay~=1);


% Set init conditions of B, Y, and T [B0 Y0 T0] based on experimental setups
Binit = nanmean(Conjugation.ByExperiment.Reads_Spark{Idx}(t_nixOne(2:4),1,:),'all');
Yinit = nanmean(Conjugation.ByExperiment.Reads_Spark{Idx}(t_nixOne(2:4),2,:),'all');
Tinit = 0;
initial_conditions = [Binit Yinit Tinit]; 

% Set up hourly time vector over 6 days
t_spark = Conjugation.ByExperiment.Spark_Time{Idx}; % For batch solver, only give hours between dilutions
dilution_days = unique(t_spark.FitDay)'; % For batch solver, enter todal days diluted
num_replicates = size(Conjugation.ByExperiment.Reads_Spark{Idx},3);

dilution_days(1) = [];


% Sample latin hypercube: first line samples linearly
params_LHS = lhsu(param_bounds_lower,param_bounds_upper,number_guesses);

% Use ALL THREE of these lines to sample better in log space
% params_LHS = lhsu(log(param_bounds_lower),log(param_bounds_upper),number_guesses);
% params_LHS = exp(params_LHS);
% params_LHS(isnan(params_LHS)) = 0; %If any of your params were 0, you'll get nans that need to be corrected



%Need to txform data to put reads only on 24-hour marks
fluorescence_measured = (Conjugation.ByExperiment.Reads_Spark{Idx,1}(t_nixOne,[1 2 4],:));

% Run model w/ each parameter set
fluorescence_model = zeros([size(fluorescence_measured,1) size(fluorescence_measured,2) num_replicates number_guesses]);


% Batch Solver for SPARK DATA
parfor guess = 1:number_guesses
    params_temp = params_LHS(guess,:);
    y_batch = [];
    
    for d = dilution_days
        t_measured = t_spark.FitTime(t_spark.FitDay==d,1);
        if d==2
            New_conditions = initial_conditions;
        end
        
        [~,y_add] = ode45(@(t,y) Conjugation_ODEs_wDeath(t,y,params_temp),t_measured,New_conditions);
        y_batch = vertcat(y_batch,y_add);
        
        noise_min = 0.80;
        noise_max = 1.20;
        noise = noise_min + (noise_max - noise_min).*rand(size(y_batch(end,:,:)));
        New_conditions = noise .* y_batch(end,:,:) ./ 10;
%         New_conditions = y_batch(end,:,:) ./ 10;
    end
    fluorescence_model(:,:,guess) = y_batch;

end

toc
% Manually calculate least-squares error and rank parameter sets

% Reshape measurements and model predictions for easy matrix operations
% fluorescence_model_4D = repmat(fluorescence_model,[1 1 size(fluorescence_measured,3),1]);
fluorescence_measured_4D = repmat(fluorescence_measured,[1 1 1 number_guesses]);

% Calculate RSS error between measurements and model predictions
RSS_error = (fluorescence_measured_4D(:,:,:,:) - fluorescence_model(:,:,:,:)).^2; %Only looking at error between Bact & Yeast growth, NOT TKC
RSS_error = sum(RSS_error,1:2); 
RSS_error = squeeze(RSS_error);

% Rank parameters and model predictions by RSS error
[RSS_error, RSS_sort_idx] = sort(RSS_error,'ascend');
params_LHS_ranked = params_LHS(RSS_sort_idx,:);
fluorescence_model_ranked = fluorescence_model(:,:,RSS_sort_idx);

% Plot measuremetns and top model prediction
number_top_fits_plot = 5;

% Transform data for plot into what is actually plotted in other analyses
B_meas = squeeze(fluorescence_measured(:,1,:)); % ./ AllFits.Saved_Params.K.Bact_M2('Cells_200uL');
Y_meas = squeeze(fluorescence_measured(:,2,:)); % ./ AllFits.Saved_Params.K.Yeast('Cells_200uL');
T_meas = squeeze(fluorescence_measured(:,3,:)); 

B_mod = squeeze(fluorescence_model_ranked(:,1,1:number_top_fits_plot)); % ./ AllFits.Saved_Params.K.Bact_M2('Cells_200uL');
Y_mod = squeeze(fluorescence_model_ranked(:,2,1:number_top_fits_plot)); % ./ AllFits.Saved_Params.K.Yeast('Cells_200uL');
T_mod = squeeze(fluorescence_model_ranked(:,3,1:number_top_fits_plot));

t_tot = t_spark.CumTime(t_nixOne);

i = Idx-min(rows)+1;

subplot(4,4,i)
hold on
plot(t_tot,B_meas,'r.'); % Plot bacteria measurements
plot(t_tot,Y_meas,'b.'); % Plot yeast measurements
plot(t_tot,B_mod,'r'); % Plot bacteria model
plot(t_tot,Y_mod,'b'); % Plot yeast model
ylim([0 6000])
text = 'Row ' + string(num2str(Idx)) + ', ' + string(Sample) + ' at ' + num2str(Perc_LW) + '% LW';
title(text,'Interpreter','None');
hold off

% figure 2
% subplot(4,4,i)
% hold on
% plot(t_tot,T_meas,'g.'); % Plot yeast measurements
% plot(t_tot,T_mod,'g'); % Plot yeast model
% % title('Predicted Total Transconjugants, ?=0.004');
% ylim([0 100]);
% title(['Top 5 TKC, ' Sample ' at ' num2str(Perc_LW) '% LW'],'interpreter','none')
% hold off

Conjugation.ByExperiment.LHS_Params_wDeath{Idx} = params_LHS_ranked;

end 
    

% Calculate stats on top hits: death rate included

num = 1000; %Number of top hits to calculate from
Pairs = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'};

for p = 1:4
    pair = Pairs{p};
    catAll = [];

    Idx = find(Conjugation.ByExperiment.BactYeastPair==pair &...
        Conjugation.ByExperiment.Experiment=='230221_DynMannose_Plus');
    
    for r = 1:numel(Idx)
        row = Idx(r);
        catAll = vertcat(catAll,Conjugation.ByExperiment.LHS_Params_wDeath{row}(1:num,:));
    end

    for i = 1:size(catAll,2)

        avgs(i) = nanmean(catAll(:,i));
        meds(i) = nanmedian(catAll(:,i));

        stdev(i) = std(catAll(:,i));
    %     var(i) = var(catAll(:,i));

        SEM(i) = stdev(i)/sqrt(numel(catAll(:,i)));
        Tscore{i} = tinv([0.05 0.95],numel(catAll(:,i))-1);

        CI95{i} = avgs(i) + Tscore{i}.*SEM(i);
        
    end

    for r = 1:numel(Idx)
        row = Idx(r);
        Conjugation.ByExperiment.Mean_Params_wDeath{row} = avgs;
        Conjugation.ByExperiment.Median_Params_wDeath{row} = meds;
        Conjugation.ByExperiment.StdDev_Params_wDeath{row} = stdev;
    %     Conjugation.ByExperiment.Var_Params{33:64} = var;
        Conjugation.ByExperiment.SEM_Params_wDeath{row} = SEM;
        Conjugation.ByExperiment.Tscore_Params_wDeath{row} = Tscore;
        for i = 1:14
        temp(1,i) = CI95{i}(1);
        temp(2,i) = CI95{i}(2);
        end
        Conjugation.ByExperiment.CI95_Params_wDeath{row} = temp;
        
        Perc_LW = Conjugation.ByExperiment.Perc_AA(row);
        Conjugation.ByExperiment.Mean_Params_wDeath{row}(7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100);
        Conjugation.ByExperiment.Mean_Params_wDeath{row}(8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100);
        Conjugation.ByExperiment.Median_Params_wDeath{row}(7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100);
        Conjugation.ByExperiment.Median_Params_wDeath{row}(8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100);
    end
end
end 

%% Calculate stats on top hits

num = 1000; %Number of top hits to calculate from
Pairs = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'};

for p = 1:4
    pair = Pairs{p};
    catAll = [];

    Idx = find(Conjugation.ByExperiment.BactYeastPair==pair &...
        Conjugation.ByExperiment.Experiment=='230221_DynMannose_Plus');
    
    for r = 1:numel(Idx)
        row = Idx(r);
        catAll = vertcat(catAll,Conjugation.ByExperiment.LHS_Params{row}(1:num,:));
    end

    for i = 1:size(catAll,2)

        avgs(i) = nanmean(catAll(:,i));
        meds(i) = nanmedian(catAll(:,i));

        stdev(i) = std(catAll(:,i));
    %     var(i) = var(catAll(:,i));

        SEM(i) = stdev(i)/sqrt(numel(catAll(:,i)));
        Tscore{i} = tinv([0.05 0.95],numel(catAll(:,i))-1);

        CI95{i} = avgs(i) + Tscore{i}.*SEM(i);
        
    end

    for r = 1:numel(Idx)
        row = Idx(r);
        Conjugation.ByExperiment.Mean_Params{row} = avgs;
        Conjugation.ByExperiment.Median_Params{row} = meds;
        Conjugation.ByExperiment.StdDev_Params{row} = stdev;
    %     Conjugation.ByExperiment.Var_Params{33:64} = var;
        Conjugation.ByExperiment.SEM_Params{row} = SEM;
        Conjugation.ByExperiment.Tscore_Params{row} = Tscore;
        for i = 1:14
        temp(1,i) = CI95{i}(1);
        temp(2,i) = CI95{i}(2);
        end
        Conjugation.ByExperiment.CI95_Params{row} = temp;
        
        Perc_LW = Conjugation.ByExperiment.Perc_AA(row);
        Conjugation.ByExperiment.Mean_Params{row}(7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100);
        Conjugation.ByExperiment.Mean_Params{row}(8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100);
        Conjugation.ByExperiment.Median_Params{row}(7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100);
        Conjugation.ByExperiment.Median_Params{row}(8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100);
    end
end

% end

