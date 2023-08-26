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
for iter = 1:4
for Idx = 1:16

tic
iter
Idx
params = zeros(3,22);

% For "fixed" parameters, set sampling range for solver (multipliers)
WideRange = [0.001 1 100];
TightRange = [0.75 1 1.25];

Perc_LW = Conjugation.Clumping.Perc_AA(Idx);
Sample = Conjugation.Clumping.BactYeastPair(Idx);
Exp = categorical("230221_DynMannose");

BactFcn = string(Conjugation.Clumping.BactFcn(Idx));
YeastFcn = string(Conjugation.Clumping.YeastFcn(Idx));

% Set number guesses (simulations) for model to perform
number_guesses = 1E5;


% Parameters...
% 1-6: Rb, Ry, Kb, Ky, Cb, Cy,...
% 7-14: Gb, Gy (global AA needed by b,y), Ab, Ay (secretion from other cell), kb, ky, D, Gam,...
% 17-20: Pb, Py (proximity feeding bumps for clumped yeast), Rlb, Rly (clumped growth rates),...
% 21-22: Rl (rate of clumping per coincident b,y), GamL (TKC rate for clumped cells)

% Most non-clump params carried over from conjugation batch fits
% for c = 1:16
%     params(:,c) = AllFits.Saved_Params.ConjPlotParams(Idx,c) .* TightRange;
% end
% 
% % Fixed global AA
% params(:,7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100) .* [1 1 1];
% params(:,8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100) .* [1 1 1];
% 
% % Expected TKC rate for non-clumped
% params(:,14) = 1E-4 * AllFits.Saved_Params.Fluor_Conversion.Gamma(1) .* TightRange;
% 
% % Proximity boost for amino acid secretion
% params(:,17) = [1 5 10];
% params(:,18) = [1 5 10];
% 
% % Clumped carrying capacities taken from max of clump image data
% params(:,19) = AllFits.Saved_Params.ConjPlotParams(Idx,1) .* TightRange;
% params(:,20) = AllFits.Saved_Params.ConjPlotParams(Idx,2) .* TightRange;
% 
% % Clumping rate starting blind
% params(:,21) = [.0001 .05 .01];
% 
% % Clumped TKC rate starting w/ literature value
% params(:,22) = 4E-3 * AllFits.Saved_Params.Fluor_Conversion.Gamma(1) .* [.001 1 2];
% 
% params(:,1) = params(:,1) ./ 2;


% OR use some saved version for params
params(2,:) = AllFits.Saved_Params.ClumpPlotParams(Idx,:);
% params(2,:) = compile(Idx,:);% + [-.2 0 0 800 -0.9 0.15 0 0 0 0 0 0 0 0 0 0];
% params(2,:) = Conjugation.Clumping.Mean_Params{Idx};
if BactFcn=="wtB"
    params(2,11) = 0;
end
if YeastFcn=="wtY"
    params(2,12) = 0;
end

params(1,:) = params(2,:) .* .75;
params(3,:) = params(2,:) .* 1.25;

params(:,7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100) .* [1 1 1];
params(:,8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100) .* [1 1 1];


param_bounds_lower = params(1,:);
param_bounds_upper = params(3,:);

t_nixOne = find(Conjugation.Clumping.Spark_Time{Idx}.FitDay~=1);


% Set init conditions of B, Y, and T [B0 Y0 T0] based on experimental setups
Binit = nanmean(Conjugation.Clumping.Reads_Spark{Idx}(t_nixOne(2:4),1,:),'all');
Yinit = nanmean(Conjugation.Clumping.Reads_Spark{Idx}(t_nixOne(2:4),2,:),'all');
Tinit = 0;
Linit = 0;
Lbinit = 1E-10;
Lyinit = 1E-10;
initial_conditions = [Binit Yinit Tinit Linit Lbinit Lyinit]; 

% Set up hourly time vector over 6 days
t_spark = Conjugation.Clumping.Spark_Time{Idx}; % For batch solver, only give hours between dilutions
dilution_days = unique(t_spark.FitDay)'; % For batch solver, enter todal days diluted
num_replicates = size(Conjugation.Clumping.Reads_Spark{Idx},3);

dilution_days(1) = [];


% Sample latin hypercube: first line samples linearly
% params_LHS = lhsu(param_bounds_lower,param_bounds_upper,number_guesses);

% Use ALL THREE of these lines to sample better in log space
params_LHS = lhsu(log(param_bounds_lower),log(param_bounds_upper),number_guesses);
params_LHS = exp(params_LHS);
params_LHS(isnan(params_LHS)) = 0; %If any of your params were 0, you'll get nans that need to be corrected



%Need to txform data to put reads only on 24-hour marks
fluorescence_measured = (Conjugation.Clumping.Reads_Spark{Idx,1}(t_nixOne,[1:2 4:7],:));

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
        
        try
        [~,y_add] = ode45(@(t,y) Clump_ODEs(t,y,params_temp),t_measured,New_conditions);
        y_batch = vertcat(y_batch,y_add);
        catch
        [~,y_add] = ode45(@(t,y) Clump_ODEs(t,y,(params_temp .* .75)),t_measured,New_conditions);
        
        % This block accounts for "breakaway" bact from clumps, assuming a
        % max bact/clump of 221 (converts to 4562 citrine). Any clumped
        % bact over this value get redefined as free bact 
        Lb_over = find(y_add(:,4)>4562);
        y_add(Lb_over,1) = y_add(Lb_over,1) + (y_add(Lb_over,4) - 4562);
        y_add(Lb_over,4) = 4562;

        y_batch = vertcat(y_batch,y_add);
 
        % Add noise into dilution step to account for bottleneck variation
        noise_min = 0.950;
        noise_max = 1.05;
        noise = noise_min + (noise_max - noise_min).*rand(size(y_batch(end,:,:)));
        New_conditions = noise .* y_batch(end,:,:) ./ 10;
%         New_conditions = y_batch(end,:,:) ./ 10;
        end
    end
    fluorescence_model(:,:,1,guess) = y_batch;

end

fluorescence_model(:,3,:,:) = fluorescence_model(:,3,:,:) ./ (2*AllFits.Saved_Params.Fluor_Conversion.Citrine);      
    
% Columns 1 & 2 will be FREE bact and yeast, col.s 5 & 6 CLUMPED B & Y, but
% measurements for B & Y are for TOTAL counts, so need to combine free and
% clumped versions into one vector each.
fluorescence_model(:,1,:,:) = fluorescence_model(:,1,:,:) + fluorescence_model(:,5,:,:);
fluorescence_model(:,2,:,:) = fluorescence_model(:,2,:,:) + fluorescence_model(:,6,:,:);



% % Need to reset TKC values to counts to compare with data
% fluorescence_model(:,3,:,:) = fluorescence_model(:,3,:,:) ./ AllFits.Saved_Params.Fluor_Conversion.Gamma(1);

toc

% Manually calculate least-squares error and rank parameter sets

% Reshape measurements and model predictions for easy matrix operations
fluorescence_model_4D = repmat(fluorescence_model(:,:,1,:),[1 1 size(fluorescence_measured,3) 1]);
fluorescence_measured_4D = repmat(fluorescence_measured,[1 1 1 number_guesses]);

% Calculate RSS error between measurements and model predictions
RSS_error = (fluorescence_measured_4D(:,1:4,:,:) - fluorescence_model_4D(:,1:4,:,:)).^2; 
RSS_error = sum(RSS_error,1:3); 
RSS_error = squeeze(RSS_error);

% Rank parameters and model predictions by RSS error
[RSS_error, RSS_sort_idx] = sort(RSS_error,'ascend');
params_LHS_ranked = params_LHS(RSS_sort_idx,:);
fluorescence_model_ranked = fluorescence_model(:,:,RSS_sort_idx);

% Plot measuremetns and top model prediction
number_top_fits_plot = 5;

% Transform data for plot into what is actually plotted in other analyses
B_meas = squeeze(fluorescence_measured(:,1,:)); 
Y_meas = squeeze(fluorescence_measured(:,2,:)); 
T_meas = squeeze(fluorescence_measured(:,3,:)); 
L_meas = squeeze(fluorescence_measured(:,4,:)); 
Lb_meas = squeeze(fluorescence_measured(:,5,:)); 
Ly_meas = squeeze(fluorescence_measured(:,6,:)); 

B_mod = squeeze(fluorescence_model_ranked(:,1,1:number_top_fits_plot)); % ./ AllFits.Saved_Params.K.Bact_M2('Cells_200uL');
Y_mod = squeeze(fluorescence_model_ranked(:,2,1:number_top_fits_plot)); % ./ AllFits.Saved_Params.K.Yeast('Cells_200uL');
T_mod = squeeze(fluorescence_model_ranked(:,3,1:number_top_fits_plot));
L_mod = squeeze(fluorescence_model_ranked(:,4,1:number_top_fits_plot));
Lb_mod = squeeze(fluorescence_model_ranked(:,5,1:number_top_fits_plot));
Ly_mod = squeeze(fluorescence_model_ranked(:,6,1:number_top_fits_plot));

t_tot = t_spark.CumTime(t_nixOne);

figure(1)
subplot(4,4,Idx)
hold on
plot(t_tot,B_meas,'r.'); % Plot bacteria measurements
plot(t_tot,Y_meas,'b.'); % Plot yeast measurements
plot(t_tot,B_mod,'r'); % Plot bacteria model
plot(t_tot,Y_mod,'b'); % Plot yeast model
ylim([0 7000])
text = string(Sample) + ' at ' + num2str(Perc_LW) + '% LW';
title(text,'Interpreter','None');
hold off
sgtitle('Total bact and yeast');

figure(2)
subplot(4,4,Idx)
hold on
plot(t_tot,T_meas,'g.'); % Plot yeast measurements
plot(t_tot,T_mod,'g'); % Plot yeast model
ylim([0 80]);
text = string(Sample) + ' at ' + num2str(Perc_LW) + '% LW';
title(text,'Interpreter','None');
hold off
sgtitle('Total TKC');

figure(3)
subplot(4,4,Idx)
hold on
plot(t_tot,Lb_meas,'color',[1 .3 0],'marker','.','linestyle','None'); % Plot bacteria measurements
plot(t_tot,Ly_meas,'color',[0 .5 .8],'marker','.','linestyle','None'); % Plot yeast measurements
plot(t_tot,Lb_mod,'color',[1 .3 0],'linestyle','-'); % Plot bacteria measurements
plot(t_tot,Ly_mod,'color',[0 .5 .8],'linestyle','-'); % Plot yeast model
ylim([0 6000])
text = string(Sample) + '  at ' + num2str(Perc_LW) + '% LW';
title(text,'Interpreter','None');
hold off
sgtitle('Clumped B and Y');

figure(5)
subplot(4,4,Idx)
hold on
plot(t_tot,L_meas,'k.'); % Plot bacteria measurements
plot(t_tot,L_mod,'k'); % Plot bacteria model
ylim([0 800])
text = string(Sample) + '  at ' + num2str(Perc_LW) + '% LW';
title(text,'Interpreter','None');
hold off
sgtitle('Number of clumps');


Conjugation.Clumping.LHS_Params{Idx} = params_LHS_ranked;

 
    
% end


% Calculate stats on top hits

num = 1000; %Number of top hits to calculate from
Pairs = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'};

for p = 1:4
    pair = Pairs{p};
    catAll = [];

    pIdx = find(Conjugation.Clumping.BactYeastPair==pair);
    
    for r = 1:numel(pIdx)
        row = pIdx(r);
        catAll = vertcat(catAll,Conjugation.Clumping.LHS_Params{row}(1:num,:));
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

    for r = 1:numel(pIdx)
        row = pIdx(r);
        Conjugation.Clumping.Mean_Params{row} = avgs;
        Conjugation.Clumping.Median_Params{row} = meds;
        Conjugation.Clumping.StdDev_Params{row} = stdev;
        Conjugation.Clumping.SEM_Params{row} = SEM;
        Conjugation.Clumping.Tscore_Params{row} = Tscore;
        for i = 1:14
        temp(1,i) = CI95{i}(1);
        temp(2,i) = CI95{i}(2);
        end
        Conjugation.Clumping.CI95_Params{row} = temp;
        
        Perc_LW = Conjugation.Clumping.Perc_AA(row);
        Conjugation.Clumping.Mean_Params{row}(7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100);
        Conjugation.Clumping.Mean_Params{row}(8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100);
        Conjugation.Clumping.Median_Params{row}(7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100);
        Conjugation.Clumping.Median_Params{row}(8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100);
    end
end
 

end
end