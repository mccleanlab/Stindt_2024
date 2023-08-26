%% Test param fits

CherPerCit = AllFits.Saved_Params.Fluor_Conversion.Cherry / AllFits.Saved_Params.Fluor_Conversion.Citrine;

for Idx = 1:16

number_guesses = 1;


% For "fixed" parameters, set sampling range for solver (multipliers)
Perc_LW = Conjugation.Clumping.Perc_AA(Idx);
% Sample = Conjugation.Clumping.BactYeastPair(Idx);
Sample = Conjugation.Clumping.AltNom(Idx);

CitPerCher = AllFits.Saved_Params.Fluor_Conversion.Citrine / AllFits.Saved_Params.Fluor_Conversion.Cherry;

% Parameters...
% 1-6: Rb, Ry, Kb, Ky, Cb, Cy,...
% 7-14: Gb, Gy (global AA needed by b,y), Ab, Ay (secretion from other cell), kb, ky, D, Gam,...
% 17-20: Pb, Py (proximity feeding bumps for clumped yeast), Rlb, Rly (clumped growth rates),...
% 21-22: Rl (rate of clumping per coincident b,y), GamL (TKC rate for clumped cells)

% param_bounds_lower = temp(Idx,:);
% param_bounds_lower = Conjugation.Clumping.Mean_Params{Idx};
param_bounds_lower = AllFits.Saved_Params.ClumpFits_Revised(Idx,:);

% param_bounds_lower(:,7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100);
% param_bounds_lower(:,8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100);
% 
% compile(Idx,:) = param_bounds_lower;


param_bounds_upper = param_bounds_lower;

t_nixOne = find(Conjugation.Clumping.Spark_Time{Idx}.FitDay~=1);


% Set init conditions of B, Y, and T [B0 Y0 T0] based on experimental setups
Binit = nanmean(Conjugation.Clumping.Reads_Spark{Idx}(t_nixOne(2:4),1,:),'all');
Yinit = nanmean(Conjugation.Clumping.Reads_Spark{Idx}(t_nixOne(2:4),2,:),'all');
Tinit = 0;
Linit = 0;
Lbinit = 1E-20;
Lyinit = 1E-20;
initial_conditions = [Binit Yinit Tinit Linit Lbinit Lyinit]; 

% Set up hourly time vector over 6 days
t_spark = Conjugation.Clumping.Spark_Time{Idx}; % For batch solver, only give hours between dilutions
dilution_days = unique(t_spark.FitDay)'; % For batch solver, enter todal days diluted
num_replicates = size(Conjugation.Clumping.Reads_Spark{Idx},3);

dilution_days(1) = [];


% Sample latin hypercube: first line samples linearly
params_LHS = lhsu(param_bounds_lower,param_bounds_upper,number_guesses);

%Need to txform data to put reads only on 24-hour marks
fluorescence_measured = (Conjugation.Clumping.Reads_Spark{Idx,1}(t_nixOne,[1:2 4:7],:));

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
        
        [~,y_add] = ode45(@(t,y) Clump_ODEs(t,y,params_temp),t_measured,New_conditions);

                
        % This block accounts for "breakaway" bact from clumps, assuming a
        % max bact/clump of 221 (converts to 4562 citrine). Any clumped
        % bact over this value get redefined as free bact 
        Lb_max = y_add;
        Lb_max(:,7) = y_add(:,5) - (10 * y_add(:,4) * CherPerCit);
        Lb_over = find(Lb_max(:,7) > 0);
        y_add(Lb_over,1) = y_add(Lb_over,1) + (Lb_max(Lb_over,7));
        y_add(Lb_over,5) = y_add(Lb_over,5) - (Lb_max(Lb_over,7));
        
        y_batch = vertcat(y_batch,y_add);
        
        % Add noise into dilution step to account for bottleneck variation
        noise_min = 0.950;
        noise_max = 1.05;
        noise = noise_min + (noise_max - noise_min).*rand(size(y_batch(end,:,:)));
        New_conditions = noise .* y_batch(end,:,:) ./ 10;
%         New_conditions = y_batch(end,:,:) ./ 10;
    end
    fluorescence_model(:,:,1,guess) = y_batch;
       
fluorescence_model(:,3,:,:) = fluorescence_model(:,3,:,:) ./ (2*AllFits.Saved_Params.Fluor_Conversion.Citrine);    
    
% Columns 1 & 2 will be FREE bact and yeast, col.s 5 & 6 CLUMPED B & Y, but
% measurements for B & Y are for TOTAL counts, so need to combine free and
% clumped versions into one vector each.
fluorescence_model(:,1,:,:) = fluorescence_model(:,1,:,:) + fluorescence_model(:,5,:,:);
fluorescence_model(:,2,:,:) = fluorescence_model(:,2,:,:) + fluorescence_model(:,6,:,:);


% Need to reset TKC values to counts to compare with data
% fluorescence_model(:,3,:,:) = fluorescence_model(:,3,:,:) ./ AllFits.Saved_Params.Fluor_Conversion.Gamma(1);


% Manually calculate least-squares error and rank parameter sets

% Reshape measurements and model predictions for easy matrix operations
fluorescence_model_4D = repmat(fluorescence_model(:,:,1,:),[1 1 size(fluorescence_measured,3) 1]);
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
ylim([0 10000])
xlabel('t (hr)');
ylabel('Normalized FU');
text = string(Sample) + ' at ' + num2str(Perc_LW) + '% LW';
% title(text,'Interpreter','None');
title(text);
hold off
% sgtitle('Total bact and yeast');
% 
% figure(2)
% subplot(4,4,Idx)
% hold on
% plot(t_tot,T_meas,'g.'); % Plot yeast measurements
% plot(t_tot,T_mod,'g'); % Plot yeast model
% ylim([0 80]);
% text = string(Sample) + ' ' + num2str(Perc_LW) + '% LW, Gam=' +...
%     param_bounds_lower(14) +...
%     ' GamL=' + param_bounds_lower(22);
% title(text,'Interpreter','None');
% hold off
% sgtitle('Total TKC, using GamL (clumped B & Y');
% 
figure(3)
subplot(4,4,Idx)
hold on
plot(t_tot,L_meas,'k.'); % Plot bacteria measurements
plot(t_tot,L_mod,'k'); % Plot bacteria model
ylim([0 800])
xlabel('t (hr)');
ylabel('Number of clumps');
text = string(Sample) + ' ' + num2str(Perc_LW) + '% LW {\it R_c}=' + param_bounds_lower(21);
% title(text,'Interpreter','None');
title(text);
hold off
% sgtitle('Number of clumps');

figure(4)
subplot(4,4,Idx)
hold on
plot(t_tot,Lb_meas,'color',[1 .3 0],'marker','.','linestyle','None'); % Plot bacteria measurements
plot(t_tot,Ly_meas,'color',[0 .5 .8],'marker','.','linestyle','None'); % Plot yeast measurements
plot(t_tot,Lb_mod,'color',[1 .3 0],'linestyle','-'); % Plot bacteria measurements
plot(t_tot,Ly_mod,'color',[0 .5 .8],'linestyle','-'); % Plot yeast model
ylim([0 5000])
xlabel('t (hr)');
ylabel('Clumped cells');
text = string(Sample) + '  at ' + num2str(Perc_LW) + '% LW';
% title(text,'Interpreter','None');
title(text);
hold off
% sgtitle('Clumped B and Y');

end