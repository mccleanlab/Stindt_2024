clearvars; clc;

%% Load and format data

row = 3;

fluorescence_measured = CompetitionRS.Data{row};
t_measured = CompetitionRS.Time_Vector{row};

% Calculate initial conditions
number_initial_measurements = 3;
initial_conditions = mean(fluorescence_measured(1:number_initial_measurements,:,:),[3 1]);

% Run model w/ parameter sets from latin hypercube

% Set number guesses (simulations) for model to perform
number_guesses = 100000;

% Set param bounds
param_bounds_lower = [.4 0.2 CompetitionRS.Kb(row)-500 1200 -1 -1];
param_bounds_upper = [.9 0.5 CompetitionRS.Kb(row)+500 1500 2 2];

% Sample latin hypercube 
params_LHS = lhsu(param_bounds_lower,param_bounds_upper,number_guesses);

% % Sample latin hypercube in log space then correct (alternative) 
% % Better when sampling over many orders of magnitude (?)
% params_LHS = lhsu(log(param_bounds_lower),log(param_bounds_upper),number_guesses);
% params_LHS = exp(params_LHS);


% Run model w/ each parameter set
fluorescence_model = zeros([size(fluorescence_measured,1) size(fluorescence_measured,2) 1 number_guesses]);
parfor guess = 1:number_guesses

    params_temp = params_LHS(guess,:);
    [~,y_temp] = ode45(@(t,y) Competition_ODEs(t,y,params_temp),t_measured,initial_conditions);
    fluorescence_model(:,:,guess) = y_temp;

end

%% Manually calculate least-squares error and rank parameter sets

% Reshape measurements and model predictions for easy matrix operations
fluorescence_model_4D = repmat(fluorescence_model,[1 1 size(fluorescence_measured,3),1]);
fluorescence_measured_4D = repmat(fluorescence_measured,[1 1 1 number_guesses]);

% Calculate RSS error between measurements and model predictions
RSS_error = (fluorescence_measured_4D - fluorescence_model_4D).^2;
RSS_error = sum(RSS_error,1:3);
RSS_error = squeeze(RSS_error);

% Rank parameters and model predictions by RSS error
[RSS_error, RSS_sort_idx] = sort(RSS_error,'ascend');
params_LHS_ranked = params_LHS(RSS_sort_idx,:);
CompetitionRS.LHS_params{row} = params_LHS_ranked;
fluorescence_model_ranked = fluorescence_model(:,:,RSS_sort_idx);

% Plot measuremetns and top model prediction
number_top_fits_plot = 5;

close all; hold on
plot(t_measured,squeeze(fluorescence_measured(:,1,:)),'r.'); % Plot bacteria measurements
plot(t_measured,squeeze(fluorescence_measured(:,2,:)),'b.'); % Plot yeast measurements
plot(t_measured,squeeze(fluorescence_model_ranked(:,1,1:number_top_fits_plot)),'r'); % Plot bacteria model
plot(t_measured,squeeze(fluorescence_model_ranked(:,2,1:number_top_fits_plot)),'b'); % Plot yeast model

%% Plot distribution of top params and all params
figure('position',[100 100 1500 500]); hold on
for jj = 1:6
    subplot(2,6,jj);
    histogram(params_LHS_ranked(:,jj),'Normalization','pdf','FaceColor','k','FaceAlpha',0.25);
    xlim([param_bounds_lower(jj),param_bounds_upper(jj)]);
    
    subplot(2,6,6+jj);
    histogram(params_LHS_ranked(1:100,jj),'Normalization','pdf','FaceColor','g','FaceAlpha',0.25);
    xlim([param_bounds_lower(jj),param_bounds_upper(jj)]);

end

%% Plot distribution of top params for all WT bact datasets
figure('position',[100 100 1500 500]); hold on
n = 1000;
params = {'Rb','Ry','Kb','Ky','Cb','Cy'};

for ii = 1:3
setparams = [CompetitionRS.Rb(ii) CompetitionRS.Ry(ii) CompetitionRS.Kb(ii) CompetitionRS.Ky(ii)];
newparams = [CompetitionRS.New_Rb(ii) CompetitionRS.New_Ry(ii) CompetitionRS.New_Kb(ii) CompetitionRS.New_Ky(ii)];
param_bounds_lower = [.4 0.2 CompetitionRS.Kb(ii)-500 1200 -1 -1];
param_bounds_upper = [.9 0.5 CompetitionRS.Kb(ii)+500 1500 2 2];


    for jj = 1:6

        param = params{jj};

        subplot(3,6,(ii-1)*6+jj);
        hold on
        histogram(CompetitionRS.LHS_params{ii,1}(1:1000,jj),'Normalization','pdf','FaceColor','r','FaceAlpha',0.1);
        histogram(CompetitionRS.LHS_params{ii,1}(1:100,jj),'Normalization','pdf','FaceColor','y','FaceAlpha',0.1);
        histogram(CompetitionRS.LHS_params{ii,1}(1:10,jj),'Normalization','pdf','FaceColor','g','FaceAlpha',0.1);
        xlim([param_bounds_lower(jj),param_bounds_upper(jj)]);
        legend('Top 10%','Top 1%','Top 0.1%');

        if jj<5
%         n = xline(newparams(jj),'color',[1 .5 0],'label','Model Median');
%         n.LabelHorizontalAlignment = 'left';
% 
        xline(setparams(jj),'Color',[0 0 1],'label','Monoculture Fit');
%         elseif jj>=5
%         med = nanmedian(CompetitionRS.LHS_params{ii,1}(:,jj));
%         xline(med,'Color',[1 .5 0],'label',['Median of all = ' num2str(med)]);
        end

        tit = [param CompetitionRS.Bacteria(ii)];
        title(tit,'Interpreter','none');

    end
end


%% Consolidate shit

C_Combined_10 = vertcat(CompetitionRS.LHS_params{1,1}(1:10,5:6),CompetitionRS.LHS_params{2,1}(1:10,5:6),...
    CompetitionRS.LHS_params{3,1}(1:10,5:6));
Cb_10 = nanmean(C_Combined_10(:,1));
Cy_10 = nanmean(C_Combined_10(:,2));
CompetitionRS.Cb_10(:) = Cb_10;
CompetitionRS.Cy_10(:) = Cy_10;

C_Combined_100 = vertcat(CompetitionRS.LHS_params{1,1}(1:100,5:6),CompetitionRS.LHS_params{2,1}(1:100,5:6),...
    CompetitionRS.LHS_params{3,1}(1:100,5:6));
Cb_100 = nanmean(C_Combined_100(:,1));
Cy_100 = nanmean(C_Combined_100(:,2));
CompetitionRS.Cb_100(:) = Cb_100;
CompetitionRS.Cy_100(:) = Cy_100;

C_Combined_1000 = vertcat(CompetitionRS.LHS_params{1,1}(1:1000,5:6),CompetitionRS.LHS_params{2,1}(1:1000,5:6),...
    CompetitionRS.LHS_params{3,1}(1:1000,5:6));
Cb_1000 = nanmean(C_Combined_1000(:,1));
Cy_1000 = nanmean(C_Combined_1000(:,2));
CompetitionRS.Cb_1000(:) = Cb_1000;
CompetitionRS.Cy_1000(:) = Cy_1000;

C_Combined_100k = vertcat(CompetitionRS.LHS_params{1,1}(:,5:6),CompetitionRS.LHS_params{2,1}(:,5:6),...
    CompetitionRS.LHS_params{3,1}(:,5:6));
Cb_100k = nanmean(C_Combined_100k(:,1));
Cy_100k = nanmean(C_Combined_100k(:,2));
CompetitionRS.Cb_100k(:) = Cb_100k;
CompetitionRS.Cy_100k(:) = Cy_100k;

for i = 1:3
% CompetitionRS.C_Combined{i,1} = C_Combined;
% 
CompetitionRS.C_Combined_10{i,1} = C_Combined_10;
CompetitionRS.C_Combined_100{i,1} = C_Combined_100;
CompetitionRS.C_Combined_1000{i,1} = C_Combined_1000;
CompetitionRS.C_Combined_100k{i,1} = C_Combined_100k;

% CompetitionRS.New_Rb(i) = nanmedian(CompetitionRS.LHS_params{i,1}(1:100,1));
% CompetitionRS.New_Ry(i) = nanmedian(CompetitionRS.LHS_params{i,1}(1:100,2));
% CompetitionRS.New_Kb(i) = nanmedian(CompetitionRS.LHS_params{i,1}(1:100,3));
% CompetitionRS.New_Ky(i) = nanmedian(CompetitionRS.LHS_params{i,1}(1:100,4));

end


%% Check out combined C terms 
Newparams = {'Cb','Cy'};
Num = [10 100 1000 100000];
NumS = {'10','100','1000','100k'};

for jj = 1:2
    Newparam = Newparams{jj};

    for k = 1:3
    N = Num(k);
    Ns = NumS{k};

    Col = ['C_Combined_' Ns];
        subplot(2,3,(jj-1)*3+k);
        hold on

        histogram(CompetitionRS.LHS_params{1,1}(1:N,jj+4),'Normalization','pdf','FaceColor','r','FaceAlpha',0.1);
        histogram(CompetitionRS.LHS_params{2,1}(1:N,jj+4),'Normalization','pdf','FaceColor','y','FaceAlpha',0.1);
        histogram(CompetitionRS.LHS_params{3,1}(1:N,jj+4),'Normalization','pdf','FaceColor','g','FaceAlpha',0.1);
%         med1 = nanmedian(CompetitionRS.LHS_params{1,1}(1:N,jj+4));
%         med2 = nanmedian(CompetitionRS.LHS_params{2,1}(1:N,jj+4));
%         med3 = nanmedian(CompetitionRS.LHS_params{3,1}(1:N,jj+4));
        mean = nanmean(CompetitionRS.(Col){3,1}(:,jj));
%         xline(med1,'Color','r','label',['Median = ' num2str(med1)]);
%         xline(med2,'Color',[1 .5 0],'label',['Median = ' num2str(med2)]);
%         xline(med3,'Color','g','label',['Median = ' num2str(med3)]);
        xline(mean,'Color','b','label',['Mean of all = ' num2str(mean)]);
        legend('WT M2','WT TR2','WT TR6')

        tit = [Newparam ', Top ' Ns ' Hits for each strain'];
        title(tit,'Interpreter','none');
    end
end