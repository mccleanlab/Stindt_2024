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

% TKC sweep, given a set of params
clear Gam
clear params_LHS
Gam = linspace(-6,-2,100);
Gam = 10.^(Gam);
% saved_params = Conjugation.ByExperiment.Mean_Params_wDeath{Idx};
saved_params = AllFits.Saved_Params.ConjPlotParams(z,:);

for i = 1:numel(Gam)
    params_LHS(i,:) = saved_params;
end

params_LHS(:,14) = Gam;% .* AllFits.Saved_Params.Fluor_Conversion.Gamma(1);

params_LHS(:,7) = AllFits.Saved_Params.G.L('Molar') * (Perc_LW/100);
params_LHS(:,8) = AllFits.Saved_Params.G.W('Molar') * (Perc_LW/100);


% Set init conditions of B, Y, and T [B0 Y0 T0] based on experimental setups
Binit = nanmean(Conjugation.ByExperiment.Reads_Spark{Idx}(2:4,1,:),'all');
Yinit = nanmean(Conjugation.ByExperiment.Reads_Spark{Idx}(2:4,2,:),'all');
Tinit = 0;
initial_conditions = [Binit Yinit Tinit]; 

% Set up hourly time vector over 6 days
t_spark = Conjugation.ByExperiment.Spark_Time{Idx}; % For batch solver, only give hours between dilutions
dilution_days = unique(t_spark.FitDay)'; % For batch solver, enter todal days diluted
num_replicates = size(Conjugation.ByExperiment.Reads_Spark{Idx},3);

%Need to txform data to put reads only on 24-hour marks
fluorescence_measured = (Conjugation.ByExperiment.Reads_Spark{Idx}(:,[1 2 4],:));

% Run model w/ each parameter set
fluorescence_model = zeros([size(fluorescence_measured,1) size(fluorescence_measured,2) numel(Gam)]);


% Batch Solver for SPARK DATA
for guess = 1:numel(Gam)
    params_temp = params_LHS(guess,:);
    y_batch = [];
    
    for d = dilution_days
        t_measured = t_spark.FitTime(t_spark.FitDay==d,1);
        if d==1
            New_conditions = initial_conditions;
        end
        
        [~,y_add] = ode45(@(t,y) Conjugation_ODEs_wDeath(t,y,params_temp),t_measured,New_conditions);
        y_batch = vertcat(y_batch,y_add);
        
        New_conditions = y_batch(end,:,:) ./ 10;
    end
    fluorescence_model(:,:,guess) = y_batch;
end


fluorescence_model(:,3,:) = fluorescence_model(:,3,:) ./ (2*AllFits.Saved_Params.Fluor_Conversion.Citrine(1));


% Get day-end TKCs

for d = dilution_days
    row = find(t_spark.FitDay==d,1,'last');
%     EndPoints(d,1,:) = d;
    EndPoints(d,:,:) = fluorescence_model(row,:,:);
end

for g = 1:numel(Gam)
    TKCs(g,:) = EndPoints(:,3,g);
end

TKCs = flipud(TKCs);

GamLab = round(flip(Gam),2,"significant");
GamLab(11) = 0.004;
GamLab(60) = 4E-5;
GamLab(79) = 7E-6;
GamLab([2:10 12:19 21:29 31:39 41:49 51:59 61:69 71:78 80:89 91:99]) = nan;

text = string(Sample) + ' at ' + num2str(Perc_LW) + '% LW';
load('MyColormap.mat')

maxTKCs = zeros(1,16);
maxTKCs([1 2 9 10]) = [6 4 1 1];
maxColor = 20;

HighlightRow = round((size(mymap,1) * maxTKCs(z)/maxColor));
if HighlightRow > 0
    mymap(HighlightRow,:) = [0 1 1];
end



subplot(4,4,z)
H = heatmap(TKCs,'ColorMethod','None','Colormap',mymap);
H.YDisplayLabels = GamLab;
H.YLabel = 'IDC rate ({\it ?})';
H.XLabel = 'Day';
H.Title = (text);
% H.NodeChildren(3).Title.Interpreter = 'latex';
caxis([0 20]);

orig_state = warning('off','all'); %Temporarily turning off warnings bc the tiff import generates one for every file...it's a lot 

warning('off','MATLAB:structOnObject')
S = struct(H); % Undocumented
ax = S.Axes;    % Undocumented
H.GridVisible = 'off';
yline(ax,11,'color',[0.5 0.50 0.5],'LineWidth',2); % see footnotes [1,2]
yline(ax,60,'color',[1 0.40 0.0980],'LineWidth',2); % see footnotes [1,2]
yline(ax,79,'color',[1 0.40 0.0980],'LineWidth',2); % see footnotes [1,2]

% sgtitle('Predicted total transconjugants per 100uL for range of gamma values')
warning(orig_state)

end

