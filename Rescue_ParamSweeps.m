%% Sweeping [AA]

Pairs = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wty'};
altPairs = {'E_{cross} S_{cross}','E_{cross} S','E S_{cross}','E S'};

clear EndPoints TKCs
CherMin = AllFits.Saved_Params.Fluor_Conversion.Cherry;
CitMin = AllFits.Saved_Params.Fluor_Conversion.Citrine;

for Idx = 1:4

number_guesses = 1;
CherPerCit = AllFits.Saved_Params.Fluor_Conversion.Cherry / AllFits.Saved_Params.Fluor_Conversion.Citrine;

% For "fixed" parameters, set sampling range for solver (multipliers)
% Sample = Pairs{Idx};
Sample = categorical(cellstr(altPairs{Idx}));


% TKC sweep, given a set of params
clear UH L
clear params_LHS

expUH = 1E-4 .* [.15 .1 .05 0];
UH = linspace(0,max(expUH),20);

expL = AllFits.Saved_Params.G.L('Molar') .* [.15 .1 .05 0];
L = linspace(0,max(expL),20);

saved_params = AllFits.Saved_Params.RescueParams(Idx,:);
for i = 1:numel(UH)
    for j = 1:numel(L)
        row = numel(UH)*(i-1) + j;
        params_LHS(row,:) = saved_params;
        params_LHS(row,8) = UH(i);
        params_LHS(row,7) = L(j);
    end
end

t_nixOne = find(Conjugation.Clumping.Spark_Time{1}.FitDay~=1);
rows_nixed = find(Conjugation.Clumping.Spark_Time{1}.FitDay==1);
num_nixed = numel(rows_nixed);

% Set init conditions of B, Y, and T [B0 Y0 T0] based on experimental setups
Binit = nanmean(Conjugation.Clumping.Reads_Spark{1}(t_nixOne(2:4),1,:),'all');
Yinit = nanmean(Conjugation.Clumping.Reads_Spark{1}(t_nixOne(2:4),2,:),'all');
Tinit = 0;
Linit = 0;
Lbinit = 1E-20;
Lyinit = 1E-20;
initial_conditions = [Binit Yinit Tinit Linit Lbinit Lyinit]; 

% Set up hourly time vector over 6 days
t_spark = Conjugation.Clumping.Spark_Time{1}; % For batch solver, only give hours between dilutions
dilution_days = unique(t_spark.FitDay)'; % For batch solver, enter todal days diluted
num_replicates = size(Conjugation.Clumping.Reads_Spark{1},3);

dilution_days(1) = [];


%Need to txform data to put reads only on 24-hour marks
fluorescence_measured = (Conjugation.Clumping.Reads_Spark{1,1}(t_nixOne,[1:2 4:7],:));

% Run model w/ each parameter set
fluorescence_model = zeros([size(fluorescence_measured,1) size(fluorescence_measured,2) numel(UH)*numel(L)]);

% Batch Solver for SPARK DATA
for guess = 1:size(params_LHS,1)
    params_temp = params_LHS(guess,:);
    y_batch = [];
    
    for d = dilution_days
        t_measured = t_spark.FitTime(t_spark.FitDay==d,1);
        if d==2
            New_conditions = initial_conditions;
        end
        
        [~,y_add] = ode45(@(t,y) Rescue_ODEs(t,y,params_temp),t_measured,New_conditions);

                
        % This block accounts for "breakaway" bact from clumps, assuming a
        % max bact/clump of 10. Any clumped bact over this value get
        % redefined as free bact 
        Lb_max = y_add;
        Lb_max(:,7) = y_add(:,5) - (10 * y_add(:,4) * CherPerCit);
        Lb_over = find(Lb_max(:,7) > 0);
        y_add(Lb_over,1) = y_add(Lb_over,1) + (Lb_max(Lb_over,7));
        y_add(Lb_over,5) = y_add(Lb_over,5) - (Lb_max(Lb_over,7));
        
        % Get rid of nonsense TKC values
        nonsense = find((y_add(:,1)<CherMin | y_add(:,2)<CitMin) |...
            (y_add(:,4)<CherMin | y_add(:,5)<CitMin));
        y_add(nonsense,3) = 0;
        
        y_batch = vertcat(y_batch,y_add);
        
        New_conditions = y_batch(end,:,:) ./ 10;
    end
    fluorescence_model(:,:,guess) = y_batch;
end
fluorescence_model(:,3,:) = fluorescence_model(:,3,:) ./ (2*AllFits.Saved_Params.Fluor_Conversion.Citrine);    
    
% Columns 1 & 2 will be FREE bact and yeast, col.s 5 & 6 CLUMPED B & Y, but
% measurements for B & Y are for TOTAL counts, so need to combine free and
% clumped versions into one vector each.
fluorescence_model(:,1,:,:) = fluorescence_model(:,1,:,:) + fluorescence_model(:,5,:,:);
fluorescence_model(:,2,:,:) = fluorescence_model(:,2,:,:) + fluorescence_model(:,6,:,:);


% Get day-end TKCs

for d = dilution_days
    Drow = find(t_spark.FitDay==d,1,'last');
    Drow = Drow - num_nixed;
    EndPoints(d,:,:) = fluorescence_model(Drow,:,:);
end

for i = 1:numel(UH)
    for j = 1:numel(L)
        iter = numel(UH)*(i-1) + j;
        TKCs(j,i) = EndPoints(6,3,iter);
    end
end

TKCs = flipud(TKCs);

LLabY = round(flip(L),2,"significant");
UHLabX = round(UH,2,"significant");

maxC = 5000;

text = string(Sample);
load('MyColormap.mat')

    
subplot(2,2,Idx)
H = heatmap(TKCs,'ColorMethod','None','Colormap',mymap);
H.YDisplayLabels = LLabY;
H.XDisplayLabels = UHLabX;
H.YLabel = '[L] (molar)';
H.XLabel = '[UH] (molar)';
H.Title = (text);
% H.NodeChildren(3).Title.Interpreter = 'latex';
caxis([0 maxC]);

orig_state = warning('off','all'); %Temporarily turning off warnings bc the tiff import generates one for every file...it's a lot 

warning('off','MATLAB:structOnObject')
S = struct(H); % Undocumented
ax = S.Axes;    % Undocumented
H.GridVisible = 'off';

for r = 1:4
    MUH = expUH(r);
    ML = expL(r);
    
    diffUH = abs(MUH - UHLabX);
    diffL = abs(ML - LLabY);
    
    [~,colLine] = min(diffUH);
    [~,rowLine] = min(diffL);
    UHLabX(colLine) = MUH;
    LLabY(rowLine) = ML;
    
    xline(ax,colLine,'color',[0.2 1 0.8],'LineWidth',2);
    yline(ax,rowLine,'color',[0.2 1 0.8],'LineWidth',2);
end

% sgtitle('Predicted transconjugants per 100uL for range of limited amino acids at Day 6')
warning(orig_state)
end
