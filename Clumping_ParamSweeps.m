%% Fitting Proximity Term

clear EndPoints B

Idx = 8

number_guesses = 1;
CherPerCit = AllFits.Saved_Params.Fluor_Conversion.Cherry(1) / AllFits.Saved_Params.Fluor_Conversion.Citrine(1);


% For "fixed" parameters, set sampling range for solver (multipliers)
Perc_LW = Conjugation.Clumping.Perc_AA(Idx);
% Sample = Conjugation.Clumping.BactYeastPair(Idx);
Sample = Conjugation.Clumping.AltNom(Idx);

% Parameters time bounds: Rb, Ry, Kb, Ky, Cb, Cy, Gb (global AA needed by Bact),
% Gy, Ab, Ay, kb, ky, D, Gam. Middle row is actually fixed

% TKC sweep, given a set of params
clear P
clear params_LHS
P = linspace(0,3,100);
PL = 10.^(P);
% saved_params = Conjugation.Clumping.Mean_Params{Idx};
saved_params = AllFits.Saved_Params.ClumpFits_Revised(Idx,:);
for i = 1:numel(PL)
    params_LHS(i,:) = saved_params;
end
params_LHS(:,17) = PL;% .* AllFits.Saved_Params.Fluor_Conversion.Gamma(1);


% Set init conditions of B, Y, and T [B0 Y0 T0] based on experimental setups
Binit = nanmean(Conjugation.Clumping.Reads_Spark{Idx}(2:4,1,:),'all');
Yinit = nanmean(Conjugation.Clumping.Reads_Spark{Idx}(2:4,2,:),'all');
Tinit = 0;
Linit = 0;
Lbinit = 1E-20;
Lyinit = 1E-20;
initial_conditions = [Binit Yinit Tinit Linit Lbinit Lyinit]; 

% Set up hourly time vector over 6 days
t_spark = Conjugation.Clumping.Spark_Time{Idx}; % For batch solver, only give hours between dilutions
dilution_days = unique(t_spark.FitDay)'; % For batch solver, enter todal days diluted
num_replicates = size(Conjugation.Clumping.Reads_Spark{Idx},3);

%Need to txform data to put reads only on 24-hour marks
fluorescence_measured = (Conjugation.Clumping.Reads_Spark{Idx}(:,[1:2 4:7],:));

% Run model w/ each parameter set
fluorescence_model = zeros([size(fluorescence_measured,1) size(fluorescence_measured,2) numel(PL)]);


% Batch Solver for SPARK DATA
for guess = 1:numel(PL)
    params_temp = params_LHS(guess,:);
    y_batch = [];
    
    for d = dilution_days
        t_measured = t_spark.FitTime(t_spark.FitDay==d,1);
        if d==1
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
    row = find(t_spark.FitDay==d,1,'last');
    EndPoints(d,:,:) = fluorescence_model(row,:,:);
end

for g = 1:numel(PL)
    B(g,:) = EndPoints(:,1,g);
end

B = flipud(B);


PLab = round(flip(PL),2,"significant");

PLab([2:10 11:20 22:33 35:43 45:53 55:66 68:76 78:86 88:99]) = nan;
% PLab(25) = 1000;

text = string(Sample) + ' at ' + num2str(Perc_LW) + '% LW';
load('MyColormap.mat')

subplot(1,1,1)
H = heatmap(B,'ColorMethod','None','Colormap',mymap);
H.YDisplayLabels = PLab;
H.YLabel = 'Proximity Multiplier ({\it P})';
H.XLabel = 'Day';
H.Title = (text);
% H.NodeChildren(3).Title.Interpreter = 'latex';
caxis([1000 7000]);

orig_state = warning('off','all'); %Temporarily turning off warnings bc the tiff import generates one for every file...it's a lot 

warning('off','MATLAB:structOnObject')
S = struct(H); % Undocumented
ax = S.Axes;    % Undocumented
H.GridVisible = 'off';
yline(ax,44,'color',[1 0.40 0.0980],'LineWidth',2); % see footnotes [1,2]
% yline(ax,25,'color',[0.4660 1 0.1880],'LineWidth',2); % see footnotes [1,2]

% sgtitle('Predicted bacterial growth for crossB_wtY given a range of proximity factors','Interpreter','None')
warning(orig_state)



%% Gam on Gam action

clear EndPoints TKCs

for Idx = 1:16

number_guesses = 1;
CherPerCit = AllFits.Saved_Params.Fluor_Conversion.Cherry(1) / AllFits.Saved_Params.Fluor_Conversion.Citrine(1);

% For "fixed" parameters, set sampling range for solver (multipliers)
Perc_LW = Conjugation.Clumping.Perc_AA(Idx);
% Sample = Conjugation.Clumping.BactYeastPair(Idx);
Sample = Conjugation.Clumping.AltNom(Idx);

% Parameters time bounds: Rb, Ry, Kb, Ky, Cb, Cy, Gb (global AA needed by Bact),
% Gy, Ab, Ay, kb, ky, D, Gam. Middle row is actually fixed

% TKC sweep, given a set of params
clear Gam GamL
clear params_LHS
Gam = linspace(-6,-3,15);
Gam = 10.^(Gam);
GamL = Gam;

saved_params = AllFits.Saved_Params.ClumpFits_Revised(Idx,:);
for i = 1:numel(Gam)
    for j = 1:numel(GamL)
        row = numel(GamL)*(i-1) + j;
        params_LHS(row,:) = saved_params;
        params_LHS(row,14) = Gam(i);
        params_LHS(row,22) = GamL(j);
    end
end

t_nixOne = find(Conjugation.Clumping.Spark_Time{Idx}.FitDay~=1);
rows_nixed = find(Conjugation.Clumping.Spark_Time{Idx}.FitDay==1);
num_nixed = numel(rows_nixed);

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


%Need to txform data to put reads only on 24-hour marks
fluorescence_measured = (Conjugation.Clumping.Reads_Spark{Idx,1}(t_nixOne,[1:2 4:7],:));

% Run model w/ each parameter set
fluorescence_model = zeros([size(fluorescence_measured,1) size(fluorescence_measured,2) numel(GamL)*numel(Gam)]);

% Batch Solver for SPARK DATA
for guess = 1:size(params_LHS,1)
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
%         noise_min = 0.950;
%         noise_max = 1.05;
%         noise = noise_min + (noise_max - noise_min).*rand(size(y_batch(end,:,:)));
%         New_conditions = noise .* y_batch(end,:,:) ./ 10;
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

for i = 1:numel(Gam)
    for j = 1:numel(GamL)
        iter = numel(GamL)*(i-1) + j;
        TKCs(i,j) = EndPoints(6,3,iter);
    end
end

TKCs = flipud(TKCs);

GamLabY = round(flip(GamL),2,"significant");
GamLabX = round(GamL,2,"significant");

AvgTKCs = [54.3 5.7 0 0 31 18.7 18.3 14 11 .7 0 0 31.3 12.3 12 10.3];
maxC = 100;

text = string(Sample) + ' at ' + num2str(Perc_LW) + '% LW';
load('MyColormap.mat')

HighlightRows = round((size(mymap,1) * AvgTKCs(Idx)/maxC)) + [-1 0 1];
if sum(HighlightRows) < 3
    HighlightRows(HighlightRows<1) = 1;
end
mymap(HighlightRows,:) = [0 1 1] .* [1; 1; 1];

subplot(4,4,Idx)
H = heatmap(TKCs,'ColorMethod','None','Colormap',mymap);
H.YDisplayLabels = GamLabY;
H.XDisplayLabels = GamLabX;
H.YLabel = 'Free IDC Rate ({\it ?})';
H.XLabel = 'Clumped IDC Rate ({\it ?_c})';
H.Title = (text);
% H.NodeChildren(3).Title.Interpreter = 'latex';
caxis([0 maxC]);

orig_state = warning('off','all'); %Temporarily turning off warnings bc the tiff import generates one for every file...it's a lot 

warning('off','MATLAB:structOnObject')
S = struct(H); % Undocumented
ax = S.Axes;    % Undocumented
set(gca,"FontSize",8);


H.GridVisible = 'off';
% yline(ax,47,'color',[0.8500 0.3250 0.0980],'LineWidth',2); % see footnotes [1,2]
% yline(ax,rowLine,'color',[0.4660 1 0.1880],'LineWidth',2); % see footnotes [1,2]

% sgtitle('Predicted total transconjugants per 100uL for range of gamma values at Day 6')
warning(orig_state)
end
