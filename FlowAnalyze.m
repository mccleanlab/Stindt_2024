clearvars; clc; close all
%% Draw and save gate (comment out if loading previously saved gate)
channels_to_gate = {'FSC-A','BL2-A','Include'}; % Specify pairs of channels for gating ;'FSC-A','VL1-A','Include'
channels_to_scale = {'linear','log'}; % Specify scale for each pair of channels ;'linear','log'
Yeast_gate_Lightning = draw_gate(channels_to_gate, channels_to_scale); % Draw and save gate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The above code will allow you to draw a gate on an .fcs file selected by
% UI prompt, first on a plot of FSCA-A vs SCC-A (with linear scaling) then
% on a plot of FSC-A vs FSC-H (with log scaling). It will automatically
% save the gate as a .mat file to the current folder. I recommend that you
% run this block of code once to draw gates on whatever .fcs file you
% choose and after that just comment out this section and load the saved
% gate, as shown below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test gates on several samples

gate_to_test = Yeast_gate_Lightning; % Enter which gate you wanna check out
show = test_gate(gate_to_test); 


%% Load previously saved gate (comment out if running for first time)
% BactGate = load('BactGate.mat');
% BactGate = BactGate.gate_out;
% YeastGate = load('YeastGate.mat');
% YeastGate = YeastGate.gate_out;

%% Add to previously saved gates (comment out unless doing fancy stuff)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code allows you to (1) add more gates to already saved
% gates (eg if you realized you need another level of selection), or (2)
% apply gates from one file to another file. In the latter case, be careful
% that experimental conditions (eg voltages) between files match
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_gate = Yeast_gate_Lightning; % Enter which pre-loaded gate you wanna use

dlgTitle    = 'How to Apply Gate';
dlgQuestion = 'Would you like to SUBTRACT events falling within this gate, or INCLUDE them?';
choice = questdlg(dlgQuestion,dlgTitle,'Subtract','Include','Include');
row = size(load_gate,1);
load_gate{row,5} = choice;

channels_to_gate = {'FSC-A', 'VL1-A','Include'}; % Specify pairs of channels for gating
channels_to_scale = {'linear','log'}; % Specify scale for each pair of channels
Cut_gate_Lightning = update_gate(channels_to_gate, channels_to_scale, load_gate); % Test desired 

%% Load measurments and labels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here each folder contains all .fcs files for a given plate (keep
% other .fcs files, eg, those used for calibration, in another folder).
% You can easily add folder paths as new lines or comment out folders as
% shown above. I apply labels to the measurements as they are imported by
% placing a .xlsx plate map in each of the above folders and using the
% function load_fcs('map','plate') below. Each plate map contains labels
% for each well of a 96 well plate(though you can leave empty wells blank).
% You can add an arbitray number of labels to the plate map so long as they
% include the map_* label flag and follow the 96 well format shown in the
% example plate map.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify folders from which to load .fcs files
% clearvars -except Y1 BactGate BactDblsGate;

for s = 1
    plt = ['Plate' num2str(s)];
%     bact_folder_list.Plate1 = {
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day1_Bact';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day2_Bact';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day3_Bact';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day4_Bact';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day5_Bact';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day6_Bact'};%...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day8M_Both';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day10M_Bact';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day12M_Bact'};

%     bact_folder_list.Plate2 = {
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230118 Opto2Dyn2/Flow Data/Day1_Dyn2_Bact';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230118 Opto2Dyn2/Flow Data/Day2_Dyn2_Bact';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230118 Opto2Dyn2/Flow Data/Day3_Dyn2_Bact';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230118 Opto2Dyn2/Flow Data/Day4_Dyn2_Bact'};

    yeast_folder_list.Plate1 = {
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day1_Yeast';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day2_Yeast';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day3_Yeast';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day4_Yeast';...
        '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day5_Yeast';...
        '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day6_Yeast'};...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day8M_Both';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day10M_Yeast';...
%         '/Users/kevinlauterjung/Documents/McClean/Time Courses/230321 CRISPR/Flow/Day12M_Yeast'};

% Loop through folders and load labelled measurements from .fcs files
%     for b = 1:size(bact_folder_list.(plt),1)
%         b
%         dataB_temp = load_fcs('map','plate','folder',bact_folder_list.(plt){b}); % Load and label measurements    
%         dataB_temp = add_gate(dataB_temp,Bact_gate_Lightning); % Apply gate to measurements
%         dataB_temp = format_fcsdat(dataB_temp); % Convert measurements to table
% %         dataB_temp((dataB_temp.Bact_gate==0),:) = [];
%         dataB.(plt){b} = dataB_temp; % Collect table
%     end
       
    
    for y = 1:size(yeast_folder_list.(plt),1)
        y
        dataY_temp = load_fcs('map','plate','folder',yeast_folder_list.(plt){y}); % Load and label measurements
        dataY_temp = add_gate(dataY_temp,Yeast_gate_Lightning); % Apply gate to measurements
        dataY_temp = add_gate(dataY_temp,Cut_gate_Lightning); % Apply gate to measurements
        dataY_temp = format_fcsdat(dataY_temp); % Convert measurements to table
        dataY_temp((dataY_temp.Yeast_gate_Lightning_net==0),:) = [];
        dataY.(plt){y} = dataY_temp; % Collect table
    end

    % Convert collected tables into single big table
    
%     dataB.(plt) = vertcat(dataB.(plt){:});
    dataY.(plt) = vertcat(dataY.(plt){:});

end
%% Convert relevant categories to numbers, input any extra info

catsToNum = {'day','time','TKC','DrawVol','Perc_W','Perc_U'}; %Enter any categories that should be numerical (aka doubles...default is categorical)

for s = 1
    plt = ['Plate' num2str(s)];
    
    for i = 1:numel(catsToNum)
        cat = catsToNum{i};
        
%         dataB.(plt).(cat) = double(string(dataB.(plt).(cat))); %These ones might apply to many experiments
        dataY.(plt).(cat) = double(string(dataY.(plt).(cat)));

    end
    
end

%% Create Summary Table(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the data as you like but here's a useful trick. Use grpstats() to
% do a calculation on some subset of your data. You can merge the resulting
% table back in to the primary data table so long as your keywords match
% up. Here I calculate basal expression for each strain (and replicate)
% then merge that number back into the main table and use it to calculate
% fold change. This avoids lots of complicated looping.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get rid of nonsense measurements
% data(data.BL2A<=0,:) = [];
% data(data.YL2A<=0,:) = [];

for s = 1
    plt = ['Plate' num2str(s)];
    
%     dataB.(plt).Bact_gate_net = double(dataB.(plt).Bact_gate_net);
%     dataB.(plt).Bact_gate_net(dataB.(plt).Bact_gate_net == 0) = nan;
    
%     dataB.(plt).Bact_clumps_net = double(dataB.(plt).Bact_clumps_net);
%     dataB.(plt).Bact_clumps_net(dataB.(plt).Bact_clumps_net == 0) = nan;
% 
%     dataB.(plt).TKC(isnan(dataB.(plt).TKC)) = -1; %Needed to include wells w/o TKC values. Section below fixes after the fact
%     dataB.(plt).Bubbles(isundefined(dataB.(plt).Bubbles)) = categorical("N");

        stats = {'day','well','Bacteria','BactFcn','Yeast','YeastFcn','BactYeastPair',...
        'Cytometer','SC_M9','DrawVol','TKC','Perc_W','Perc_U','Bubbles'}; %cycCat,'Contaminated',

%     dataBsum_fc.(plt) = grpstats(dataB.(plt),(stats),...
%         {'numel','nanmedian','nanmean'},'DataVars',{'Bact_gate_net','BL1A','BL2A','YL2A','VL1A'});%'VL1A',

    dataY.(plt).Yeast_gate_net = double(dataY.(plt).Yeast_gate_Lightning_net);
    dataY.(plt).Yeast_gate_net(dataY.(plt).Yeast_gate_net == 0) = nan;
    
    dataY.(plt).Cut_gate_net = double(dataY.(plt).Cut_gate_Lightning_net);
    dataY.(plt).Cut_gate_net(dataY.(plt).Cut_gate_net == 0) = nan;
% 
    dataY.(plt).TKC(isnan(dataY.(plt).TKC)) = -1;
    dataY.(plt).Bubbles(isundefined(dataY.(plt).Bubbles)) = categorical("N");

    dataYsum_fc.(plt) = grpstats(dataY.(plt),(stats),...
        {'numel','nanmedian'},'DataVars',{'Yeast_gate_net','Cut_gate_net','BL1A','BL2A','YL2A','VL1A'}); %'TKC_gate_net'

    % I rename these "B" and "Y" for bact and yeast, respectively
%     B.(plt) = dataBsum_fc.(plt);
    Y.(plt) = dataYsum_fc.(plt);

end


%% Modify and calculate

for s = 1
    plt = ['Plate' num2str(s)];
%     cycCat = cycledCats{s};

    %Fix TKC values that are wrong
%     dataB.(plt).TKC(dataB.(plt).TKC == -1) = nan;
    dataY.(plt).TKC(dataY.(plt).TKC == -1) = nan;
%     B.(plt).TKC(B.(plt).TKC == -1) = nan;
    Y.(plt).TKC(Y.(plt).TKC == -1) = nan;
    
%     dataB.(plt).(cycCat)(dataB.(plt).(cycCat) == -1) = nan;
%     dataY.(plt).(cycCat)(dataY.(plt).(cycCat) == -1) = nan;
%     B.(plt).(cycCat)(B.(plt).(cycCat) == -1) = nan;
%     Y.(plt).(cycCat)(Y.(plt).(cycCat) == -1) = nan;
%     B.(plt)(B.(plt).replicate=="0",:) = [];
%     Y.(plt)(Y.(plt).replicate=="0",:) = [];
%     B.(plt).Properties.RowNames = {};
    Y.(plt).Properties.RowNames = {};

%     B.(plt).DrawVol(:) = 20;
    Y.(plt).DrawVol(:) = 20;
    
end


%% Calculate fold change DAPI
%     Y.Plate1.Properties.VariableNames(end-1:end) = {'VL1A_fc'};
%     Y.Plate1 = join(dataY.Plate1,Y.Plate1);
%     Y.Plate1.VL1A_fc = dataY.Plate1.VL1A./dataY.Plate1.VL1A_fc;
%     data.Plate1.VL1A_fc = data.Plate1.VL1A./data.Plate1.VL1A_fc;


%% Combination station

for s = 1
    plt = ['Plate' num2str(s)];
    
    B.(plt).Voltage(:) = categorical("Bact");

    Y.(plt).Voltage(:) = categorical("Yeast");

    B.(plt).TKCVol(:) = 100;
    Y.(plt).TKCVol(:) = 100;

    Bconds.(plt) = B.(plt).Properties.VariableNames;
    Yconds.(plt) = Y.(plt).Properties.VariableNames;

    for i = 1:numel(Bconds.(plt))
        Bcond = Bconds.(plt){i};     
        if ~any(ismember(Yconds.(plt),Bcond))
            Y.(plt).(Bcond)(:) = nan;
        end        
    end
    
    for i = 1:numel(Yconds.(plt))
        Ycond = Yconds.(plt){i};        
        if ~any(ismember(Bconds.(plt),Ycond))
            B.(plt).(Ycond)(:) = nan;
        end        
    end

    
    C.(plt) = vertcat(Y.(plt),B.(plt));

end
%%
C.Plate1.Type(C.Plate1.Bacteria=='AR M2' | C.Plate1.Bacteria=='WT M2' | C.Plate1.Yeast=='yMM1585') = categorical("Cis");
C.Plate1.Type(C.Plate1.Bacteria=='AR TR2' | C.Plate1.Bacteria=='WT TR2' | C.Plate1.Yeast=='yMM1720') = categorical("Trans");

%% Normalization pt. I

% Determine baseline noise, based on "blanks": check before normalizing counts
cytometers = unique(C.Plate1.Cytometer);

for i = 1:numel(cytometers)
    
    cyt = cytometers(i);

    % Get baseline values for blank correction
    YBlanks = find(C.Plate1.Yeast=="None" & C.Plate1.Voltage=="Yeast" & C.Plate1.day~=0 & C.Plate1.Cytometer==cyt); 
    YBlankRds{i} = C.Plate1.numel_Yeast_gate_net(YBlanks);
    Ybaseline(i) = nanmean(YBlankRds{i})
        
    BBlanks = find(C.Plate1.Bacteria=="None" & C.Plate1.Voltage=="Bact" & C.Plate1.day~=0 & C.Plate1.Cytometer==cyt);
    BBlankRds{i} = C.Plate1.numel_Bact_gate_net(BBlanks);
    Bbaseline(i) = nanmean(BBlankRds{i})
   
%     ClumpBlankRds = C.Plate1.numel_Bact_clumps_net(BBlanks);
%     Clumpbaseline = nanmean(ClumpBlankRds)

end

% Ybaseline(3) = max(YBlankRds{3});
% Bbaseline(3) = min(BBlankRds{3});
%% Normalization pt. II

% Apply baselines to normalize counts, then calculate cells per 100uL

cytometers = unique(C.Plate1.Cytometer);

for i = 1:numel(cytometers)
    s = 1;
    plt = ['Plate' num2str(s)];
    
    cyt = cytometers(i);

    C.(plt).TKC(C.(plt).TKC==0) = 0.1; %This helps set the noise floor for TKC values
    
    idx = find(C.Plate1.Cytometer==cyt);

    C.(plt).BPer100Mic(idx) = 1000.*(C.(plt).numel_Bact_gate_net(idx) - Bbaseline(i))./C.(plt).DrawVol(idx); %Total bact
%     C.(plt).BPer100Mic(C.(plt).day==0) = (50*100*1.5).*(C.(plt).numel_Bact_gate_net(C.(plt).day==0) -...
%         Bbaseline)./C.(plt).DrawVol(C.(plt).day==0); %Total bact
    C.(plt).BCorrection(idx) = Bbaseline(i);
    
    C.(plt).YPer100Mic(idx) = 1000.*(C.(plt).numel_Yeast_gate_net(idx) - Ybaseline(i))./C.(plt).DrawVol(idx);
%     C.(plt).YPer100Mic(C.(plt).day==0) = (50*100).*(C.(plt).numel_Yeast_gate_net(C.(plt).day==0) -...
%         Ybaseline)./C.(plt).DrawVol(C.(plt).day==0);
    C.(plt).YCorrection(idx) = Ybaseline(i);
%     C.(plt).TKCCorrection(:) = TKCbaseline;
        
%     C.(plt).CutPer100Mic(idx) = 1000.*(C.(plt).numel_Cut_gate_net(idx) - Ybaseline(i))./C.(plt).DrawVol(idx);

    C.(plt).BPer100Mic(C.(plt).BPer100Mic<0) = 0; % If normalization resulted in neg values, set to zero
    C.(plt).Bmax_cyt(idx) = max(C.(plt).BPer100Mic(idx));
    C.(plt).Bnormalized_cyt(idx) = C.(plt).BPer100Mic(idx) ./ C.(plt).Bmax_cyt(idx);
    
    C.(plt).YPer100Mic(C.(plt).YPer100Mic<0) = 0; % If normalization resulted in neg values, set to zero
    C.(plt).CutPer100Mic(C.(plt).YPer100Mic<0) = 0; % If normalization resulted in neg values, set to zero
    C.(plt).Ymax_cyt(idx) = max(C.(plt).YPer100Mic(idx));
    C.(plt).Ynormalized_cyt(idx) = C.(plt).YPer100Mic(idx) ./ C.(plt).Ymax_cyt(idx);

end
    
    C.(plt).Bmax(:) = max(C.(plt).BPer100Mic(C.(plt).day~=0));
    C.(plt).Bnormalized = C.(plt).BPer100Mic ./ C.(plt).Bmax;
    
    C.(plt).Ymax(:) = max(C.(plt).YPer100Mic(C.(plt).day~=0));
    C.(plt).Ynormalized = C.(plt).YPer100Mic ./ C.(plt).Ymax;

%% Cell Ratios

% Calculate cell ratios by first identifying both instances of each sample
% (bact and yeast voltages)

for s = 1
    plt = ['Plate' num2str(s)];


% Start by scanning all wells, days 
    wells = unique(C.(plt).well);
    days = unique(C.(plt).day);

    for i = 1:numel(wells)
        well = wells(i);

        for j = 1:numel(days)
            day = days(j);

            % Find both instances of that sample
            idx = find(C.(plt).day==day & C.(plt).well==well);

            % We can only do this for cocultures, so need to ignore monos
            if numel(idx)==2

                % Get the bact and yeast voltage locations            
                Bidx = find(C.(plt).well==well & C.(plt).day==day & C.(plt).Voltage=="Bact");
                Yidx = find(C.(plt).well==well & C.(plt).day==day & C.(plt).Voltage=="Yeast");

                % Compute
%                 C.(plt).BYnormRatio(idx) = ...
%                     C.(plt).Bnormalized(Bidx) / C.(plt).Ynormalized(Yidx);
%                 C.(plt).BYRatio(idx) = ...
%                     C.(plt).BPer100Mic(Bidx) / C.(plt).YPer100Mic(Yidx);
                C.(plt).BYnormRatio_cyt(idx) = ...
                    C.(plt).Bnormalized_cyt(Bidx) / C.(plt).Ynormalized_cyt(Yidx);

            end
        end
    end
clear wells days
end

% C.(plt).BYnormRatio(C.(plt).BYnormRatio==Inf | C.(plt).BYnormRatio==0) = nan;
% C.(plt).BYRatio(C.(plt).BYRatio==Inf | C.(plt).BYRatio==0) = nan;
C.(plt).BYnormRatio_cyt(C.(plt).BYnormRatio_cyt==Inf | C.(plt).BYnormRatio_cyt==0) = nan;

%% TKC rates

for s = 1
    plt = ['Plate' num2str(s)];

    % Compute TKC rates per each cell type AND the cell ratios
    C.(plt).TKCPerY = (C.(plt).TKC)./(C.(plt).YPer100Mic);
    C.(plt).TKCPerB = (C.(plt).TKC)./(C.(plt).BPer100Mic);
    
%     C.(plt).TKCPerRatioNorm = (C.(plt).TKC)./(C.(plt).BYnormRatio);
%     C.(plt).TKCPerRatioNorm(C.(plt).TKCPerRatioNorm==Inf) = nan;

    C.(plt).TKCPerYnorm_cyt = (C.(plt).TKC)./(C.(plt).Ynormalized_cyt);
    C.(plt).TKCPerBnorm_cyt = (C.(plt).TKC)./(C.(plt).Bnormalized_cyt);
%     
    C.(plt).TKCPerRatioNorm_cyt = (C.(plt).TKC)./(C.(plt).BYnormRatio_cyt);
    C.(plt).TKCPerRatioNorm_cyt(C.(plt).TKCPerRatioNorm_cyt==Inf) = nan;

%     C.(plt).TKCPerRatio = (C.(plt).TKC)./(C.(plt).BYRatio);
%     C.(plt).TKCPerRatio(C.(plt).TKCPerRatio==Inf) = nan;
%     
    C.(plt).CutFreq = C.(plt).numel_Cut_gate_net ./ C.(plt).numel_Yeast_gate_net;
    C.(plt).CutFreq(C.(plt).YPer100Mic==0) = 0;
%     C.(plt).CutFreq = C.(plt).CutPer100Mic ./ C.(plt).YPer100Mic;
    
%     C.(plt).TKCPerClumpRatio = (C.(plt).TKC)./(C.(plt).ClumpsPerBact);
%     C.(plt).TKCPerClumpRatio(C.(plt).TKCPerClumpRatio==Inf) = nan;

% %     C.(plt).TKCPerY = 100.*(C.(plt).TKC./C.(plt).YPer100Mic)./C.(plt).TKCVol;
% %     C.(plt).TKCPerB = 100.*(C.(plt).TKC./C.(plt).BPer100Mic)./C.(plt).TKCVol;
% %     C.(plt).TKCPerRatio = 100.*(C.(plt).TKC./C.(plt).BYnormRatio)./C.(plt).TKCVol;
%  
% 
%     C.(plt).TKCPerCells(idx) = (100.*C.(plt).TKC(idx)./C.(plt).TKCVol(idx))./...
%         (C.(plt).BPer100Mic(Bidx) + C.(plt).YPer100Mic(Yidx));
%     C.(plt).TKCPerCellsNorm(idx) = (100.*C.(plt).TKC(idx)./C.(plt).TKCVol(idx))./...
%         (C.(plt).Bnormalized(Bidx) + C.(plt).Ynormalized(Yidx));
%     
%     C.(plt).TKCPerCells(idx) = C.(plt).TKC(idx)./(C.(plt).BPer100Mic(Bidx) + C.(plt).YPer100Mic(Yidx));
%     C.(plt).TKCPerCellsNorm(idx) = C.(plt).TKC(idx)./(C.(plt).Bnormalized(Bidx) + C.(plt).Ynormalized(Yidx));
% 
%     C.(plt).TKCPerCells(idx) = C.(plt).TKC(idx)./(C.(plt).BPerCol(Bidx) + C.(plt).YPerCol(Yidx));
%     C.(plt).TKCPerCellsNorm(idx) = C.(plt).TKC(idx)./(C.(plt).Bnormalized(Bidx) + C.(plt).Ynormalized(Yidx));

%     C.(plt).TKCPerY(C.(plt).TKCPerY==Inf) = nan;
%     C.(plt).TKCPerB(C.(plt).TKCPerB==Inf) = nan;
%     C.(plt).TKCPerB(C.(plt).TKCPerRatio==Inf) = nan;

%     C.(plt).TKCPerCells(C.(plt).TKCPerCells==Inf) = nan;
%     C.(plt).TKCPerCellsNorm(C.(plt).TKCPerCellsNorm==Inf) = nan;

% C.Plate1.CeruleanRatio = C.Plate1.numel_Cerulean_gate_net ./ C.Plate1.numel_Yeast_gate_net;
    
%     C.Plate1.time = double(string(C.Plate1.time));
%     C.Plate2.time = double(string(C.Plate2.time));
end

% C.Plate1.TKCPerBFP = (C.Plate1.TKC)./(C.Plate1.numel_TKC_gate_net);
% C.Plate1.BFPPerYeast = (C.Plate1.numel_TKC_gate_net)./(C.Plate1.numel_Yeast_gate_net);



