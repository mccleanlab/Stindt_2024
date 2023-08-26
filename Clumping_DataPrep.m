%% Data Prep Spark: Adjust Times

Source = AllData.Plate2;
Ends = EndPoints.Plate2;
dayTimes = unique(EndPoints.Plate2.time);

Source.FitDay(Source.time<=dayTimes(1)) = 1;
Source.FitTime(Source.FitDay==1) = (Source.time(Source.FitDay==1) - min(Source.time(Source.FitDay==1)))./3600;

for i = 2:numel(dayTimes)
    
    Source.FitDay(Source.time>dayTimes(i-1) & Source.time<=dayTimes(i)) = i;
    Source.FitTime(Source.FitDay==i) = (Source.time(Source.FitDay==i) - min(Source.time(Source.FitDay==i)))./3600;
    
end

Source = movevars(Source, 'FitDay', 'Before', 'temp');
Source = movevars(Source, 'FitTime', 'Before', 'temp');
Source = removevars(Source, {'tDiff','checkTime'});

%% Data Prep Spark: TKC

wells = unique(Source.well);

for w = 1:numel(wells)
    well = wells(w);
    
    for d = 1:numel(dayTimes)
        t = dayTimes(d);
        
        Source.TKC(Source.well==well & Source.time==t) = Ends.TKC(Ends.well==well & Ends.time==t);
        
    end
end

Source.TKC(~ismember(Source.time,dayTimes)) = nan;


%% Data Prep Spark: Reshape Everything

% (Similar to the next block used for flow data)

Pairs = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'};

Percs = [15 10 5 0];
AA = categorical("LW");
PercTit = 'PercLW';

plt = 'Plate1';
Exp = categorical("230221_DynMannose_Plus");

data = [];

for p = 1:numel(Pairs)
    Pair = categorical(cellstr(Pairs{p}));
    
    for a = 1:numel(Percs)
        G = Percs(a);
        
        wells = unique(Source.well(Source.BactYeastPair==Pair & Source.(PercTit)==G));

        Idx = find(Source.BactYeastPair==Pair & Source.(PercTit)==G);

        data0 = table();
        data0.Experiment = Exp;
        data0.Bacteria = Source.Bacteria(Idx(1));
        data0.BactFcn = Source.BactFcn(Idx(1));
        data0.Yeast = Source.Yeast(Idx(1));
        data0.YeastFcn = Source.YeastFcn(Idx(1));
        data0.BactYeastPair = Pair;
        data0.Perc_AA = Source.(PercTit)(Idx(1));
        data0.Lim_AA = AA;

        for w = 1:numel(wells)
            well = wells(w);

            Widx = find(Source.well==well & Source.BactYeastPair==Pair & Source.(PercTit)==G);

            Reads_Spark(:,1,w) = Source.Cherry(Widx);
            Reads_Spark(:,2,w) = Source.Citrine(Widx);       
            Reads_Spark(:,3,w) = Source.TKC(Widx); % Could've used Yidx too, should be repeated
            data0.Reads_Spark{1} = Reads_Spark;
            
        end
        
        Times = unique(Source.time);
        start = min(Times);
        FitTime = table();

        for t = 1:numel(Times)
            time = Times(t);
            
            Tidx = find(Source.time==time);
            
            FitTime.FitTime(t) = Source.FitTime(Tidx(1));           
            FitTime.FitDay(t) = Source.FitDay(Tidx(1));
            FitTime.CumTime(t) = (Source.time(Tidx(1))-start) / 3600;
            data0.Spark_Time{1} = FitTime;
            
        end

        data = vertcat(data,data0);

        clear Reads_Spark FitTime
    end
end

    
% Conjugation.Clumping = vertcat(Conjugation.Clumping,data);

%% Do you want to make TKC data continuous?

for row = 1:16
    
    data = Conjugation.Clumping.Reads_Spark{row,1}(:,3,:);
    data = squeeze(data);
    pos = find(~isnan(data(:,1)));
    cols = size(Conjugation.Clumping.Reads_Spark{row,1},3);
    
    for c = 1:cols
        for i = 1:numel(pos)
            r = pos(i);

            if i == 1
                data(1:r,c) = linspace(0,data(r,c),r);
            else
                prev = pos(i-1);
                start = data(prev,c)/10;
                data(prev+1:r,c) = linspace(start,data(r,c),r-prev);
            end
        end
        
        Conjugation.Clumping.Reads_Spark{row,1}(:,4,c) = data(:,c);
        
    end
    
end


%% Consolidate clump data from image analyses

% This code adds columns for 1) Clumps per total yeast, 2) Clumped yeast
% per total yeast, and 3) Clumped bact per total clumps. The "per"
% (frequency) nature of each term is necessary bc the dilutions change and
% generally there's no easy way to normalize an image range to the entire
% population of the original well. Thus, these fractional terms allow
% normalizing ODE outputs to Spark data 

% Note that this depends on still having TKC vector with nans in it!

Pairs = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'};
Percs = [15 10 5 0];

for p = 1:4
    pair = Pairs{p};
    
    for l = 1:4
        perc = Percs(l);

        destRow = find(Conjugation.Clumping.BactYeastPair==pair &...
            Conjugation.Clumping.Perc_AA==perc);
        destTKC = Conjugation.Clumping.Reads_Spark{destRow,1}(:,3,:);
        destTKC = squeeze(destTKC);
        destPos = find(~isnan(destTKC(:,1)));
        
        data = zeros(max(destPos),3,3);
        
        for day = 1:6
            
            wellList = unique(dataS.Well(dataS.BactYeastPair==pair & dataS.PercLW==perc &...
                dataS.Day==day & dataS.Mannose=='N' & dataS.Dilution~=1));
            
            if numel(wellList)==2
                wellList(3) = wellList(1);
            elseif numel(wellList)==1
                wellList(2:3) = wellList(1);
            elseif numel(wellList)==0
                wellList(1:3) = categorical("1");
            end
                
            
            pos = destPos(day);

            for w = 1:numel(wellList)
                
                
                well = wellList(w);

                Idx = find(dataS.Well==well & dataS.BactYeastPair==pair & dataS.PercLW==perc &...
                dataS.Day==day & dataS.Mannose=='N' & dataS.Dilution~=1);
                
                if ~isempty(Idx)
                data(pos,1,w) = nanmean(dataS.ClumpsPerY(Idx));
                data(pos,2,w) = nanmean(dataS.BactPerClump(Idx));
                data(pos,3,w) = nanmean(dataS.YfractClumped(Idx));
                else
                data(pos,1,w) = nan;
                data(pos,2,w) = nan;
                data(pos,3,w) = nan;
                end
                
            end
            
        end
        
        for i = 1:numel(destPos)
            r = destPos(i);

            for c = 1:3
                for w = 1:numel(wellList)
                    if i == 1
                        data(1:r,c,w) = linspace(0,data(r,c,w),r);
                    else
                        prev = destPos(i-1);
                        start = data(prev,c,w)/10;
                        data(prev+1:r,c,w) = linspace(start,data(r,c,w),r-prev);
                    end
                end
            end            
        end 
        
        Conjugation.Clumping.Reads_Spark{destRow}(:,5:7,:) = data(:,1:3,:);
        
    end
end

%% Convert freqs of clumps, clump-bact, and clump-yeast via total yeast

CherPerCit = AllFits.Saved_Params.Fluor_Conversion.Cherry / AllFits.Saved_Params.Fluor_Conversion.Citrine;

for i = 1:16
    
    % Num_Clumps_per_Tot_Yeast * Tot Yeast = num clumps in correct units:
    Conjugation.Clumping.Reads_Spark{i}(:,8,:) = Conjugation.Clumping.Reads_Spark{i}(:,5,:) .* ...
        Conjugation.Clumping.Reads_Spark{i}(:,2,:);

    % Clump_Bact_per_Clump * Num_Clumps = clump-bact in correct units:
    % (units of cherry after converting from citrine)
    Conjugation.Clumping.Reads_Spark{i}(:,9,:) = Conjugation.Clumping.Reads_Spark{i}(:,6,:) .* ...
        Conjugation.Clumping.Reads_Spark{i}(:,8,:) .* (CherPerCit);

    % Clump_Yeast_per_Tot_Yeast * Tot_Yeast = clump-yeast in correct units:
    Conjugation.Clumping.Reads_Spark{i}(:,10,:) = Conjugation.Clumping.Reads_Spark{i}(:,7,:) .* ...
        Conjugation.Clumping.Reads_Spark{i}(:,2,:);

    Conjugation.Clumping.Reads_Spark{i}(:,5:7,:) = [];
    
end



%% Get clump carrying capacities

for i = 1:size(dataS,1)
    maxYsize(i) = max(dataS.Areas{i});
end

maxYsize = sort(maxYsize,'descend');
BperClump = sort(dataS.BactPerClump,'descend');

KLy = round(sqrt(maxYsize(2)/20));
KLb = round(BperClump(2));

AllFits.Saved_Params.K.KLb('Cells_Clump') = KLb;
AllFits.Saved_Params.K.KLy('Cells_Clump') = KLy;

AllFits.Saved_Params.K.KLb('Fluor') = KLb * AllFits.Saved_Params.Fluor_Conversion.Cherry(1);
AllFits.Saved_Params.K.KLy('Fluor') = KLy * AllFits.Saved_Params.Fluor_Conversion.Citrine(1);
       