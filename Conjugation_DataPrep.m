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


%% Do you want to make TKC data continuous?

for row = 65:80
    
    data = Conjugation.ByExperiment.Reads_Spark{row,1}(:,3,:);
    data = squeeze(data);
    pos = find(~isnan(data(:,1)));
    cols = size(Conjugation.ByExperiment.Reads_Spark{row,1},3);
    
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
        
        Conjugation.ByExperiment.Reads_Spark{row,1}(:,4,c) = data(:,c);
        
    end
end

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

    
% Conjugation.ByExperiment = vertcat(Conjugation.ByExperiment,data);




%% Data Prep: Flow

Pairs = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'};

Percs = [15 10 5 0];
AA = categorical("LW");
PercTit = 'Perc_LW';

plt = 'Plate1';
Exp = categorical("230118_Dynamics");

data = [];

for p = 1:numel(Pairs)
    Pair = categorical(cellstr(Pairs{p}));
    
    for a = 1:numel(Percs)
        G = Percs(a);
        
        wells = unique(C.(plt).well(C.(plt).BactYeastPair==Pair & C.(plt).(PercTit)==G));

        Idx = find(C.(plt).BactYeastPair==Pair & C.(plt).(PercTit)==G);

        data0 = table();

        data0.Experiment = Exp;
        data0.Bacteria = C.(plt).Bacteria(Idx(1));
        data0.BactFcn = C.(plt).BactFcn(Idx(1));
        data0.Yeast = C.(plt).Yeast(Idx(1));
        data0.YeastFcn = C.(plt).YeastFcn(Idx(1));
        data0.BactYeastPair = Pair;
        data0.Perc_AA = C.(plt).(PercTit)(Idx(1));
        data0.Lim_AA = AA;


        for w = 1:numel(wells)
            well = wells(w);

            Bidx = find(C.(plt).Voltage=='Bact' & C.(plt).well==well & C.(plt).day~=0);
            Yidx = find(C.(plt).Voltage=='Yeast' & C.(plt).well==well & C.(plt).day~=0);

            Reads_Per_100Mic(:,1,w) = C.(plt).BPer100Mic(Bidx);
            Reads_Per_100Mic(:,2,w) = C.(plt).YPer100Mic(Yidx);       
            Reads_Per_100Mic(:,3,w) = C.(plt).TKC(Bidx); % Could've used Yidx too, should be repeated
            data0.Reads_Per_100Mic{1} = Reads_Per_100Mic;

            Reads_Normalized(:,1,w) = C.(plt).Bnormalized(Bidx);
            Reads_Normalized(:,2,w) = C.(plt).Ynormalized(Yidx);       
            Reads_Normalized(:,3,w) = C.(plt).TKC(Bidx); % Could've used Yidx too, should be repeated
            data0.Reads_Normalized{1} = Reads_Normalized;

        end
        
        data = vertcat(data,data0);

    end
end
       
% Conjugation.ByExperiment = vertcat(Conjugation.ByExperiment,data);

%% Consolidate experiments per condition

data = [];

% Deal with this unique case...
for i = 1:16
    Conjugation.ByExperiment.Reads_Per_100Mic{i}(5:6,:,:) = nan;
    Conjugation.ByExperiment.Reads_Normalized{i}(5:6,:,:) = nan;
end

Pairs = unique(Conjugation.ByExperiment.BactYeastPair);
Percs = unique(Conjugation.ByExperiment.Perc_AA);

for i = 1:numel(Pairs)
    pair = Pairs(i);
    for j = 1:numel(Percs)
        data0 = table();
        
        perc = Percs(j);
        
        Idx = find(Conjugation.ByExperiment.Perc_AA==perc & Conjugation.ByExperiment.BactYeastPair==pair);
        
        data0.BactYeastPair = pair;
        data0.Perc_AA = perc;
        data0.Lim_AA = categorical("LW");
        data0.Reads_Per_100Mic{1} = cat(3,Conjugation.ByExperiment.Reads_Per_100Mic{Idx});
        data0.Reads_Normalized{1} = cat(3,Conjugation.ByExperiment.Reads_Normalized{Idx});
        
        data = vertcat(data,data0);
        
    end
end

% Conjugation.Consolidated = data;
