% clearvars -except data rawinfo
% close all
% clc

% Before running, decide whether you're going to import sample info. If so,
% make a new sheet in your Excel file (should be called 'Sheet1' by default) 
% and enter sample info there (see example for formatting). Can also enter
% experiment info here, if desired:

Experiment = {'230321'}; % Enter some kinda title
numPlates = 2; % How many plates are included in the raw data?...
    %Should be based on actual Spark reads (so even if a plate gets read twice, count as 2, for now)
nP = numPlates;

% Importing sample info?

 dlgTitle    = 'Import Option';
 dlgQuestion = ({'Are you importing sample info?','(Label info for each plate as "Sheet#")'});
 choice = questdlg(dlgQuestion,dlgTitle,'Yes','Nah','Old', 'Yes');

% Load raw data

[file, folder] =  uigetfile('.xlsx','Select data file(s)','MultiSelect','on');
opts = detectImportOptions([folder file]);
opts = setvartype(opts,'char');
opts.DataRange = 'A1';
sheetData = 'SparkControl magellan Sheet 1';
%
rawdata = readtable([folder file], opts,'Sheet',sheetData);
if choice == 'Yes'
    for i = nP
%         if i == 3 %Only for double Spark reads/plate, comment out if loop for almost all cases!
%             i=2;
%         end
        sheet = ['Sheet' num2str(i)];
        sheetInfo{i} = sheet;
        plt = ['Plate' num2str(i)];
        AllInfo.(plt) = readtable([folder file], opts,'Sheet',sheetInfo{i});
    end
end


% Create well name combos
R = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
C = string(1:12);
[c, r] = ndgrid(1:numel(C),1:numel(R));
wellMat = [R(r(:)).' C(c(:)).'];
wellList = join(wellMat);
wellList = strrep(wellList,' ','');

% Organize rawdata into table 'data'
clearvars data well idxList data0
orig_state = warning('off','all'); %Temporarily turn off warnings

% The following finds all measurements for Spark data that's exported as
% "amend to existing file" in Excel. It does NOT apply to single Spark
% "experiments" (continuous readings) exported once as a spreadsheet
IdxReads = find(strcmp('SPARK',rawdata{:,1})); 
totReads = size(IdxReads,1)/nP;
warning('on')
if rem(size(IdxReads,1),nP)~=0
    warning('Your total read count is not a multiple of your number of plates. This could mean your raw data requires editing, and that your plates will be incorrectly assigned to raw reads!')
end

% Define time 0
ext = '(?<hour>\d+):(?<min>\d+):(?<sec>\d+)';
% time = regexp(rawdata{IdxReads(1)-4,1},ext,'names');
% t0 = str2double(time{1,1}.hour)*3600 + str2double(time{1,1}.min)*60 + str2double(time{1,1}.sec);

% Use the first date to figure out all the day #s (doesn't account for
% changing months)
exd = '(?<year>\d+)-(?<month>\d+)-(?<day>\d+)';
% day = regexp(rawdata{IdxReads(1)-4,1},exd,'names');
% dayMod = str2double(day{1,1}.day) - 1; 

% Enter the reading types exactly as they're listed in the Spark data
channelList.Plate1 = {'OD600','Cherry','Citrine','BFP'}; 
chans.Plate1 = size(channelList.Plate1,2);
channelList.Plate2 = {'mCherry','mVenus','mTurq','OD600','OD700'}; 
chans.Plate2 = size(channelList.Plate2,2);

%For each measurement, there are as many instances of "<>" as channels
%read. The following loop goes through all measurements and finds each
%instance of a plate read, and assigns it to a column based on channel
%(column1 = channel1, etc.)

idxDat(1,:) = (find(strcmp('<>',rawdata{1:IdxReads(1),1})))';

for j = 2:length(IdxReads)
    clear idx0
    idx0 = nan(1,5);
    idx0 = (find(strcmp('<>',rawdata{IdxReads(j-1):IdxReads(j),1}))+IdxReads(j-1)-1)';
    if numel(idx0) < size(idxDat,2)
        colB = numel(idx0)+1;
        colE = size(idxDat,2);
        idx0(1,colB:colE) = nan;
    elseif numel(idx0) > size(idxDat,2)
        row = size(idxDat,1);
        colE = numel(idx0);
        colB = size(idxDat,2)+1;
        idxDat(row,colB:colE) = nan;
    end
    idxDat = vertcat(idxDat,idx0);
end


% And then we split those idexes into different plates in a struct. Each
% struct (named eg "Plate1") has every read as a row, each channel as a
% column

for i = 1:nP
    plt = ['Plate' num2str(i)];
    IdxDat.(plt) = idxDat((1:totReads)*nP-(nP-i),:);
    
    DelCol = find(isnan(IdxDat.(plt)(1,:)));
    IdxDat.(plt)(:,DelCol) = [];
end


%% Multiple Spark reads per plate? (not true for most experiments!)

DblRd = [4]; %Ok, so which plates were read twice? Enter as vector based on actual reads (not position in Fluent)

for i = DblRd
    plt = ['Plate' num2str(i)];
    pltDbl = ['Plate' num2str(i+1)];
    
    lastCol = size(IdxDat.(plt),2);
    IdxDat.(plt)(:,lastCol+1) = IdxDat.(pltDbl)(:,1);
end

%% Bring it in! Compile data

orig_state = warning('off','all'); %Temporarily turn off warnings

Rds = 1:nP;

if exist('DblRd')==1
    Rds(DblRd+1) = [];
end

for s = Rds
    s
    
    plt = ['Plate' num2str(s)];
    IdxPlt = IdxDat.(plt);
%     chans = size(IdxPlt,2);
    chans = size(channelList.(plt),2);

    clear data 
    data = table();

    for j = 1:size(IdxPlt,1) 
        tbar = j/size(IdxPlt,1);
        waitbar(tbar)
        
        clearvars time day

        %Enter the time (same for every well)
        day = regexp(rawdata{IdxPlt(j,chans)+9,1},exd,'names');
        day = str2double(day{1,1}.day);
        
            %Accounts for changing months, eg 9/30 to 10/2
            if day < dayMod 
                dayfull = day - dayMod + lastday - 1;
            else
                dayfull = day - dayMod - 1;
                lastday = day;
            end

        time = regexp(rawdata{IdxPlt(j,chans)+9,1},ext,'names');
        timefull = (str2double(time{1,1}.hour)*3600 + str2double(time{1,1}.min)*60 + ...
            str2double(time{1,1}.sec))+(86400*(dayfull))-t0;

        for k = 1:length(wellList) % For every well...
            row = k+length(wellList)*(j-1);

            %Enter the time (same for every well)
            data.time(row) = timefull;
            %Enter the temp
            data.temp(row) = 30;
            %Enter the well name
            data.well(row) = wellList(k);
            %Enter the day
            data.day(row) = dayfull;
            
%             %Get rid of the following
%             if j == size(IdxPlt,1)
%                 cut = 0;
%             else
%                 cut = any(strcmp(rawdata{IdxPlt(j,1):IdxPlt(j+1,1),1},'CUT'));
%             end
%             data.CutSite(row) = cut;
% 
% 
            for i = 1:chans % For each channel...

                channel = channelList.(plt){i};

                %Get the location of that well for the given time point:
                idxR = find(strcmp(wellMat(k,1),rawdata{IdxPlt(j,i):IdxPlt(j,i)+9,1}));
                idxC = find(strcmp(wellMat(k,2),rawdata{IdxPlt(j,i),:}));

                %Enter the channel's reading
                data.(channel)(row) = str2double(rawdata{IdxPlt(j,i)+idxR-1,idxC});
%                 data.GFP(row) = str2double(rawdata{IdxGFP(j,1)+idxR-1,idxC});

            end
        end
    end
    data = movevars(data, 'well', 'Before', 'time');
    
    AllData.(plt) = data;
end
    
warning(orig_state)

%% Import sample info to table

for s = 2
    plt = ['Plate' num2str(s)];
    clear data rawinfo
    data = AllData.(plt);
    rawinfo = AllInfo.(plt);

    if choice ~= 'Nah'

        InfoIdx = find(~cellfun(@isempty,rawinfo{1,:})); %Get non-empty column numbers from infosheet 
        InfoCols = rawinfo{1,InfoIdx}; %Get titles for those columns
        testInfo = str2double(rawinfo{3,InfoIdx}); %Determine which data are numeric
        TextIdx = find(isnan(testInfo) & ~strcmp(InfoCols,'Well'));
        NumIdx = find(~isnan(testInfo));

        for i = 1:length(wellList)
            well = string(wellList{i});
            idxInfo = find(strcmp(well,rawinfo{:,1}));

            for j = TextIdx
                Col = InfoCols{j};
                data.(Col)(data.well == well) = (rawinfo{idxInfo,j});
            end
            for k = NumIdx
                Col = InfoCols{k};
                data.(Col)(data.well == well) = str2double(rawinfo{idxInfo,k});
                
            end
        end
    end 
    
%     data.BactYeastPair = categorical(data.BactYeastPair);
%     data.BactFcn = categorical(data.BactFcn);
%     data.YeastFcn = categorical(data.YeastFcn);
%     data.Bacteria = categorical(data.Bacteria);
%     data.Yeast = categorical(data.Yeast);

    AllData.(plt) = data;
    AllData.(plt).well = categorical(AllData.(plt).well);

end
%% This section tests time intervals between reads. Especially useful for verifying distribution over multiple plates 

for s = 1
    plt = ['Plate' num2str(s)];
    for i = 2:size(AllData.(plt),1)
        AllData.(plt).tDiff(1) = 0;
        AllData.(plt).tDiff(i) = AllData.(plt).time(i) - AllData.(plt).time(i-1);
    end
    AllData.(plt).checkTime(AllData.(plt).tDiff > 2350) = 1;
    
    %This spits out rows in your data table to check bc the time gaps are unusually large
    clear i n
    i = find(AllData.(plt).checkTime == 1);
    n = numel(i);
    RowsToCheck(1:n,s) = i; 
    
    %This gets the "end point" times for each day, assuming all time gaps
    %are day-end changes
    for j = 1:n
        DayEnds(j,s) = AllData.(plt).time(RowsToCheck(j,s)-1);
        DayEnds(j+1,s) = max(AllData.(plt).time);
    end
end
DayEndsTime = DayEnds./3600; %Converts times to hour format so you can check against raw data. CHECK FIRST BEFORE DOING ENDPOINTS!

%% Copy out end-point measurements (mostly for TKC stuff)

% First check your "DayEnds" variable to make sure every time point is the
% end of each day
% DayEnds = DayEndsTime.*3600; %Converts back in case of changes made (e.g. rows deleted that weren't actually endpoints)

for s = 1
    plt = ['Plate' num2str(s)];
    EndPoints.(plt) = table();
    Data = AllData.(plt);
    
    for i = 1:size(DayEnds,1)
        EndPoints0 = table();
        EndPoints0 = Data(Data.time == DayEnds(i,s),:);
        EndPoints.(plt) = vertcat(EndPoints.(plt),EndPoints0);
    end
    
    EndPoints.(plt).TKC(:) = nan; %Just getting the table variable going, still have to copy in values
    EndPoints.(plt) = movevars(EndPoints.(plt), 'TKC', 'Before', 'OD600');
end


%% Calculate End Point TKC values

%Do this after making EndPoints table(s) and transferring TKC raw data in

for s = 1
    plt = ['Plate' num2str(s)];
    clear maxY maxR maxB
    maxY = max(EndPoints.(plt).Citrine);
    maxR = max(EndPoints.(plt).Cherry);
    maxB = max(EndPoints.(plt).BFP);
    
    EndPoints.(plt).TKC(EndPoints.(plt).TKC==-1) = nan;
    
    EndPoints.(plt).TKCpCitrine = 2*(EndPoints.(plt).TKC ./ EndPoints.(plt).Citrine);
    EndPoints.(plt).TKCpCherry = 2*(EndPoints.(plt).TKC ./ EndPoints.(plt).Cherry);
    EndPoints.(plt).TKCpBFP = 2*(EndPoints.(plt).TKC ./ EndPoints.(plt).BFP);
    
    EndPoints.(plt).NormCher = EndPoints.(plt).Cherry ./ maxR;
    EndPoints.(plt).NormCit = EndPoints.(plt).Citrine ./ maxY;
    EndPoints.(plt).NormBFP = EndPoints.(plt).BFP ./ maxB;
    
    EndPoints.(plt).BYRatio = EndPoints.(plt).Cherry(:) ./ EndPoints.(plt).Citrine(:);
    EndPoints.(plt).NormBYRatio = EndPoints.(plt).NormCher(:) ./ EndPoints.(plt).NormCit(:);
    
    EndPoints.(plt).BFPpCit = EndPoints.(plt).BFP(:) ./ EndPoints.(plt).Citrine(:);
    EndPoints.(plt).NormBFPpCit = EndPoints.(plt).NormBFP(:) ./ EndPoints.(plt).NormCit(:);

end

%% Contamination?

%This section is only to flag samples demonstrated elsewhere eg flow data
%as contaminated

EndPoints.Plate1.Contamination(EndPoints.Plate1.well == "A9" | EndPoints.Plate1.well == "A10" | ...
    EndPoints.Plate1.well == "B10" | EndPoints.Plate1.well == "C10" | EndPoints.Plate1.well == "D9" | ...
    EndPoints.Plate1.well == "E9" | EndPoints.Plate1.well == "E12" | EndPoints.Plate1.well == "F9" | ...
    EndPoints.Plate1.well == "G9") = 1;


