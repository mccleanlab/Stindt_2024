%% Select all files to load
[file, folder] =  uigetfile('.tif','Select files','MultiSelect','on');

%% Convert file info into struct (fancy table)
data = [];
data0 = table();
for i = 1:length(file)
    C = strsplit(file{i},{'_','.'});
    data0.folder = cellstr(folder);
    data0.filename = file(i);
    data0.Day = categorical(C(1));
    data0.Dilution = categorical(C(2));
    data0.Well = string(C(3));
    data0.Img = double(string(C(5)));
    data = vertcat(data,data0);
end

%% Convert Stuff (could do this in Finder too...)


data.Day(data.Day=='D1') = '1';
data.Day(data.Day=='D2') = '2';
data.Day(data.Day=='D3') = '3';
data.Day(data.Day=='D4') = '4';
data.Day(data.Day=='D5') = '5';
data.Day(data.Day=='D6') = '6';
data.Day = double(string(data.Day));

data.Dilution(data.Dilution=='9Dil') = '9';
data.Dilution(data.Dilution=='10Dil') = '10';
data.Dilution(data.Dilution=='50Dil') = '50';
data.Dilution(data.Dilution=='500Dil') = '500';
data.Dilution = double(string(data.Dilution));

data.Well = erase(data.Well,"Well");
data.Well = categorical(data.Well);

%% Repeat for undiluted samples
[file, folder] =  uigetfile('.tif','Select files','MultiSelect','on');

%% Repeat for undiluted samples
datanew = [];
data0 = table();
for i = 1:length(file)
    C = strsplit(file{i},{'_','.'});
    data0.folder = cellstr(folder);
    data0.filename = file(i);
    data0.Day = categorical(C(1));
    data0.Dilution = 1;
    data0.Mannose = categorical(C(3));
    data0.Well = string(C(4));
    data0.Img = double(string(C(6)));
    datanew = vertcat(datanew,data0);
end


%% Convert Stuff (could do this in Finder too...)


datanew.Day(datanew.Day=='D1') = '1';
datanew.Day(datanew.Day=='D2') = '2';
datanew.Day(datanew.Day=='D3') = '3';
datanew.Day(datanew.Day=='D4') = '4';
datanew.Day(datanew.Day=='D5') = '5';
datanew.Day(datanew.Day=='D6') = '6';
datanew.Day = double(string(datanew.Day));

datanew.Well = erase(datanew.Well,"Well");
datanew.Well = categorical(datanew.Well);

data = vertcat(data,datanew);

%% Add info

% Load raw data

[file, folder] =  uigetfile('.xlsx','Select data file(s)','MultiSelect','on');
opts = detectImportOptions([folder file]);
opts = setvartype(opts,'char');
opts.DataRange = 'A1';

% Create well name combos
R = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
C = string(1:12);
[c, r] = ndgrid(1:numel(C),1:numel(R));
wellMat = [R(r(:)).' C(c(:)).'];
wellList = join(wellMat);
wellList = strrep(wellList,' ','');

sheets = {'ManMin','ManPlus','Dilution'};

for h = 1:numel(sheets)
    s = sheets{h};
    rawinfo.(s) = readtable([folder file], opts,'Sheet',s);
    
    InfoIdx = find(~cellfun(@isempty,rawinfo.(s){1,:})); %Get non-empty column numbers from infosheet 
    InfoCols = rawinfo.(s){1,InfoIdx}; %Get titles for those columns
    testInfo = str2double(rawinfo.(s){3,InfoIdx}); %Determine which data are numeric
    TextIdx = find(isnan(testInfo) & ~strcmp(InfoCols,'Well'));
    NumIdx = find(~isnan(testInfo));

    for i = 1:length(wellList)
        well = string(wellList{i});
        idxInfo = find(strcmp(well,rawinfo.(s){:,1}));

%         for j = TextIdx
%             Col = InfoCols{j};
%             
%             if h==1
%                 data.(Col)(data.Well == well & data.Mannose=='Min') = (rawinfo.(s){idxInfo,j});
%                 data.Plate(data.Well == well & data.Mannose=='Min') = categorical(cellstr(s));
%             elseif h==2
%                 data.(Col)(data.Well == well & data.Mannose=='Plus') = (rawinfo.(s){idxInfo,j});
%                 data.Plate(data.Well == well & data.Mannose=='Plus') = categorical(cellstr(s));
%             elseif h==3
%                 data.(Col)(data.Well == well & data.Dilution~=1) = (rawinfo.(s){idxInfo,j});
%                 data.Plate(data.Well == well & data.Dilution~=1) = categorical(cellstr(s));
%             end
%         end
        
        for k = NumIdx
            Col = InfoCols{k};
            if h==1
                data.(Col)(data.Well == well & data.Mannose=='Min') = str2double(rawinfo.(s){idxInfo,k});
            elseif h==2
                data.(Col)(data.Well == well & data.Mannose=='Plus') = str2double(rawinfo.(s){idxInfo,k});
            elseif h==3
                data.(Col)(data.Well == well & data.Dilution~=1) = str2double(rawinfo.(s){idxInfo,k});
            end
        end
    end
end

%     data.BactYeastPair = categorical(data.BactYeastPair);
%     data.BactFcn = categorical(data.BactFcn);
%     data.YeastFcn = categorical(data.YeastFcn);
%     data.Bacteria = categorical(data.Bacteria);
%     data.Yeast = categorical(data.Yeast);

%% Import and analyze
orig_state = warning('off','all'); %Temporarily turning off warnings bc the tiff import generates one for every file...it's a lot 

for j = 1:size(data,1)
%     tic
    time = j/(size(data,1));
    waitbar(time)
        
        clearvars -except data file folder orig_state j

        
    %Import each Tiff file, binarize, fill, for both green then red channels
    path = [data.folder{j} data.filename{j}];

    t = Tiff(path);
    
    rawR = read(t);  
    factR = imbinarize(rawR,'adaptive','ForegroundPolarity','bright','Sensitivity',0.4);
    fillR = imfill(factR,'holes');

    PropsR = regionprops('table',fillY,'Area');
    
    setDirectory(t,2);
    rawY = read(t);  
    factY = imbinarize(rawY,'adaptive','ForegroundPolarity','bright','Sensitivity',0.4);
    fillY = imfill(factY,'holes');
    
    PropsY = regionprops('table',fillY,'Area');
    Areas = PropsY.Area(PropsY.Area>10);
    
    if numel(Areas) > 0 && numel(Areas) < 1000
        data.Areas{j} = Areas(:);
        data.NumEvents(j) = categorical("InRange");
    elseif numel(Areas) < 3 || numel(PropsR) < 3
        data.Areas{j} = [];
        data.NumEvents(j) = categorical("FewEvents");
    elseif numel(Areas) > 1000
        data.Areas{j} = [];
        data.NumEvents(j) = categorical("Oversaturated");
    end

    
%     toc
end
warning(orig_state)

% Make a 'select' (or 'small') version of your table, with only the entries
% that had data w/in range
dataS = data(data.NumEvents=='InRange',:);

%% Import and analyze
orig_state = warning('off','all'); %Temporarily turning off warnings bc the tiff import generates one for every file...it's a lot 

for j = 242:size(dataS,1)
    tic
    j
    time = j/(size(dataS,1));
    waitbar(time)
        
        clearvars -except dataS file folder orig_state j

        
    %Import each Tiff file, binarize, fill, for both green then red channels
    path = [dataS.folder{j} dataS.filename{j}];

    t = Tiff(path);
    setDirectory(t,2);
    rawY = read(t);  
    factY = imbinarize(rawY,'adaptive','ForegroundPolarity','bright','Sensitivity',0.4);
    fillY = imfill(factY,'holes');
    
    PropsY = regionprops('table',fillY,'Area','BoundingBox');
    PropsY = PropsY(PropsY.Area>10,:);
    
    setDirectory(t,1);
    rawR = read(t);
    factR = imbinarize(rawR,'adaptive','ForegroundPolarity','bright','Sensitivity',0.5);
    
    PropsR = regionprops('table',factR,'Centroid');
    
    for p = 1:size(PropsY,1)
        B = PropsY.BoundingBox(p,:);
        PropsY.Robj(p) = 0;
        
        for k = 1:size(PropsR,1)
            C = PropsR.Centroid(k,:);
            
            test = (C(1)>B(1) && C(1)<B(1)+B(3) && C(2)>B(2) && C(2)<B(2)+B(4));
            
            if test==1
                PropsY.Robj(p,k) = 1;
%             else
%                 PropsY.Robj(p,k) = 0;
            end
                        
        end
        
        PropsY.Robjs(p) = sum(PropsY.Robj(p,:));
        
    end
        
   dataS.Robjs{j} = PropsY.Robjs(:);
   
   toc
end
warning(orig_state)


%% Filter out any that only have a couple of events in either channel

%This is important to get rid of optical defects (which screw with
%thresholding, making it seem there are no cells) and overly-dense images,
%which also show up as no cells (thresholding finds all intensity values
%the same throughout)

for j = 1:size(dataS,1);

    time = j/(size(dataS,1));
    waitbar(time)

    path = [dataS.folder{j} dataS.filename{j}];

    t = Tiff(path);
%     setDirectory(t,2);
    rawR = read(t);  
    factR = imbinarize(rawR,'adaptive','ForegroundPolarity','bright','Sensitivity',0.5);
    fillR = imfill(factR,'holes');

    PropsR = regionprops('table',fillR,'Area');
    if numel(PropsR) < 3
        dataS.Kill(j) = categorical("Y");
    end

    setDirectory(t,2);
    rawY = read(t);  
    factY = imbinarize(rawY,'adaptive','ForegroundPolarity','bright','Sensitivity',0.4);
    fillY = imfill(factY,'holes');

    PropsY = regionprops('table',fillY,'Area');
    if numel(PropsY) < 3
        dataS.Kill(j) = categorical("Y");
    end
    
end
    
%% Check any out

for j = [1158 2188 2444]

    path = [dataS.folder{j} dataS.filename{j}];

    t = Tiff(path);
%     setDirectory(t,2);
    rawR = read(t);  
    factR = imbinarize(rawR,'adaptive','ForegroundPolarity','bright','Sensitivity',0.5);
    fillR = imfill(factR,'holes');

    setDirectory(t,2);
    rawY = read(t);  
    factY = imbinarize(rawY,'adaptive','ForegroundPolarity','bright','Sensitivity',0.4);
    fillY = imfill(factY,'holes');

    
subplot(2,2,1);
imshow(factY);
title('Binarized Y');

subplot(2,2,2);
imshow(fillY);
title('Filled Y');

subplot(2,2,3);
imshow(factR);
title('Binarized R');

subplot(2,2,4);
imshow(fillR);
title('Filled R');

sgtitle(['Row ' num2str(j) ', ' dataS.filename{j}],'Interpreter','None')

pause

end

%% Counts

for j = 1:size(dataS,1)
    dataS.Ycount(j) = numel(dataS.Areas{j});
    
    dataS.YeastWBact(j) = numel(find(dataS.Robjs{j}>0));
    dataS.TotBactOnYeast(j) = sum(dataS.Robjs{j});
   
    dataS.AvgYeastArea(j) = mean(dataS.Areas{j});
    dataS.MedianYeastArea(j) = median(dataS.Areas{j});
end

dataS.AvgBactPerYeast = dataS.TotBactOnYeast ./ dataS.YeastWBact;
dataS.FreqYeastWBact = dataS.YeastWBact ./ dataS.Ycount;
%% Normalize sizes to Yevent counts

for i = 1:size(dataS,1)
    clear a
    s = numel(dataS.Areas{i});
    val = dataS.Day(i);
    a(1:s) = val;
    dataS.DayA{i} = a;
end

%% Bring in Tecan data

% This assumes Tecan data saved as table "EndPoints", from
% TecAnalyzeSparkMulti_KL. Modify the first variables conds based on
% which data from the Spark table that you want paired with image data

% conds = ({'TKC','OD600','Cherry','Citrine','BFP'});
conds = Cols;

for j = 1:size(dataS,1)
    
    clear row
%     b = dataS.Bacteria(j);
%     y = dataS.Yeast(j);
    d = dataS.Day(j);
%     a = dataS.PercLW(j);
    w = dataS.WellSource(j);
    m = dataS.Mannose(j);
    
%     row = find(EndPoints.Plate23.Bacteria==b & EndPoints.Plate23.Yeast==y & EndPoints.Plate23.Mannose==m &...
%         EndPoints.Plate23.PercLW==a & EndPoints.Plate23.day==d & EndPoints.Plate23.well==w);
    row = find(EndPoints.Plate23.Mannose==m & EndPoints.Plate23.day==d & EndPoints.Plate23.well==w);
    
    if numel(row)==1
        for i = 1:numel(conds)
            cond = conds{i};
            dataS.(cond)(j) = EndPoints.Plate23.(cond)(row);
        end
    else
        dataS.Check(j) = 1;
    end
    
end

%% Set everything else to nans
% CombinedIdx = vertcat(HighColIdx, LoColIdx);
OthersIdx = 1:size(data,1);
OthersIdx(CombinedIdx) = [];

for k = 1:numel(OthersIdx)
    j = OthersIdx(k);
        for i = 1:numel(Bconds)         
            cond = Bconds{i};
            data.(cond)(j) = nan;
        end
        
        for i = 1:numel(Yconds)
            Ycond = Yconds{i};
            data.(Ycond)(j) = nan;
        end
        
end

%% Get counts of clumped yeast

% First decide how to define single yeast vs. clumped. Smallest yeast are
% around 20 px^2. If scaling this up (2 similarly sized yeast together),
% need to square theh increase, too, e.g. 2 yeast is 80 px^2, 3 is 180 px^2

base = 20;
dataS.YbaselineArea(:) = base;

for i = 1:size(dataS,1)
    Areas = dataS.Areas{i};
    biggest = max(Areas);
    maxcount = round(sqrt(biggest/base));
    
    ClumpCount(Areas > 0 & Areas < base*4) = 1;
    
    if maxcount>1
        for c = 2:maxcount
            ClumpCount(Areas > base*(c-1)^2 & Areas <= base*c^2) = c-1;
        end

        ClumpCount(Areas > base*c^2) = c;
    end
    
    dataS.YclumpCount{i} = ClumpCount';
    clear ClumpCount
end

%% Modify your clumped bact counts to be reflective of new yeast-size definition

for i = 1:size(dataS,1)
    yeasts = dataS.YclumpCount{i};
    bact = dataS.Robjs{i};
    bact(yeasts>=2 & bact==0) = 1; %By definition, everyclump has 1 bact
    dataS.Bclumped(i) = nansum(bact(yeasts>=2));
end


%% Add up clumpy yeast per image, get % clumped

% for i = 1:size(dataS,1)
%     ClumpCount = dataS.YclumpCount{i};
%     
%     dataS.Yclumped(i) = sum(ClumpCount(ClumpCount>=2));
%     dataS.Ysinglets(i) = sum(ClumpCount(ClumpCount<2));
%     dataS.Ytot(i) = dataS.Yclumped(i) + dataS.Ysinglets(i);
%     dataS.YfractClumped(i) = dataS.Yclumped(i)/dataS.Ytot(i);
%     dataS.NumClumps(i) = numel(ClumpCount(ClumpCount>=2));
%     dataS.ClumpsPerY(i) = dataS.NumClumps(i) / dataS.Ytot(i);
%     clear ClumpCount
% 
% end

% dataS.AvgYperClump = dataS.Yclumped ./ dataS.NumClumps;
dataS.BactPerClump = dataS.Bclumped ./ dataS.NumClumps;
dataS.BactPerClump(isinf(dataS.BactPerClump)) = 0;
% dataS.BactPerClump(isnan(dataS.BactPerClump)) = 0;