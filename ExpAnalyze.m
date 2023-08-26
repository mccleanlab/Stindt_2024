%% Use this to import previously saved data, otherwise skip this section!!

[file, folder] =  uigetfile('.tif','Select files','MultiSelect','on'); %First window is to select source images, folder
% load(uigetfile); %Second window allows you to select data file

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
    data0.Cells = categorical(C(2));
    data0.Perc_LW = double(string(C(3)));
    data0.Colony = double(string(C(4)));
    data = vertcat(data,data0);
end

%% Convert Days (could do this in Finder too...)


data.Day(data.Day=='Day0') = '0';
data.Day(data.Day=='Day1') = '1';
data.Day(data.Day=='Day2') = '2';
data.Day(data.Day=='Day3') = '3';
data.Day(data.Day=='Day4') = '4';
data.Day(data.Day=='Day5') = '5';
data.Day(data.Day=='Day6') = '6';
data.Day = double(string(data.Day));

%% Expand information into full names. Hardcoding joy
data.Bacteria(data.Cells=="AM" | data.Cells=="AM15" | data.Cells=="AM16") = categorical("AR M2");
data.Bacteria(data.Cells=="WM" | data.Cells=="WM15" | data.Cells=="WM16") = categorical("WT M2");
data.Bacteria(data.Cells=="AT" | data.Cells=="AT17" | data.Cells=="AT16") = categorical("AR TR6");
data.Bacteria(data.Cells=="WT" | data.Cells=="WT17" | data.Cells=="WT16") = categorical("WT TR6");

data.BactFcn(data.Cells=="AM" | data.Cells=="AM15" | data.Cells=="AM16" |...
    data.Cells=="AT" | data.Cells=="AT17" | data.Cells=="AT16") = categorical("crossB");
data.BactFcn(data.Cells=="WM" | data.Cells=="WM15" | data.Cells=="WM16" |...
    data.Cells=="WT" | data.Cells=="WT17" | data.Cells=="WT16") = categorical("wtB");
data.BactFcn(data.Cells=="15" | data.Cells=="16" | data.Cells=="17") = categorical("None");

data.Yeast(data.Cells=="15" | data.Cells=="AM15" | data.Cells=="WM15") = categorical("yMM1585");
data.Yeast(data.Cells=="16" | data.Cells=="AM16" | data.Cells=="WM16" | data.Cells=="AT16" |...
    data.Cells=="WT16") = categorical("yMM1636");
data.Yeast(data.Cells=="17" | data.Cells=="AT17" | data.Cells=="WT17") = categorical("yMM1720");

data.YeastFcn(data.Cells=="15" | data.Cells=="AM15" | data.Cells=="WM15" |...
    data.Cells=="17" | data.Cells=="AT17" | data.Cells=="WT17") = categorical("crossY");
data.YeastFcn(data.Cells=="16" | data.Cells=="AM16" | data.Cells=="WM16" | data.Cells=="AT16" |...
    data.Cells=="WT16") = categorical("wtY");
data.Yeast(data.Cells=="AM" | data.Cells=="WM" | data.Cells=="AT" | data.Cells=="WT") = categorical("None");

data.BactYeastPair(data.Cells=="AM15" | data.Cells=="AT17") = categorical("crossB_crossY");
data.BactYeastPair(data.Cells=="WM15" | data.Cells=="WT17") = categorical("wtB_crossY");
data.BactYeastPair(data.Cells=="AM16" | data.Cells=="AT16") = categorical("crossB_wtY");
data.BactYeastPair(data.Cells=="WM16" | data.Cells=="WT16") = categorical("wtB_wtY");
data.BactYeastPair(data.Cells=="AM" | data.Cells=="AT") = categorical("crossB");
data.BactYeastPair(data.Cells=="15" | data.Cells=="17") = categorical("crossY");
data.BactYeastPair(data.Cells=="16") = categorical("wtY");
data.BactYeastPair(data.Cells=="WM" | data.Cells=="WT") = categorical("wtB");

data.Bacteria(isundefined(data.Bacteria)) = categorical("None");
data.BactFcn(isundefined(data.BactFcn)) = categorical("None");
data.Yeast(isundefined(data.Yeast)) = categorical("None");
data.YeastFcn(isundefined(data.YeastFcn)) = categorical("None");

%% Radial metrics: import, circle select, intensities, Li coloc.
orig_state = warning('off','all'); %Temporarily turning off warnings bc the tiff import generates one for every file...it's a lot 


%Do first: enter conversions between pixels and real colony length
%How to get: open image(s) in Fiji, then select Image->Properties (or hit
%Shift+Cmd+P). The "Pixel Width" field gives length/pixel. 

% data.MicronsPerPixel(data.Day == 0) = 5.1600;
data.MicronsPerPixel(data.Day ~= 0) = 8.0625;

% Radii = [];
% Radii0 = table();
% cols = {'filename','Bacteria','BactFcn','Yeast','YeastFcn','BactYeastPair','Day','Perc_LW','RadiusSelect',...
%     'Center','RadiusPixels'}; %,'Pixels2Microns','RadiusMillimeters'
% J = 2001+round(rand(1,4)*928);   



for j = 601
    j
    tic
    time = j/(size(data,1));
    waitbar(time)
        
        clearvars -except data file folder orig_state j fix

        
    %Import each Tiff file, binarize, fill, for both green then red channels
    path = [data.folder{j} data.filename{j}];

    [rawG factG fillG rawR factR fillR] = TiffImport(path);

    
    %Find features of colony: center, radius, etc. by passing "filled"
    %versions of binarized image to the ColID
    
    [Center Radius CircSelect] = ColID(fillG, fillR);
    data.Center{j} = Center;
    data.RadiusPixels(j) = Radius;
    data.CircleSelect(j) = CircSelect;

    
    % Convert radius pixels to real length units
    
    data.RadiusMillimeters(j) = data.RadiusPixels(j) * data.MicronsPerPixel(j) / 1000;

    
    
    % "BGdrop" matrixes are the pixel intensities after correcting (dropping off) background
    
    BGdropGreen = double(rawG) .* factG; 
    BGdropRed = double(rawR) .* factR;
%     Red2Green = (BGdropRed)./(BGdropGreen); %Ratio of all red pixels to green pixels
        
    % Make the zeros beyond the radius 'nan' instead
    for row = 1:size(BGdropGreen,1)
        for col = 1:size(BGdropRed,2)
            if round(sqrt((row-Center(1))^2+(col-Center(2))^2)) > Radius
                BGdropGreen(row,col) = "NaN";
                BGdropRed(row,col) = "NaN";
            end
        end
    end
    
    %Basic metrics on colony: total red/green, etc.
        data.TotGreen(j) = nansum(BGdropGreen,'all');     %Some colony-wide metrics
        data.TotRed(j) = nansum(BGdropRed,'all');
        data.AvgGreen(j) = mean(BGdropGreen,'all','omitnan');
        data.AvgRed(j) = mean(BGdropRed,'all','omitnan');
    
        
        
    %Take line scans of colony intensities, get info per radius. This
    %function takes a lot of time, so skip if you don't want these
 
    RadialMetrics = RadialScan(Center, Radius, CircSelect, BGdropGreen, BGdropRed);
    
    
    
    %Transfer those metrics into data table
    if ~isempty(RadialMetrics)
        Mets = RadialMetrics.Properties.VariableNames;
        for i = 1:numel(Mets)
            met = char(Mets(i));
            data.(met)(j) = RadialMetrics.(met)(1);
        end
    end
    
    
    
    
    % Li colocalization for entire colony
        
    if CircSelect ~= "Not enough info"

        DiffG = (BGdropGreen - data.AvgGreen(j));
        DiffR = (BGdropRed - data.AvgRed(j));
        DiffProd = DiffG .* DiffR;
        Pos = DiffProd > 0;

        data.ICQ(j) = (nnz(Pos) / numel(DiffProd)) - 0.5;
        data.LiSum(j) = nansum(DiffR, 'all') * nansum(DiffG, 'all');

    else

        data.ICQ(j) = nan;
        data.LiSum(j) = nan;

    end
    
    
    
    toc
end
warning(orig_state)

% data.CircleSelect = cellstr(data.CircleSelect);
        
        % This breaks out the intensity profiles per radius and saves them
        % to a new table. Each row of "Radii" is thus 1 pixel further from
        % the center of the circle, and each "profile" is a "shell" of
        % intensities at that radius. Takes a lot of RAM to make, so leave
        % commented out unless you really wanna go crazy
        
%         for r = 1:size(IRt,1)
%             
%             for c = 1:numel(cols)
%             col = cols{c};
%             Radii0.(col)(1:r) = data.(col)(j);
%             end
%         
%         Radii0.rPixels(r) = r;
%         Radii0.RedProfiles{r} = IRt(r,:);
%         Radii0.GreenProfiles{r} = IGt(r,:);
%         Radii0.R2GProfiles{r} = IRatt(r,:);
%         end
%         
%         Radii = vertcat(Radii,Radii0);
%     end
%     toc
%     clearvars -except Radii Radii0 folder file data cols orig_state
% end


%% Check out specific circles if wanted

%Ignore this section if you just wanna crank data

%This block uses current TIFFs (not saved to table), good for working
%through above section
subplot(1,3,1)
imshow(BGdropGreen) %Can also enter in a specific row from data table
title('Green')
viscircles(Center,Radius)

subplot(1,3,2)
imshow(BGdropRed)
title('Red')
viscircles(Center,Radius)

subplot(1,3,3)
imshow(Red2Green)
title('R:G Ratio')
viscircles(Center,Radius)

%% Scan through colonies

%This block lets you go through everything that has been saved to table. It
%displays in the title 1) the row of the image, 2) which channel was used
%to make the circle, and 3) the experimental condition

%While going through, make sure to note which circles need to be drawn
%manually, as well as which channel would be best for defining the circle

    %Can select good examples or bad ones, up to you
%     GoodIdx = find(data.CircleSelect~="Not enough info" & data.Day~=0);
%     CheckIdx = find(data.CircleSelect=="Not enough info");
    
    for i = 1:numel(fix) %Make sure to enter here which one you're using
    
%     j = GoodIdx(i);
%     j = CheckIdx(i);
    j = fix(i);
    
    path = [data.folder{j} data.filename{j}];
    [rawG factG fillG rawR factR fillR] = TiffImport(path);
    
    BGdropGreen = double(rawG) .* factG; 
    BGdropRed = double(rawR) .* factR;
    Combo = (BGdropRed)+(BGdropGreen); %Ratio of all red pixels to green pixels

    Center = data.Center{j};
    Radius = data.RadiusPixels(j);
    CircSel = data.CircleSelect(j);
    
    subplot(1,3,1);
    imshow(BGdropGreen); %Can also enter in a specific row from data table
    title('Green');
    viscircles(Center,Radius);

    subplot(1,3,2);
    imshow(BGdropRed);
    title('Red');
    viscircles(Center,Radius);

    subplot(1,3,3);
    imshow(Combo);
    title('R + G');
    viscircles(Center,Radius);
    
    row = num2str(j);
    sgtitle([row CircSel data.filename{j}],'Interpreter','None')
    
    pause
    end
    
%% Flag any that had issues

    %Find some in there that need fixing after all? Manually enter those (rows) into a
    %variable called "fix", then run the following
    
    %You can do this several ways, but I make a fresh "fix" variable for
    %each channel, then run them separately
    
    for i = 1:numel(fix)
        row = fix(i);
        if data.CircleSelect(row)=="Not enough info"
            strjoin(['Womp womp you wrote down this row wrong',num2str(row)])
        else
            data.CircleSelect(row) = categorical("Redo");
        end
    end
    
%% Fix bad circles by drawing new ones

FixIdx = find((data.CircleSelect=="Redo" | data.CircleSelect=="Check circles") & data.Day~=0);

%This will display green, red, and combined channels, and allow you to
%manually draw and modify a circle to whichever works best (probably
%combined), then saves the center and radius, also labeling the Circle
%Selection as "Manual"

%First, though, you need to choose which channel to draw a circle onto
dlgTitle    = 'Circle channel';
dlgQuestion = 'Are you drawing circles on Green, Red, or Combo channel?';
choice = questdlg(dlgQuestion,dlgTitle,'Combo','Green','Red','Combo');
choice = categorical(cellstr(choice));

for i=1:numel(fix)
    clear drawC
%     j = FixIdx(i);
%     j = CheckIdx(i);
    j = fix(i);
    
    path = [data.folder{j} data.filename{j}];
    [rawG factG fillG rawR factR fillR] = TiffImport(path);
    
    BGdropGreen = double(rawG) .* factG; 
    BGdropRed = double(rawR) .* factR;
%     Red2Green = (BGdropRed)./(BGdropGreen); %Ratio of all red pixels to green pixels

    %Takes all the '1's of green, adds to those of red
    fillRat = fillR + fillG;
    
    %Display channels
    subplot(1,3,1);
    imshow(BGdropGreen); %Can also enter in a specific row from data table
    title('Green');
    viscircles(data.Center{j},data.RadiusPixels(j));

    subplot(1,3,2);
    imshow(BGdropRed);
    title('Red');
    viscircles(data.Center{j},data.RadiusPixels(j));

    subplot(1,3,3);
    imshow(fillRat);
    title('Combined fill');
    
    row = num2str(j);
    sgtitle([row data.CircleSelect(j) data.filename{j} 'Draw and modify circle, then hit any button to proceed'],'Interpreter','None')

    %Draw and modify circle
    if choice=='Combo'
        drawC = drawcircle(subplot(1,3,3));
        Circ = categorical("ManualCombo");
        
    elseif choice == 'Green'
        drawC = drawcircle(subplot(1,3,1)); 
        Circ = categorical("ManualGreen");

    elseif choice == 'Red'
        drawC = drawcircle(subplot(1,3,2));
        Circ = categorical("ManualRed");
    end

    pause
   
    if ~isempty(drawC)
        data.Center{j} = drawC.Center;
        data.Radius(j) = drawC.Radius; 
        data.CircleSelect(j) = Circ;
    end
     
end    
    
%% Get metrics with manually-drawn circles

orig_state = warning('off','all'); %Temporarily turning off warnings bc the tiff import generates one for every file...it's a lot 

ManIdx = find(data.CircleSelect == "ManualCombo" | data.CircleSelect == "ManualGreen" | data.CircleSelect == "ManualRed");
% AllGoodIdx = find(data.CircleSelect ~= "Not enough info");

for k = 1:numel(fix)
    k
    tic
    time = k/(numel(ManIdx));
    waitbar(time)
        
%     clearvars -except data AllGoodIdx k orig_state ManIdx

%     j = ManIdx(k);
    j = fix(k);
%     j = AllGoodIdx(k);
%     j = 795;
    
    
    %Import each Tiff file, binarize, fill, for both green then red channels
    path = [data.folder{j} data.filename{j}];

    [rawG factG fillG rawR factR fillR] = TiffImport(path);

    
    %For manually drawn circles, details already saved, need to PULL them
    %from 'data', rather than save to 'data', as above
    
    Center = data.Center{j};
    Radius = data.RadiusPixels(j);
    CircSelect = data.CircleSelect(j);

    
    % Convert radius pixels to real length units
    
    data.RadiusMillimeters(j) = data.RadiusPixels(j) * data.MicronsPerPixel(j) / 1000;

    
    
    % "BGdrop" matrixes are the pixel intensities after correcting (dropping off) background
    
    BGdropGreen = double(rawG) .* factG; 
    BGdropRed = double(rawR) .* factR;
%     Red2Green = (BGdropRed)./(BGdropGreen); %Ratio of all red pixels to green pixels
        
    % Make the zeros beyond the radius 'nan' instead
    for row = 1:size(BGdropGreen,1)
        for col = 1:size(BGdropRed,2)
            if round(sqrt((row-Center(1))^2+(col-Center(2))^2)) > Radius
                BGdropGreen(row,col) = "NaN";
                BGdropRed(row,col) = "NaN";
            end
        end
    end
    
    %Basic metrics on colony: total red/green, etc.
        data.TotGreen(j) = nansum(BGdropGreen,'all');     %Some colony-wide metrics
        data.TotRed(j) = nansum(BGdropRed,'all');
        data.AvgGreen(j) = mean(BGdropGreen,'all','omitnan');
        data.AvgRed(j) = mean(BGdropRed,'all','omitnan');
    
        
        
    %Take line scans of colony intensities, get info per radius. This
    %function takes a lot of time, so skip if you don't want these
 
    RadialMetrics = RadialScan(Center, Radius, CircSelect, BGdropGreen, BGdropRed);
    
    
    
    %Transfer those metrics into data table
    if ~isempty(RadialMetrics)
        Mets = RadialMetrics.Properties.VariableNames;
        for i = 1:numel(Mets)
            met = char(Mets(i));
            data.(met)(j) = RadialMetrics.(met)(1);
        end
    end
    
    
    
    
    % Li colocalization for entire colony
        
    if CircSelect ~= "Not enough info"

        DiffG = (BGdropGreen - data.AvgGreen(j));
        DiffR = (BGdropRed - data.AvgRed(j));
        DiffProd = DiffG .* DiffR;
        Pos = DiffProd > 0;

        data.ICQ(j) = (nnz(Pos) / numel(DiffProd)) - 0.5;
        data.LiSum(j) = nansum(DiffR, 'all') * nansum(DiffG, 'all');

    else

        data.ICQ(j) = [];
        data.LiSum(j) = [];

    end
    
    
    
    toc
end
warning(orig_state)

%% Get AUC/area

data.TotGreenPerArea = data.TotGreen ./ (pi .* (data.RadiusMillimeters).^2);
data.TotRedPerArea = data.TotRed ./ (pi .* (data.RadiusMillimeters).^2);

%% Make radius variable that matches size of each sample's radial metrics

for j = 1:size(data,1)
    
    s = numel(data.RadialAvgR{j}); %Arbitrary choice of metric to base size on, but all should be the same size
    data.RadiusPixelsV{j} = (1:s);
    data.RadiusMillimetersV{j} = (1:s) .* data.MicronsPerPixel(j) / 1000;
end

%% If you have time course data that correspond to these images, port that into 'data'

% This assumes flow data saved as consolidated summary table "C", from
% nomoflowjo_KL. Modify the first 2 variables Bconds and Yconds based on
% which data from the summary table that you want paired with image data

% Note that which variable to assign to B vs Y cond depends mostly on where
% the valuable info is for that variable in your flow data
Bconds = ({'TKC','numel_Bact_gate_net','nanmean_YL2A','nanmedian_YL2A','BPer100Mic','BCorrection',...
    'Bmax','Bnormalized','BYRatio','BYnormRatio','TKCPerB','TKCPerRatio','TKCPerRatioNorm'});

Yconds = ({'numel_Yeast_gate_net','nanmean_BL1A','nanmedian_BL1A','nanmean_VL1A','nanmedian_VL1A','YPer100Mic','YCorrection',...
    'Ymax','Ynormalized','TKCPerY','BFPbaseline','FoldBFP'});

HighColIdx = find(data.Colony>3 & data.Day~=0);
LoColIdx = find(data.Colony<4 & data.Day==6);
CombinedIdx = vertcat(HighColIdx, LoColIdx);

for k = 1:numel(CombinedIdx)
%     j = HighColIdx(k);
%     j = LoColIdx(k);
    j = CombinedIdx(k);
    
    b = data.Bacteria(j);
    y = data.Yeast(j);
    c = data.Colony(j);
    d = data.Day(j);
    a = data.Perc_LW(j);
    
    rows = find(C.Plate1.Bacteria==b & C.Plate1.Yeast==y & C.Plate1.colony==c &...
        C.Plate1.Perc_LW==a);
    
    if numel(rows)==2
        Brow = find(C.Plate1.Bacteria==b & C.Plate1.Yeast==y & C.Plate1.colony==c &...
            C.Plate1.Perc_LW==a & C.Plate1.Voltage=="Bact" & C.Plate1.day==d);
        Yrow = find(C.Plate1.Bacteria==b & C.Plate1.Yeast==y & C.Plate1.colony==c &...
            C.Plate1.Perc_LW==a & C.Plate1.Voltage=="Yeast" & C.Plate1.day==d);
        
        for i = 1:numel(Bconds)
            Bcond = Bconds{i};
            data.(Bcond)(j) = C.Plate1.(Bcond)(Brow);
        end
        
        for i = 1:numel(Yconds)
            Ycond = Yconds{i};
            data.(Ycond)(j) = C.Plate1.(Ycond)(Yrow);
        end

        
    elseif C.Plate1.Voltage(rows)=="Bact"
        
        for i = 1:numel(Bconds)
            Bcond = Bconds{i};
            data.(Bcond)(j) = C.Plate1.(Bcond)(rows);
        end

    elseif C.Plate1.Voltage(rows)=="Yeast"
        
        for i = 1:numel(Yconds)
            Ycond = Yconds{i};
            data.(Ycond)(j) = C.Plate1.(Ycond)(rows);
        end
        
%     elseif isempty(rows)==1
        
    end
end

%% Set everything else to nans
% CombinedIdx = vertcat(HighColIdx, LoColIdx);
OthersIdx = 1:size(data,1);
OthersIdx(CombinedIdx) = [];

for k = 1:numel(OthersIdx)
    j = OthersIdx(k);
        for i = 1:numel(Bconds)         
            Bcond = Bconds{i};
            data.(Bcond)(j) = nan;
        end
        
        for i = 1:numel(Yconds)
            Ycond = Yconds{i};
            data.(Ycond)(j) = nan;
        end
        
end

%% Threshold tinkering

orig_state = warning('off','all'); %Temporarily turning off warnings bc the tiff import generates one for every file...it's a lot 

for j = 131
    tic

    %Import each Tiff file, binarize, fill, for both green then red channels
    path = [data.folder{j} data.filename{j}];
    
    t = Tiff(path);
    rawR = read(t);
    factR = imbinarize(rawR); %Works fine for BG correct
    factR3 = imbinarize(rawR,'adaptive','ForegroundPolarity','bright','Sensitivity',0.3);
    factR2 = imbinarize(rawR,'adaptive','ForegroundPolarity','bright','Sensitivity',0.2);
    factR2_5 = imbinarize(rawR,'adaptive','ForegroundPolarity','bright','Sensitivity',0.25);

    fillR = imfill(factR,'holes');
    
    setDirectory(t,2);
    rawG = read(t);  
    factG = imbinarize(rawG);
    factG3 = imbinarize(rawG,'adaptive','ForegroundPolarity','bright','Sensitivity',0.3);
    factG2 = imbinarize(rawG,'adaptive','ForegroundPolarity','bright','Sensitivity',0.2);
    factG2_5 = imbinarize(rawG,'adaptive','ForegroundPolarity','bright','Sensitivity',0.25);
    fillG = imfill(factG,'holes');

    subplot(2,4,1);
    imshow(factR);
    title('Otsu');

    subplot(2,4,2);
    imshow(factR2);
    title('Sens=0.2','Interpreter','None');

    subplot(2,4,3);
    imshow(factR2_5);
    title('Sens=0.25','Interpreter','None');

    subplot(2,4,4);
    imshow(factR3);
    title('Sens=0.3','Interpreter','None');

    subplot(2,4,5);
    imshow(factG);
    title('Otsu','Interpreter','None');

    subplot(2,4,6);
    imshow(factG2);
    title('Sens=0.2','Interpreter','None');

    subplot(2,4,7);
    imshow(factG2_5);
    title('Sens=0.25','Interpreter','None');

    subplot(2,4,8);
    imshow(factG3);
    title('Sens=0.3','Interpreter','None');

    sgtitle([num2str(j) ' ' data.filename{j}],'Interpreter','None')

    %For manually drawn circles, details already saved, need to PULL them
    %from 'data', rather than save to 'data', as above
    
%     Center = data.Center{j};
%     Radius = data.RadiusPixels(j);
%     CircSelect = data.CircleSelect(j);
% 
%     
%     % Convert radius pixels to real length units
%     
%     data.RadiusMillimeters(j) = data.RadiusPixels(j) * data.MicronsPerPixel(j) / 1000;
% 
%     
%     
%     % "BGdrop" matrixes are the pixel intensities after correcting (dropping off) background
%     
%     BGdropGreen = double(rawG) .* factG; 
%     BGdropRed = double(rawR) .* factR;
% %     Red2Green = (BGdropRed)./(BGdropGreen); %Ratio of all red pixels to green pixels
%         
%     % Make the zeros beyond the radius 'nan' instead
%     for row = 1:size(BGdropGreen,1)
%         for col = 1:size(BGdropRed,2)
%             if round(sqrt((row-Center(1))^2+(col-Center(2))^2)) > Radius
%                 BGdropGreen(row,col) = "NaN";
%                 BGdropRed(row,col) = "NaN";
%             end
%         end
%     end
%     
%     %Basic metrics on colony: total red/green, etc.
%         data.TotGreen(j) = nansum(BGdropGreen,'all');     %Some colony-wide metrics
%         data.TotRed(j) = nansum(BGdropRed,'all');
%         data.AvgGreen(j) = mean(BGdropGreen,'all','omitnan');
%         data.AvgRed(j) = mean(BGdropRed,'all','omitnan');
%     
%         
%         
%     %Take line scans of colony intensities, get info per radius. This
%     %function takes a lot of time, so skip if you don't want these
%  
%     RadialMetrics = RadialScan(Center, Radius, CircSelect, BGdropGreen, BGdropRed);
%     
%     
%     
%     %Transfer those metrics into data table
%     if ~isempty(RadialMetrics)
%         Mets = RadialMetrics.Properties.VariableNames;
%         for i = 1:numel(Mets)
%             met = char(Mets(i));
%             data.(met)(j) = RadialMetrics.(met)(1);
%         end
%     end
%     
%     
%     
%     
%     % Li colocalization for entire colony
%         
%     if CircSelect ~= "Not enough info"
% 
%         DiffG = (BGdropGreen - data.AvgGreen(j));
%         DiffR = (BGdropRed - data.AvgRed(j));
%         DiffProd = DiffG .* DiffR;
%         Pos = DiffProd > 0;
% 
%         data.ICQ(j) = (nnz(Pos) / numel(DiffProd)) - 0.5;
%         data.LiSum(j) = nansum(DiffR, 'all') * nansum(DiffG, 'all');
% 
%     else
% 
%         data.ICQ(j) = [];
%         data.LiSum(j) = [];
% 
%     end
%     
%     
%     
%     toc
end
warning(orig_state)

%% BFP tinkering

orig_state = warning('off','all'); %Temporarily turning off warnings bc the tiff import generates one for every file...it's a lot 

for k = 1:numel(scan)
    j = scan(k);

    %Import each Tiff file, binarize, fill, for both green then red channels
    path = [data.folder{j} data.filename{j}];
    
    t = Tiff(path);
    rawR = read(t);
    factR4 = imbinarize(rawR,'adaptive','ForegroundPolarity','bright','Sensitivity',0.4);
    fillR = imfill(factR,'holes');
    
    setDirectory(t,2);
    rawG = read(t);  
    factG = imbinarize(rawG);

    setDirectory(t,3);
    rawB = read(t); 
    factB = imbinarize(rawB);
    factB5 = imbinarize(rawB,'adaptive','ForegroundPolarity','bright','Sensitivity',0.5);
    
    rawB2R = rawB ./ rawR;
    factB2R = imbinarize(rawB2R);
    
    factB2RFlip = abs(factB2R - 1);
    
    subplot(2,3,1);
    imshow(factR4);
    title('Red');

    subplot(2,3,2);
    imshow(factG);
    title('Green','Interpreter','None');
    
    subplot(2,3,3);
    imshow(factB5);
    title('Blue');

    subplot(2,3,4);
    imshow(rawB2R,'DisplayRange',[0 10]);
    title('Raw BFP/RFP');

    subplot(2,3,5);
    imshow(factB2R);
    title('Bin. BFP/RFP');

    subplot(2,3,6);
    imshow(factB2RFlip);
    title('Flipped Bin. BFP/RFP');



    sgtitle([num2str(j) ' ' data.filename{j}],'Interpreter','None')
    
   
    pause
end

%% BFP Metrics

orig_state = warning('off','all'); %Temporarily turning off warnings bc the tiff import generates one for every file...it's a lot 

%Set your thresholding sensitivity! Runs from 0 to 1
sens = 0.5;
% data.Bsens(:) = sens;

for j = 1:size(data,1)
    j
    tic
    time = j/(size(data,1));
    waitbar(time)
        
    clearvars -except data file folder orig_state j sens

    %Import each Tiff file, binarize, fill, for both green then red channels
    path = [data.folder{j} data.filename{j}];
    
    t = Tiff(path);
%     rawR = read(t);
%     factR = imbinarize(rawR);
%     fillR = imfill(factR,'holes');
%     
%     setDirectory(t,2);
%     rawG = read(t);  
%     factG = imbinarize(rawG);

    setDirectory(t,3);
    rawB = read(t); 
    factB = imbinarize(rawB,'adaptive','ForegroundPolarity','bright','Sensitivity',sens);
    
%     B2R = rawB./rawR;
%     B2R = double(B2R);
%     
%     B2G = rawB./rawG;
%     B2G = double(B2G);

    Center = data.Center{j};
    Radius = data.RadiusPixels(j);
    CircSelect = data.CircleSelect(j);
%     
%     % "BGdrop" matrixes are the pixel intensities after correcting (dropping off) background
%     
%     BGdropGreen = double(rawG) .* factG; 
%     BGdropRed = double(rawR) .* factR;
    BGdropBlue = double(rawB) .* factB;
%   
% 
%     % Make the zeros beyond the radius 'nan' instead
%     for row = 1:size(BGdropGreen,1)
%         for col = 1:size(BGdropRed,2)
%             if round(sqrt((row-Center(1))^2+(col-Center(2))^2)) > Radius
%                 BGdropGreen(row,col) = "NaN";
%                 BGdropRed(row,col) = "NaN";
%                 BGdropBlue(row,col) = "NaN";
%                 B2R(row,col) = nan;
%                 B2G(row,col) = nan;
%             end
%         end
%     end
%     
          
%     PropsB = regionprops('table',factB,'Area','Centroid');
%     
%     %Nix events outside radius, too
%     for p = 1:size(PropsB,1)
%         Ctrd = PropsB.Centroid(p,:);
%         
%         if round(sqrt((Ctrd(1)-Center(1))^2+(Ctrd(2)-Center(2))^2)) > Radius
%             PropsB.Kill(p) = 1;
%         end
%     end
%     
%     PropsB(PropsB.Kill==1,:) = [];
%     
%     data.BAreas{j} = PropsB.Area;
%     data.BAreasFilt{j} = PropsB.Area(PropsB.Area>1 & PropsB.Area<100);
%     data.BNumEvents(j) = numel(PropsB.Area);
%     data.BNumEventsFilt(j) = numel(PropsB.Area(PropsB.Area>1 & PropsB.Area<100));
%     data.BAvgSize(j) = mean(PropsB.Area);
%     data.BAvgSizeFilt(j) = mean(PropsB.Area(PropsB.Area>1 & PropsB.Area<100));
%     data.BMedSize(j) = median(PropsB.Area);
%     data.BMedSizeFilt(j) = median(PropsB.Area(PropsB.Area>1 & PropsB.Area<100));
%     data.BCentroids{j} = PropsB.Centroid;
    
    %Basic metrics on colony: total red/green, etc.
%         data.TotBlue(j) = nansum(BGdropBlue,'all');    
%         data.AvgBlue(j) = mean(BGdropBlue,'all','omitnan');
%         data.TotB2R(j) = nansum(B2R,'all');     
%         data.AvgB2R(j) = mean(B2R,'all','omitnan');
%         data.TotB2G(j) = nansum(B2G,'all');     
%         data.AvgB2G(j) = mean(B2G,'all','omitnan');
% 
%         
%     %Take line scans of colony intensities, get info per radius. This
%     %function takes a lot of time, so skip if you don't want these
%  
    if CircSelect ~= "Not enough info"
        theta = linspace(0, 360, 4*pi*(Radius)); % More than needed to avoid gaps.
        xA = Center(1) + (Radius) * cosd(theta); %"Radius + 10" used to give a little extra buffer region
        yA = Center(2) + (Radius) * sind(theta);
        xy = round([xA', yA']);
        xy = unique(xy,'rows');
        xA = xy(:, 1);
        yA = xy(:, 2);

        
        % Take intensity profiles from center to circumference
        
        I = zeros(length(xy),round(Radius));
        for coord = 1 : length(xy)
            IB0 = improfile(BGdropBlue,...
                [Center(1) xA(coord)],[Center(2) yA(coord)])';
            IB(coord,1:length(IB0)) = IB0;
        end


        %Here's your intensity profile for each. Rows are "r" (radius),
        %columns are theta (angle). As a check, notice that the top row of
        %each (r = 0) has the same value for every column.
        
        IBt = IB'; % Blue
        
        %But the closer to the center this is, the more duplicate reads
        %there are (the center intensity gets logged again for every
        %circumference read, e.g.). The interval of "correct" reads to keep
        %is closer to C/(2*pi*r), where C is the outer circumference, whose
        %size is based on the number of xy coordinates.
        
        R = size(IBt,1);
        IBtM = IBt;
        
        IBtM(1,2:size(IBtM,2)) = nan; %First row manually, just take 1 reading

        for r = 2:(size(IBt,1)-1)
            idx = zeros(1,size(IBt,2)); %Set an empty idx at first

            if r <= R/2
                int = round(R/r); 
                idx(1:int:end) = 1; %Idx every ith to keep
                IBtM(r,~idx) = nan; %Set all else in row to nan
            elseif r > R/2
                int = round(1/(1 - (r/R))); %Get the mirror fraction around 1/2
                idx(2:int:end) = 1; %Idx every ith to keep
                IBtM(r,idx==1) = nan; %Now ONLY set ith to nan
            end
            
        end
%             
%         %Now do calculations on those profiles (mean, variance, SD) per
%         %radius
%         
% %         AvgIB = table();
% %         TotIB = table();
        for n = 1:size(IBt,1)
%             AvgIB(n) = nanmean(IBtM(n,:));
            TotIB(n) = nansum(IBtM(n,:));
        end
%         
%         data.RadialAvgB{j} = AvgIB;
        data.RadialTotB{j} = TotIB;
%         
    else 
        data.RadialAvgB{j} = [];
    end
     
    
    toc
end
warning(orig_state)


%% Combination station

for i = 1:size(newdata,1)
    file = newdata.filename(i);
    
    idx = find(strcmp(data.filename,file));
    
    if isempty(idx)
        newdata.add(i) = 1;
    end
end
