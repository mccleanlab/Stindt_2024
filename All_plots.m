% This file contains plot code at the top and data prep at the bottom

%% F1C: Growth data examples

AllData.Plate1 = Complete(1).Plate1.Tecan;

maxCh = max(AllData.Plate1.Cherry(AllData.Plate1.day>1));
maxCit = max(AllData.Plate1.Citrine(AllData.Plate1.day>1));

clear g;
figure('position',[0 0 600 300]);

g(1,1) = gramm('x',AllData.Plate1.time,'y',AllData.Plate1.Cherry./maxCh,...
    'linestyle',AllData.Plate1.YeastFcn,'size',AllData.Plate1.YeastFcn,...
    'subset',(AllData.Plate1.Bacteria=='LeuA, TrpR' & ...
    AllData.Plate1.Yeast~='yMM1636' & AllData.Plate1.LW_Percent==15));
g(1,1).stat_summary('type','std');
g(1,1).set_title('Tecan Fluorescence Data');
g(1,1).set_names('x','Day','y','Normalized FU','linestyle','Coculture Pair','row','LW%');
g(1,1).set_order_options('linestyle',{'CrossfeederY','None'},'size',{'None','CrossfeederY'});
g(1,1).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,1).set_color_options('map',[1 0 0],'legend','expand');
g(1,1).no_legend();

g(1,1).update('y',AllData.Plate1.Citrine./maxCit,'linestyle',AllData.Plate1.BactFcn,...
    'size',AllData.Plate1.BactFcn,'subset',(AllData.Plate1.Yeast=='yMM1585' & ...
    AllData.Plate1.Bacteria~='wtB' & AllData.Plate1.LW_Percent==15));
g(1,1).stat_summary('type','std');
g(1,1).set_names('x','Day','y','Normalized FU','linestyle','Coculture Pair');
g(1,1).set_order_options('linestyle',{'CrossfeederB','None'},'size',{'None','CrossfeederB'})
g(1,1).set_color_options('map',[.65 .5 0],'legend','expand');
g(1,1).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,1).no_legend();


g(1,2) = gramm('x',Complete(1).All.day,'y',Complete(1).All.Bnormalized,...
    'linestyle',Complete(1).All.YeastFcn,'size',Complete(1).All.YeastFcn,...
    'subset',(Complete(1).All.BactFcn=="crossB" & Complete(1).All.day~=0 &...
    Complete(1).All.YeastFcn~="wtY" & Complete(1).All.Assay=="Dynamics" &...
     Complete(1).All.Plate==3 & Complete(1).All.Perc_AA==15));
g(1,2).facet_grid(Complete(1).All.Perc_AA,[]);
g(1,2).stat_summary('type','std');
g(1,2).set_title('Flow Count Data');
g(1,2).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair','row','LW%','column',[]);
g(1,2).set_order_options('row',[15 10 5 0],...
    'linestyle',{'crossY','None'},'size',{'None','crossY'});
g(1,2).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,2).set_color_options('map',[1 0 0],'legend','expand');
g(1,2).no_legend();

g(1,2).update('y',Complete(1).All.Ynormalized,'linestyle',Complete(1).All.BactFcn,...
    'size',Complete(1).All.BactFcn,'subset',(Complete(1).All.YeastFcn=="crossY" & ...
    Complete(1).All.BactFcn~="wtB" & Complete(1).All.Assay=="Dynamics" &...
    Complete(1).All.day~=0 & Complete(1).All.Plate==3 & Complete(1).All.Perc_AA==15));
g(1,2).stat_summary('type','std');
g(1,2).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair');
g(1,2).set_order_options('linestyle',{'crossB','None'},'size',{'None','crossB'})
g(1,2).set_color_options('map',[.65 .5 0],'legend','expand');
g(1,2).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,2).no_legend();


g.draw();



%% F2A: Batch culture summary "heatmap"

% This works by manually generating an RGB image based on growth values

Crgb = [];

Pairs = {'CrossB_CrossY','CrossB_wtY','wtB_CrossY','wtB_wtY'};
Percs = [15 10 5 0];

for lw = 1:4
    perc = Percs(lw);
    Crgb1 = [];

    for p = 1:4
        pair = Pairs{p};
        

        Bcount = nanmean(Complete(1).Plate3.Bnormalized(Complete(1).Plate3.day==6 &...
            Complete(1).Plate3.BactYeastPair==pair & Complete(1).Plate3.LW_Percent==perc));
        Ycount = nanmean(Complete(1).Plate3.Ynormalized(Complete(1).Plate3.day==6 &...
            Complete(1).Plate3.BactYeastPair==pair & Complete(1).Plate3.LW_Percent==perc));
        
        if perc==15 && p==4
            Bsave = Bcount;
            Ysave = Ycount;
        end

        box = ones(50^2);
        
        % Make all values in upper triangle of a box the mean of bacterial
        % values. The lower corner becomes zeros. The multiplier of 4 is
        % just to make colormap brighter
        bact = triu(box,1) .* Bcount .* 4; 
        
        % Do the same thing for yeast with lower corner of box 
        yeast = tril(box,-1) .* Ycount .* 4;
        
        % Replicate boxes to R, G, B coordinates
        Brgb = repmat(bact,1,1,3);
        Brgb(:,:,2:3) = 0; % Set green, blue channels to zilch to make red

        Yrgb = repmat(yeast,1,1,3);
        Yrgb(:,:,3) = 0; % Set blue channel to zilch to make yellow

        % Combine the bact & yeast boxes (thanks to those zeros)
        Crgb0 = Yrgb + Brgb;
        
        % Add columns for each pairing
        Crgb1 = horzcat(Crgb1,Crgb0);
    end
    
    % Save out that %LW set as cell
    Crgb{lw} = Crgb1;

end

% Stack the %LW rows 
Crgb = vertcat(Crgb{:});

imshow(Crgb);

%% F2A Legends


Rleg = zeros(100,3);
Yleg = Rleg;

Rleg(1:25,1) = linspace(0,1,25);
Rleg(26:100,1) = 1;

Yleg(1:25,1) = linspace(0,1,25);
Yleg(26:100,1) = 1;
Yleg(:,2) = Yleg(:,1);

figure();
hr = heatmap([0 0.1 0.25; 1 0.08 0],'Colormap',Rleg);

figure();
hy = heatmap([0 0.1 0.25; 1 0.08 0],'Colormap',Yleg);

%% F2B: Batch D:R ratios and TKC

C.Plate1 = Complete(1).Plate3;


clear g;
figure('position',[100 100 600 300]);
    
g(1,1) = gramm('x',C.Plate1.day,'y',C.Plate1.BYnormRatio,...
    'color',C.Plate1.BactYeastPair,'subset',(C.Plate1.YeastFcn~="None" | C.Plate1.BactFcn~="None") &...
    C.Plate1.LW_Percent==10);
g(1,1).stat_summary('geom','line');
g(1,1).geom_point();
g(1,1).axe_property('Yscale','log','Ylim',[10^-5 10^4]);
g(1,1).set_title('Normalized donor-to-recipient ratio');
g(1,1).set_names('x','Day','y',["Normalized bacterial count /","Normalized yeast count"],...
    'color','Coculture Pair','row','LW%','column',[]);
g(1,1).set_order_options('color',{'CrossB_CrossY','CrossB_wtY','wtB_CrossY','wtB_wtY'});

g(1,2) = gramm('x',C.Plate1.day,'y',C.Plate1.TKC,...
    'color',C.Plate1.BactYeastPair,'subset',(C.Plate1.YeastFcn~="None" | C.Plate1.BactFcn~="None") &...
    C.Plate1.LW_Percent==10);
g(1,2).stat_summary('geom','line');
g(1,2).geom_point('dodge',1);
g(1,2).axe_property('Yscale','log');
g(1,2).set_title('Raw IDC Counts');
g(1,2).set_names('x','Day','y','CFU','color','Coculture Pair','row','LW%','column',[]);
g(1,2).set_order_options('color',{'CrossB_CrossY','CrossB_wtY','wtB_CrossY','wtB_wtY'});

% g.set_title('D:R Ratio and TKC Counts');
g.set_point_options('base_size',3);

g.draw();


%% F2C: Batch log-log fits, D:R vs. TKC

clear g;
figure('position',[100 100 400 300]);

g(1,1) = gramm('x',Complete(1).All.logBYnormRatio,'y',Complete(1).All.logTKC,'color',Complete(1).All.BactYeastPair,...
    'subset',(Complete(1).All.YeastFcn~="None" & Complete(1).All.BactFcn~="None" &...
    Complete(1).All.day==7 & Complete(1).All.Assay=="Dynamics" & Complete(1).All.TKC>=2 &...
    Complete(1).All.Plate==3));
g(1,1).geom_point();
% g(1,1).set_title('TKC by D:R ratio, day 7');
g(1,1).set_names('x','ln(D:R)','y','ln(IDC)','column','','row','Day ','color','Coculture Pair');
g(1,1).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

g(1,1).update('color',Complete(1).All.Plate);
g(1,1).stat_glm('disp_fit','true');


g.draw();

for i = 1:size(g.results.stat_glm,1)
    g.results.stat_glm(i).text_handle.Position = [.4 .9];
    g.results.stat_glm(i).text_handle.FontSize = 12;
end


%% F2D: Batch log-log fit slopes

clear g;
figure();

g(1,1) = gramm('x',Complete(1).All.logBYnormRatio,'y',Complete(1).All.logTKC,'color',Complete(1).All.BactYeastPair,...
    'subset',(Complete(1).All.YeastFcn~="None" & Complete(1).All.BactFcn~="None" &...
    Complete(1).All.Plate==3 & Complete(1).All.TKC>=2));
g(1,1).facet_grid(Complete(1).All.day,[]);
g(1,1).geom_point();
% g(1,1).set_title('IDC by D:R ratio, day 7');
g(1,1).set_names('x','ln(D:R)','y','ln(IDC)','row','Day ');
g(1,1).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

g(1,1).update('color',Complete(1).All.Plate);
g(1,1).stat_glm('disp_fit','true');


g.draw();

Slopes = zeros(size(g.results.stat_glm,1),4);

for i = 1:size(g.results.stat_glm,1)
    Slopes(i) = g.results.stat_glm(i).model.Coefficients.Estimate('x1');
    Slopes(i,2) = g.results.stat_glm(i).model.Coefficients.Estimate('x1') -...
        g.results.stat_glm(i).model.Coefficients.SE('x1');
    Slopes(i,3) = g.results.stat_glm(i).model.Coefficients.Estimate('x1') +...
        g.results.stat_glm(i).model.Coefficients.SE('x1');
    Slopes(i,4) = i+1;
end

clear g;
figure('position',[0 0 300 300]);

g(1,1) = gramm('x',Slopes(:,4),'y',Slopes(:,1),'ymin',Slopes(:,2),'ymax',Slopes(:,3));
g(1,1).stat_summary('geom','bar');
g(1,1).geom_interval('geom','black_errorbar','dodge',0.8,'width',0);
% g(1,1).set_title('Fit slopes over time, ln(IDC) vs. ln(D:R)');
g(1,1).set_names('x','Day','y','Slope');


g.draw();


%% F3B: Clumps TKC

clear g;
figure('position',[100 100 600 300]);

EndPoints.Plate23 = Complete(1).Plate16.TecanData;

maxCh = max(EndPoints.Plate23.Cherry(EndPoints.Plate23.day>0));
maxCit = max(EndPoints.Plate23.Citrine(EndPoints.Plate23.day>0));

EndPoints.Plate23.ChNorm = EndPoints.Plate23.Cherry ./ maxCh;
EndPoints.Plate23.CitNorm = EndPoints.Plate23.Citrine ./ maxCit;
EndPoints.Plate23.BYnormRatio = EndPoints.Plate23.ChNorm ./ EndPoints.Plate23.CitNorm;

EndPoints.Plate23.TKCPerRatio = EndPoints.Plate23.TKC ./ EndPoints.Plate23.BYnormRatio;
EndPoints.Plate23.TKCPerCitrine = EndPoints.Plate23.TKC ./ EndPoints.Plate23.Citrine;

g(1,1) = gramm('x',EndPoints.Plate23.day,'y',EndPoints.Plate23.TKC,...
    'color',EndPoints.Plate23.BactYeastPair,...
    'subset',EndPoints.Plate23.PercLW==15 & (EndPoints.Plate23.YeastFcn~="None" & EndPoints.Plate23.BactFcn~="None"));
g(1,1).facet_grid([],EndPoints.Plate23.Mannose);
g(1,1).stat_summary('geom','line');
g(1,1).geom_point('dodge',1);
% g(1,2).axe_property('Yscale','log');
% g(1,1).set_title('Raw TKC Counts, 15% LW');
g(1,1).set_names('x','Day','y','CFU','color','Coculture Pair','column','Mannose');
g(1,1).set_order_options('column',{'N','Y'},...
    'color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

% g.set_title('Mannose Plates TKC');

g.draw();

%% F3B: Clumps TKC ALT for UWSP

clear g;
figure('position',[100 100 600 300]);

EndPoints.Plate23 = Complete(1).Plate16.TecanData;

maxCh = max(EndPoints.Plate23.Cherry(EndPoints.Plate23.day>0));
maxCit = max(EndPoints.Plate23.Citrine(EndPoints.Plate23.day>0));

EndPoints.Plate23.ChNorm = EndPoints.Plate23.Cherry ./ maxCh;
EndPoints.Plate23.CitNorm = EndPoints.Plate23.Citrine ./ maxCit;
EndPoints.Plate23.BYnormRatio = EndPoints.Plate23.ChNorm ./ EndPoints.Plate23.CitNorm;

EndPoints.Plate23.TKCPerRatio = EndPoints.Plate23.TKC ./ EndPoints.Plate23.BYnormRatio;
EndPoints.Plate23.TKCPerCitrine = EndPoints.Plate23.TKC ./ EndPoints.Plate23.Citrine;

g(1,1) = gramm('x',EndPoints.Plate23.day,'y',EndPoints.Plate23.TKC,...
    'color',EndPoints.Plate23.BactYeastPair,...
    'subset',EndPoints.Plate23.PercLW==15 & (EndPoints.Plate23.YeastFcn=="crossY" &...
    EndPoints.Plate23.BactFcn=="crossB"));
g(1,1).facet_grid([],EndPoints.Plate23.Mannose);
g(1,1).stat_summary('geom','line');
g(1,1).geom_point('dodge',1);
% g(1,2).axe_property('Yscale','log');
g(1,1).set_title('Raw TKC Counts, 15% LW');
g(1,1).set_names('x','Day','y','Yeasts with DNA','color','Coculture Pair','column','Mannose');
g(1,1).set_order_options('column',{'N','Y'},...
    'color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

g.set_title('Conjugated yeast, clumped vs. un-clumped cultures');

g.draw();
%% F3C: Clumps D:R


clear g;
h = figure('position',[0 0 500 400]);

g(1,1) = gramm('x',Complete(1).Plate16.TecanData.PercLW,'y',Complete(1).Plate16.TecanData.NormBYRatio,...
    'ymin',Complete(1).Plate16.TecanData.CImin,'ymax',Complete(1).Plate16.TecanData.CImax,...
    'color',Complete(1).Plate16.TecanData.Mannose,...
    'subset',(Complete(1).Plate16.TecanData.BactYeastPair=="crossB_wtY" & Complete(1).Plate16.TecanData.day==6));
g(1,1).stat_summary('geom','bar','width',0.5);
g(1,1).geom_interval('geom','black_errorbar','dodge',.6,'width',0);
g(1,1).axe_property('YLim',[0 1]);
% g(1,1).set_title(["D:R ratio (normalized Cherry/Citrine)","crosB_wtY, day 6"]);
g(1,1).set_names('x','LW%','y',["Normalized mCherry /","Normalized ymCitrine"],'color','Mannose');
g(1,1).set_order_options('x',[15 10 5 0],'color',{'Y','N'});

% g.set_title('Mannose Plates D:R');

g.draw();

ZoomHandle = zoom(h);
set(ZoomHandle,'Motion','horizontal')

p = figure('position',[0 400 500 400]);

g(1,1) = gramm('x',Complete(1).Plate16.TecanData.PercLW,'y',Complete(1).Plate16.TecanData.NormBYRatio,...
    'ymin',Complete(1).Plate16.TecanData.CImin,'ymax',Complete(1).Plate16.TecanData.CImax,...
    'color',Complete(1).Plate16.TecanData.Mannose,...
    'subset',(Complete(1).Plate16.TecanData.BactYeastPair=="crossB_wtY" & Complete(1).Plate16.TecanData.day==6));
g(1,1).stat_summary('geom','bar','width',0.5);
g(1,1).geom_interval('geom','black_errorbar','dodge',.6,'width',0);
g(1,1).axe_property('YLim',[0 15]);
% g(1,1).set_title(["D:R ratio (normalized Cherry/Citrine)","crosB_wtY, day 6"]);
g(1,1).set_names('x','LW%','y',["Normalized mCherry /","Normalized ymCitrine"],'color','Mannose');
g(1,1).set_order_options('x',[15 10 5 0],'color',{'Y','N'});

% g.set_title('Mannose Plates D:R');

g.draw();

ZoomHandle = zoom(p);
set(ZoomHandle,'Motion','horizontal')


%% F5C: Colony Colocalization

data = Complete(1).Plate11.ImageData;
data.logTKC = log(data.TKC);
data.logTKC(data.logTKC==Inf) = 0;
data.logTKC(data.logTKC==-Inf) = 0;

ref = nanmean(data.ICQ(data.Day==6 & (data.Bacteria~='None' | data.Yeast~='None') & data.Type=='Cis' &...
    data.Perc_LW==100));

clear g;
figure();

g(1,1) = gramm('x',data.Day,'y',data.ICQ,'ymin',data.CImin,'ymax',data.CImax,'color',data.BactYeastPair,...
    'subset',(data.Day>0 & (data.Bacteria~='None' | data.Yeast~='None') & data.Type=='Cis'));
g(1,1).facet_grid(data.Perc_LW,[]);
g(1,1).geom_point('dodge',0.6);
g(1,1).geom_interval('geom','errorbar','dodge',0.6,'width',0);
g(1,1).geom_hline('yintercept',ref);
% g(1,1).set_title('Colocalization over time, cis');
g(1,1).set_names('x','Day','y','ICQ (colocalization)','color','Coculture Pair','row','LW%','column',[]);
g(1,1).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'},'row',[100 0]);

% g.set_title('Li colocalization');

g.draw();


%% F5D: Colony Colocalization vs. TKC

data = Complete(1).Plate11.ImageData;
data.logTKC = log(data.TKC);
data.logTKC(data.logTKC==Inf) = 0;
data.logTKC(data.logTKC==-Inf) = 0;

clear g;
figure('position',[100 100 400 300]);

g(1,1) = gramm('x',data.ICQ,'y',data.TKC,'color',data.BactYeastPair,...
    'subset',(~isnan(data.TKC) & data.Type=='Cis' & data.logTKC>=0));
g(1,1).facet_grid(data.Perc_LW,[]);
g(1,1).geom_point();
% g(1,1).set_title('IDC per Colocalization, cis');
g(1,1).set_names('x','ICQ (colocalization)','y','IDC','color','Coculture Pair','row','LW%','column',[]);
g(1,1).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'},'row',[100 0]);

g(1,1).update('color',data.Perc_LW);
g(1,1).stat_glm('disp_fit','true');

% g.set_title('Comparing Li colocalization to TKC');

g.draw();

for i = 1:size(g(1,1).results.stat_glm,1)
    rsquare = g(1,1).results.stat_glm(i).model.Rsquared.Ordinary;
    g(1,1).results.stat_glm(i).text_handle.String = ['R^2 = ' num2str(rsquare)];
    
    g(1,1).results.stat_glm(i).text_handle.Position = [0.03 .7];
    g(1,1).results.stat_glm(i).text_handle.FontSize = 10;
%     text(-.2,2,['R^2 = ' num2str(rsquare)]);


end


%% F6A: Rescue growth (UH=10, crossY)

clear g;
figure('position',[100 100 400 300]);

g(1,1) = gramm('x',Complete(1).AllRescues.day,'y',Complete(1).AllRescues.Ynormalized,...
    'linestyle',Complete(1).AllRescues.BactFcn,...
    'subset',(Complete(1).AllRescues.YeastFcn=="crossY" & ...
    Complete(1).AllRescues.Perc_UH==10));
g(1,1).stat_summary('type','std');
% g(1,1).set_title('Recipent Yeast Growth');
g(1,1).set_names('x','Day','y','Normalized yeast count','linestyle','Coculture Pair');
g(1,1).set_order_options('linestyle',{'crossB','wtB','None'});
g(1,1).set_line_options('base_size',2,'styles',{'-','-.',':'});
g(1,1).set_color_options('map',[.65 .5 0],'legend','expand');
g(1,1).no_legend();

g.draw();


%% F6B: Rescue culture summary "heatmap"

% This works by manually generating an RGB image based on growth values

Crgb = [];

Pairs = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'};
Percs = [15 10 5 0];

for uh = 1:4
    perc = Percs(uh);
    Crgb1 = [];

    for p = 1:4
        pair = Pairs{p};
        

        Bcount = nanmean(Complete(1).AllRescues.Bnormalized(Complete(1).AllRescues.day==6 &...
            Complete(1).AllRescues.BactYeastPair==pair & Complete(1).AllRescues.Perc_UH==perc));
        Ycount = nanmean(Complete(1).AllRescues.Ynormalized(Complete(1).AllRescues.day==6 &...
            Complete(1).AllRescues.BactYeastPair==pair & Complete(1).AllRescues.Perc_UH==perc));


        box = ones(50^2);
        
        % Make all values in upper triangle of a box the mean of bacterial
        % values. The lower corner becomes zeros. The multiplier of 4 is
        % just to make colormap brighter
        bact = triu(box,1) .* Bcount .* 4; 
        
        % Do the same thing for yeast with lower corner of box 
        yeast = tril(box,-1) .* Ycount .* 4;
        
        % Replicate boxes to R, G, B coordinates
        Brgb = repmat(bact,1,1,3);
        Brgb(:,:,2:3) = 0; % Set green, blue channels to zilch to make red

        Yrgb = repmat(yeast,1,1,3);
        Yrgb(:,:,3) = 0; % Set blue channel to zilch to make yellow

        % Combine the bact & yeast boxes (thanks to those zeros)
        Crgb0 = Yrgb + Brgb;
        
        % Add columns for each pairing
        Crgb1 = horzcat(Crgb1,Crgb0);
    end
    
    % Save out that %LW set as cell
    Crgb{uh} = Crgb1;

end

% Stack the %UH rows 
Crgb = vertcat(Crgb{:});

imshow(Crgb);


%% F7B: Growth (Yeast Only)

AllData.Plate1 = Complete(1).Plate15.TecanEndpoints;

maxCh = max(AllData.Plate1.Cherry(AllData.Plate1.day>1));
maxCit = max(AllData.Plate1.Citrine(AllData.Plate1.day>1));

clear g;
figure('position',[100 100 400 600]);
g(1,1) = gramm('x',AllData.Plate1.time./3600,'y',AllData.Plate1.Citrine./maxCit,...
   'linestyle',AllData.Plate1.BactFcn,...
   'subset',AllData.Plate1.YeastFcn=='crisprY' & AllData.Plate1.PercU==0);
g(1,1).facet_grid(AllData.Plate1.PercW,[]);
g(1,1).stat_summary('type','std');
g(1,1).geom_vline('xintercept',509335/3600);
% g(1,1).set_title('Cutter donor');
g(1,1).axe_property('YLim',[0 1]);
g(1,1).set_names('x','t (hr)','y',["Normalized", "Yeast Fluorescence"],'row','W%','linestyle','Yeast Pairing');
g(1,1).set_order_options('row',[15 10 5 1],'linestyle',{'cutB','nocutB','None'});
g(1,1).set_line_options('base_size',2,'step_size',2,'styles',{'-','--',':'});
g(1,1).set_color_options('map',[.7 .6 0],'legend','expand');
g(1,1).no_legend();

% g.set_title('CRISPR Yeast Growth (normalized)');

g.draw();

%% F7B: Growth (Yeast Only) alt orientation

AllData.Plate1 = Complete(1).Plate15.TecanEndpoints;

maxCh = max(AllData.Plate1.Cherry(AllData.Plate1.day>1));
maxCit = max(AllData.Plate1.Citrine(AllData.Plate1.day>1));

clear g;
figure('position',[100 100 800 250]);
g(1,1) = gramm('x',AllData.Plate1.time./3600,'y',AllData.Plate1.Citrine./maxCit,...
   'linestyle',AllData.Plate1.BactFcn,...
   'subset',AllData.Plate1.YeastFcn=='crisprY' & AllData.Plate1.PercU==0);
g(1,1).facet_grid([],AllData.Plate1.PercW);
g(1,1).stat_summary('type','std');
g(1,1).geom_vline('xintercept',509335/3600);
% g(1,1).set_title('Cutter donor');
g(1,1).axe_property('YLim',[0 1]);
g(1,1).set_names('x','t (hr)','y',["Normalized", "Yeast Fluorescence"],'column','W%','linestyle','Yeast Pairing');
g(1,1).set_order_options('column',[15 10 5 1],'linestyle',{'cutB','nocutB','None'});
g(1,1).set_line_options('base_size',2,'step_size',2,'styles',{'-','--',':'});
g(1,1).set_color_options('map',[.7 .6 0],'legend','expand');
g(1,1).no_legend();

% g.set_title('CRISPR Yeast Growth (normalized)');

g.draw();

%% F7B: Growth (Yeast Only) alt orientation, no negs

AllData.Plate1 = Complete(1).Plate15.TecanEndpoints;

maxCh = max(AllData.Plate1.Cherry(AllData.Plate1.day>1));
maxCit = max(AllData.Plate1.Citrine(AllData.Plate1.day>1));

clear g;
figure('position',[100 100 800 250]);
g(1,1) = gramm('x',AllData.Plate1.day,'y',AllData.Plate1.Citrine./maxCit,...
   'linestyle',AllData.Plate1.BactFcn,...
   'subset',AllData.Plate1.YeastFcn=='crisprY' & AllData.Plate1.BactFcn~='nocutB' &...
   AllData.Plate1.PercU==0);
g(1,1).facet_grid([],AllData.Plate1.PercW);
g(1,1).stat_summary('type','std');
% g(1,1).geom_vline('xintercept',509335/3600);
% g(1,1).set_title('Cutter donor');
g(1,1).axe_property('YLim',[0 1]);
g(1,1).set_names('x','Day','y','Normalized Yeast Count','column','W%','linestyle','Yeast Pairing');
g(1,1).set_order_options('column',[15 10 5 1],'linestyle',{'cutB','nocutB','None'});
g(1,1).set_line_options('base_size',2,'step_size',2,'styles',{'-','--',':'});
g(1,1).set_color_options('map',[.7 .6 0],'legend','expand');
g(1,1).no_legend();

% g.set_title('CRISPR Yeast Growth (normalized)');

g.draw();

%% F7B: Growth (Yeast Only) alt for UWSP

AllData.Plate1 = Complete(1).Plate15.TecanEndpoints;

maxCh = max(AllData.Plate1.Cherry(AllData.Plate1.day>1));
maxCit = max(AllData.Plate1.Citrine(AllData.Plate1.day>1));

clear g;
figure('position',[100 100 800 250]);
g(1,1) = gramm('x',AllData.Plate1.day,'y',AllData.Plate1.Citrine./maxCit,...
   'linestyle',AllData.Plate1.BactFcn,...
   'subset',AllData.Plate1.YeastFcn=='crisprY' & AllData.Plate1.BactFcn~='nocutB' &...
   AllData.Plate1.PercU==0);
g(1,1).facet_grid([],AllData.Plate1.PercW);
g(1,1).stat_summary('type','std');
% g(1,1).geom_vline('xintercept',509335/3600);
% g(1,1).set_title('Cutter donor');
g(1,1).axe_property('YLim',[0 1]);
g(1,1).set_names('x','Day','y','Normalized Yeast Count','linestyle','Yeast Pairing');
g(1,1).set_order_options('column',[10 5],'linestyle',{'cutB','nocutB','None'});
g(1,1).set_line_options('base_size',2,'step_size',2,'styles',{'-','--',':'});
g(1,1).set_color_options('map',[.7 .6 0],'legend','expand');
g(1,1).no_legend();

% g.set_title('CRISPR Yeast Growth (normalized)');

g.draw();

%% F7C: CRISPR Fold depression Yeast
% 
% 
% EndPoints.Plate1 = Complete(1).Plate15.TecanEndpoints;
% 
% clear g;
% figure('position',[100 100 800 400]);
% 
% g(1,1) = gramm('x',EndPoints.Plate1.day,'y',EndPoints.Plate1.FoldDecreaseMono,...
%     'color',EndPoints.Plate1.BactYeastPair,...
%     'subset',(EndPoints.Plate1.YeastFcn~="None" & EndPoints.Plate1.PercU==0));
% g(1,1).facet_grid(EndPoints.Plate1.PercW,[]);
% g(1,1).stat_summary('geom','bar');
% g(1,1).geom_vline('xintercept',6.5);
% g(1,1).axe_property('Ylim',[-8 -.9]);
% g(1,1).set_title('Change from monoculture, 0% Ura');
% g(1,1).set_names('x','Day','y',["Fold decrease","(Citrine)"],'color','Coculture Pair','row','W%','column','U%');
% g(1,1).set_order_options('row',[15 10 5 1],...
%     'color',{'cutB_crisprY','nocutB_crisprY'});
% 
% g(1,2) = gramm('x',EndPoints.Plate1.day,'y',EndPoints.Plate1.FoldDecreaseCoCo,...
%     'color',EndPoints.Plate1.BactYeastPair,...
%     'subset',(EndPoints.Plate1.BactYeastPair=="cutB_crisprY" & EndPoints.Plate1.YeastFcn~="None" & EndPoints.Plate1.PercU==0));
% g(1,2).facet_grid(EndPoints.Plate1.PercW,EndPoints.Plate1.PercU);
% g(1,2).stat_summary('geom','bar');
% g(1,2).geom_vline('xintercept',6.5);
% g(1,2).axe_property('Ylim',[-8 -.9]);
% g(1,2).set_title('Change from coculture negative, 0% Ura');
% g(1,2).set_names('x','Day','y',["Fold decrease","(Citrine)"],'color','Coculture Pair','row','W%','column','U%');
% g(1,2).set_order_options('row',[15 10 5 1]);
% 
% g.set_title('Comparing Fold-Loss of Yeast from Cocultures');
% 
% g.draw();

%% F7C: CRISPR Fold depression Yeast


EndPoints.Plate1 = Complete(1).Plate15.TecanEndpoints;

clear g;
figure('position',[100 100 500 400]);

g(1,1) = gramm('x',EndPoints.Plate1.day,'y',EndPoints.Plate1.FoldDecreaseMono,...
    'ymin',EndPoints.Plate1.FDM_CImin,'ymax',EndPoints.Plate1.FDM_CImax,...
    'color',EndPoints.Plate1.BactYeastPair,...
    'subset',(EndPoints.Plate1.YeastFcn~="None" & EndPoints.Plate1.PercU==0));
g(1,1).facet_grid(EndPoints.Plate1.PercW,[]);
g(1,1).stat_summary('geom','point','dodge',0.4);
g(1,1).geom_interval('geom','errorbar','dodge',0.4,'width',0);
g(1,1).geom_vline('xintercept',6.5);
g(1,1).axe_property('Ylim',[-8 .9]);
% g(1,1).set_title('Change from monoculture, 0% Ura');
g(1,1).set_names('x','Day','y',["Fold decrease","(Citrine)"],'color','Coculture Pair','row','W%','column','U%');
g(1,1).set_order_options('row',[15 10 5 1],...
    'color',{'cutB_crisprY','nocutB_crisprY'});

% g.set_title('Comparing Fold-Loss of Yeast from Monoculture');

g.draw();

%% F7C alt: CRISPR Fold depression Yeast


EndPoints.Plate1 = Complete(1).Plate15.TecanEndpoints;

clear g;
figure('position',[100 100 500 400]);

g(1,1) = gramm('x',EndPoints.Plate1.day,'y',EndPoints.Plate1.FoldDecreaseMono,...
    'ymin',EndPoints.Plate1.FDM_CImin,'ymax',EndPoints.Plate1.FDM_CImax,...
    'color',EndPoints.Plate1.BactYeastPair,...
    'subset',(EndPoints.Plate1.YeastFcn~="None" & EndPoints.Plate1.BactFcn~='nocutB' &...
    EndPoints.Plate1.PercU==0));
g(1,1).facet_grid(EndPoints.Plate1.PercW,[]);
g(1,1).stat_summary('geom','bar','dodge',0.4);
g(1,1).geom_interval('geom','errorbar','dodge',0.4,'width',0,'color','b');
g(1,1).geom_vline('xintercept',6.5);
g(1,1).axe_property('Ylim',[-8 .9]);
g(1,1).set_title('Change from monoculture, 0% Ura');
g(1,1).set_names('x','Day','y',["Fold decrease","(Citrine)"],'color','Coculture Pair','row','W%','column','U%');
g(1,1).set_order_options('row',[15 10 5 1],...
    'color',{'cutB_crisprY','nocutB_crisprY'});

g.set_title('Comparing Fold-Loss of Yeast from Monoculture');

g.draw();

%% F7D: CRISPR TKC

EndPoints.Plate1 = Complete(1).Plate15.TecanEndpoints;

clear g;
figure('position',[100 100 400 300]);

g(1,1) = gramm('x',EndPoints.Plate1.day,'y',EndPoints.Plate1.TKC,...
    'color',EndPoints.Plate1.BactYeastPair,'subset',(EndPoints.Plate1.YeastFcn~="None" &...
    EndPoints.Plate1.BactFcn~="None" & EndPoints.Plate1.PercU==0));
g(1,1).facet_grid(EndPoints.Plate1.PercW,[]);
g(1,1).stat_summary('geom','line');
g(1,1).geom_point();
g(1,1).geom_vline('xintercept',6.5);
% g(1,1).set_title('Raw TKC Counts');
g(1,1).set_names('x','Day','y','CFU','color','Coculture Pair','row','W%');
g(1,1).set_order_options('row',[15 10 5 1],...
    'color',{'cutB_crisprY','nocutB_crisprY'});


% g.set_title('CRISPR TKC');

g.draw();

%% F7D: CRISPR TKC alt

EndPoints.Plate1 = Complete(1).Plate15.TecanEndpoints;

clear g;
figure('position',[100 100 400 300]);

g(1,1) = gramm('x',EndPoints.Plate1.day,'y',EndPoints.Plate1.TKC,...
    'color',EndPoints.Plate1.BactYeastPair,'subset',(EndPoints.Plate1.YeastFcn~="None" &...
    EndPoints.Plate1.BactFcn=="cutB" & EndPoints.Plate1.PercU==0 & EndPoints.Plate1.day>2 & EndPoints.Plate1.PercW==10));
% g(1,1).facet_grid(EndPoints.Plate1.PercW,[]);
g(1,1).stat_summary('geom','line');
g(1,1).geom_point();
g(1,1).geom_vline('xintercept',6.5);
% g(1,1).set_title('Raw TKC Counts');
g(1,1).set_names('x','Day','y','IDC count','color','Coculture Pair','row','W%');
g(1,1).set_order_options('row',[15 10 5 1],...
    'color',{'cutB_crisprY','nocutB_crisprY'});


% g.set_title('CRISPR TKC');

g.draw();


%% F7D: CRISPR IDC, updated samples 240320

clear g;

C.Plate1 = Complete(1).Plate17;

g(1,1) = gramm('x',C.Plate1.day,'y',C.Plate1.TKC,...
    'color',C.Plate1.Smpl_Grp,'linestyle',C.Plate1.Mannose,'subset',(C.Plate1.YeastFcn~="None" &...
    C.Plate1.BactFcn=="cutB"));
g(1,1).facet_grid(C.Plate1.PercW,[]);
% g(1,1).stat_summary();
g(1,1).stat_summary('geom','line');
g(1,1).geom_point();
g(1,1).geom_vline('xintercept',6);
g(1,1).set_names('x','Day','y','CFU','color','Sample Group','linestyle','Mannose','row','W%');
g(1,1).set_order_options('row',[15 10 5 0],...
    'color',{'A','B'},'linestyle',{'N','Y'});
g(1,1).set_line_options('base_size',2,'step_size',2,'styles',{'-','--'});


% g.set_title('IDC from cutter bacteria');

g.draw();


%% END OF BODY FIGURES













% Next up, SI figs...


















%% SF2: Batch culture growths, cis

clear g;
figure();

g(1,1) = gramm('x',Complete(1).All.day,'y',Complete(1).All.Bnormalized,...
    'linestyle',Complete(1).All.YeastFcn,'size',Complete(1).All.YeastFcn,...
    'subset',(Complete(1).All.BactFcn=="crossB" & Complete(1).All.day~=0 &...
    Complete(1).All.YeastFcn~="wtY" & Complete(1).All.Assay=="Dynamics" &...
     Complete(1).All.Plate==3));
g(1,1).facet_grid(Complete(1).All.Perc_AA,[]);
g(1,1).stat_summary('type','std');
% g(1,1).set_title('crossB_crossB');
g(1,1).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair','row','LW%','column',[]);
g(1,1).set_order_options('row',[15 10 5 0],...
    'linestyle',{'crossY','None'},'size',{'None','crossY'});
g(1,1).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,1).set_color_options('map',[1 0 0],'legend','expand');
g(1,1).no_legend();

g(1,1).update('y',Complete(1).All.Ynormalized,'linestyle',Complete(1).All.BactFcn,...
    'size',Complete(1).All.BactFcn,'subset',(Complete(1).All.YeastFcn=="crossY" & ...
    Complete(1).All.BactFcn~="wtB" & Complete(1).All.Assay=="Dynamics" &...
    Complete(1).All.day~=0 & Complete(1).All.Plate==3));
g(1,1).stat_summary('type','std');
g(1,1).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair');
g(1,1).set_order_options('linestyle',{'crossB','None'},'size',{'None','crossB'},...
    'row',[15 10 5 0])
g(1,1).set_color_options('map',[.65 .5 0],'legend','expand');
g(1,1).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,1).no_legend();

g(1,2) = gramm('x',Complete(1).All.day,'y',Complete(1).All.Bnormalized,...
    'linestyle',Complete(1).All.YeastFcn,'size',Complete(1).All.YeastFcn,...
    'subset',(Complete(1).All.BactFcn=="crossB" &...
    Complete(1).All.YeastFcn~="crossY" & Complete(1).All.Assay=="Dynamics" &...
    Complete(1).All.day~=0 & Complete(1).All.Plate==3));
g(1,2).facet_grid(Complete(1).All.Perc_AA,[]);
g(1,2).stat_summary('type','std');
% g(1,2).set_title('crossB_wtY');
g(1,2).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair','row','LW%','column',[]);
g(1,2).set_order_options('row',[15 10 5 0],'size',{'None','wtY'},...
    'linestyle',{'wtY','None'});
g(1,2).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,2).set_color_options('map',[1 0 0],'legend','expand');
g(1,2).no_legend();

g(1,2).update('y',Complete(1).All.Ynormalized,'linestyle',Complete(1).All.BactFcn,...
    'size',Complete(1).All.BactFcn,'subset',(Complete(1).All.YeastFcn=="wtY" & ...
    Complete(1).All.BactFcn~="wtB" & Complete(1).All.Assay=="Dynamics" &...
    Complete(1).All.day~=0 & Complete(1).All.Plate==3));
g(1,2).stat_summary('type','std');
g(1,2).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair');
g(1,2).set_order_options('row',[15 10 5 0],'linestyle',{'crossB','None'},'size',{'None','crossB'});
g(1,2).set_color_options('map',[.65 .5 0],'legend','expand');
g(1,2).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,2).no_legend();

g(1,3) = gramm('x',Complete(1).All.day,'y',Complete(1).All.Bnormalized,...
    'linestyle',Complete(1).All.YeastFcn,'size',Complete(1).All.YeastFcn,...
    'subset',(Complete(1).All.BactFcn=="wtB") & ...
    Complete(1).All.YeastFcn~="wtY" & Complete(1).All.Assay=="Dynamics" &...
    Complete(1).All.day~=0 & Complete(1).All.Plate==3);
g(1,3).facet_grid(Complete(1).All.Perc_AA,[]);
g(1,3).stat_summary('type','std');
% g(1,3).set_title('wtB_crossY');
g(1,3).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair','row','LW%','column',[]);
g(1,3).set_order_options('row',[15 10 5 0],'size',{'None','crossY'},...
    'linestyle',{'crossY','None'});
g(1,3).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,3).set_color_options('map',[1 0 0],'legend','expand');
g(1,3).no_legend();

g(1,3).update('y',Complete(1).All.Ynormalized,'linestyle',Complete(1).All.BactFcn,...
    'size',Complete(1).All.BactFcn,'subset',(Complete(1).All.YeastFcn=="crossY" &...
    Complete(1).All.BactFcn~="crossB" & Complete(1).All.Assay=="Dynamics" &...
    Complete(1).All.day~=0 & Complete(1).All.Plate==3));
g(1,3).stat_summary('type','std');
g(1,3).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair');
g(1,3).set_order_options('linestyle',{'wtB','None'},'size',{'None','wtB'},'row',[15 10 5 0]);
g(1,3).set_color_options('map',[.65 .5 0],'legend','expand');
g(1,3).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,3).no_legend();

g(1,4) = gramm('x',Complete(1).All.day,'y',Complete(1).All.Bnormalized,'linestyle',Complete(1).All.YeastFcn,...
    'size',Complete(1).All.YeastFcn,'subset',(Complete(1).All.BactFcn=="wtB") & ...
    Complete(1).All.YeastFcn~="crossY" & Complete(1).All.Assay=="Dynamics" &...
    Complete(1).All.day~=0 & Complete(1).All.Plate==3);
g(1,4).facet_grid(Complete(1).All.Perc_AA,[]);
g(1,4).stat_summary('type','std');
% g(1,4).set_title('wtB_wtY');
g(1,4).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair','row','LW%','column',[]);
g(1,4).set_order_options('row',[15 10 5 0],'size',{'None','wtY'},...
    'linestyle',{'wtY','None'});
g(1,4).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,4).set_color_options('map',[1 0 0],'legend','expand');
g(1,4).no_legend();

g(1,4).update('y',Complete(1).All.Ynormalized,'linestyle',Complete(1).All.BactFcn,...
    'size',Complete(1).All.BactFcn,'subset',(Complete(1).All.YeastFcn=="wtY" & ...
    Complete(1).All.BactFcn~="crossB" & Complete(1).All.Assay=="Dynamics" &...
    Complete(1).All.day~=0 & Complete(1).All.Plate==3));
g(1,4).stat_summary('type','std');
g(1,4).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair');
g(1,4).set_order_options('row',[15 10 5 0],'linestyle',{'wtB','None'},'size',{'None','wtB'});
g(1,4).set_color_options('map',[.65 .5 0],'legend','expand');
g(1,4).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,4).no_legend();


% g.set_title('Growth: All Cis Dynamics');

g.draw();

%% SF3: All_DR_TKC

clear g;
figure();

g(1,1) = gramm('x',Complete(1).All.logBYnormRatio,'y',Complete(1).All.logTKC,'color',Complete(1).All.BactYeastPair,...
    'subset',(Complete(1).All.YeastFcn~="None" & Complete(1).All.BactFcn~="None" &...
    Complete(1).All.day~=0 & Complete(1).All.Assay=="Dynamics" &...
    Complete(1).All.TKC>=2 & (Complete(1).All.Plate==3 |Complete(1).All.Plate==7)));
g(1,1).facet_grid(Complete(1).All.day,Complete(1).All.Txfer);
g(1,1).geom_point();
% g(1,1).set_title('Crossfeeding Cultures');
g(1,1).set_names('x','ln(D:R)','y','ln(IDC)','column','Transfer Type:','row','Day ');
g(1,1).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});
g(1,1).no_legend;

g(1,1).update('color',Complete(1).All.day);
g(1,1).stat_glm('disp_fit','true');
g(1,1).no_legend;


g(1,2) = gramm('x',Complete(1).All.logBYnormRatio,'y',Complete(1).All.logTKC,'color',Complete(1).All.BactYeastPair,...
    'subset',(Complete(1).All.YeastFcn~="None" & Complete(1).All.BactFcn~="None" &...
    Complete(1).All.day~=0 & Complete(1).All.Assay=="Colony" &...
    Complete(1).All.TKC>=2));
g(1,2).facet_grid(Complete(1).All.day,Complete(1).All.Txfer);
g(1,2).geom_point();
% g(1,2).set_title('Crossfeeding Colonies');
g(1,2).set_names('x','ln(D:R)','y','ln(IDC)','column','Transfer Type:','row','Day ');
g(1,2).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});
g(1,2).no_legend;

g(1,2).update('color',Complete(1).All.day);
g(1,2).stat_glm('disp_fit','true');
g(1,2).no_legend;


g(1,3) = gramm('x',Complete(1).All.logBYnormRatio,'y',Complete(1).All.logTKC,'color',Complete(1).All.BactYeastPair,...
    'subset',(Complete(1).All.YeastFcn~="None" & Complete(1).All.BactFcn~="None" &...
    Complete(1).All.day~=0 & Complete(1).All.Assay=="Rescue" &...
    Complete(1).All.TKC>=2 & Complete(1).All.TKC~=500));
g(1,3).facet_grid(Complete(1).All.day,[]);
g(1,3).geom_point();
% g(1,3).set_title(["Rescue Cultures","(cis)"]);
g(1,3).set_names('x','ln(D:R)','y','ln(IDC)','row','Day ','color','Coculture Pair');
g(1,3).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

g(1,3).update('color',Complete(1).All.day);
g(1,3).stat_glm('disp_fit','true');
g(1,3).no_legend;

% g.set_title('TKC by Cell Ratio');

g.draw();

for j = 1:size(g,2)
for i = 1:size(g(1,j).results.stat_glm,1)
    g(1,j).results.stat_glm(i).text_handle.Position = [.4 .9];
    g(1,j).results.stat_glm(i).text_handle.FontSize = 10;
end
end

%% SF3: alt 1 for dissertation 

clear g;
figure();

g(1,1) = gramm('x',Complete(1).All.logBYnormRatio,'y',Complete(1).All.logTKC,'color',Complete(1).All.BactYeastPair,...
    'subset',(Complete(1).All.YeastFcn~="None" & Complete(1).All.BactFcn~="None" &...
    Complete(1).All.day~=0 & Complete(1).All.Assay=="Dynamics" &...
    Complete(1).All.TKC>=2 & (Complete(1).All.Plate==3 |Complete(1).All.Plate==7)));
g(1,1).facet_grid(Complete(1).All.day,Complete(1).All.Txfer);
g(1,1).geom_point();
g(1,1).set_title('TKC by D:R, Dynamics Assays');
g(1,1).set_names('x','ln(D:R)','y','ln(TKC)','column','Transfer Type:','row','Day ');
g(1,1).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});
% g(1,1).no_legend;

g(1,1).update('color',Complete(1).All.day);
g(1,1).stat_glm('disp_fit','true');
g(1,1).no_legend;



g.set_title('TKC by Cell Ratio');

g.draw();

for j = 1:size(g,2)
for i = 1:size(g(1,j).results.stat_glm,1)
    g(1,j).results.stat_glm(i).text_handle.Position = [.5 .9];
    g(1,j).results.stat_glm(i).text_handle.FontSize = 12;
end
end

%% SF3: alt 2 for dissertation

clear g;
figure();

g(1,1) = gramm('x',Complete(1).All.logBYnormRatio,'y',Complete(1).All.logTKC,'color',Complete(1).All.BactYeastPair,...
    'subset',(Complete(1).All.YeastFcn~="None" & Complete(1).All.BactFcn~="None" &...
    Complete(1).All.day~=0 & Complete(1).All.Assay=="Colony" &...
    Complete(1).All.TKC>=2));
g(1,1).facet_grid(Complete(1).All.day,Complete(1).All.Txfer);
g(1,1).geom_point();
g(1,1).set_title('TKC by D:R, Colony Assays');
g(1,1).set_names('x','ln(D:R)','y','ln(TKC)','column','Transfer Type:','row','Day ');
g(1,1).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});
% g(1,1).no_legend;

g(1,1).update('color',Complete(1).All.day);
g(1,1).stat_glm('disp_fit','true');
g(1,1).no_legend;



g.set_title('TKC by Cell Ratio');

g.draw();

for j = 1:size(g,2)
for i = 1:size(g(1,j).results.stat_glm,1)
    g(1,j).results.stat_glm(i).text_handle.Position = [.5 .9];
    g(1,j).results.stat_glm(i).text_handle.FontSize = 12;
end
end

%% SF3: alt 3 for dissertation

clear g;
figure();

g(1,1) = gramm('x',Complete(1).All.logBYnormRatio,'y',Complete(1).All.logTKC,'color',Complete(1).All.BactYeastPair,...
    'subset',(Complete(1).All.YeastFcn~="None" & Complete(1).All.BactFcn~="None" &...
    Complete(1).All.day~=0 & Complete(1).All.Assay=="Rescue" &...
    Complete(1).All.TKC>=2 & Complete(1).All.TKC~=500));
g(1,1).facet_grid(Complete(1).All.day,[]);
g(1,1).geom_point();
g(1,1).set_title('TKC by D:R, Rescue Assays');
g(1,1).set_names('x','ln(D:R)','y','ln(TKC)','row','Day ');
g(1,1).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

g(1,1).update('color',Complete(1).All.day);
g(1,1).stat_glm('disp_fit','true');
g(1,1).no_legend;

g.set_title('TKC by Cell Ratio');

g.draw();

for j = 1:size(g,2)
for i = 1:size(g(1,j).results.stat_glm,1)
    g(1,j).results.stat_glm(i).text_handle.Position = [.5 .9];
    g(1,j).results.stat_glm(i).text_handle.FontSize = 12;
end
end

%% SF4 TKC per Yeast

clear g;
figure();

C.Plate1 = Complete(1).Plate3;
g(1,1) = gramm('x',C.Plate1.day,'y',C.Plate1.TKCPerY,...
    'color',C.Plate1.BactYeastPair,...
    'subset',C.Plate1.day~=0 & (C.Plate1.YeastFcn~="None" & C.Plate1.BactFcn~="None"));
g(1,1).facet_grid(C.Plate1.LW_Percent,[]);
g(1,1).stat_summary('geom','line');
g(1,1).geom_point();
g(1,1).axe_property('Yscale','log','Ylim',[1E-8 5E-2]);
g(1,1).set_title('cis Culture');
g(1,1).set_names('x','Day','y',["Transconjugant CFU/","Yeast Count"],...
    'color','Coculture Pair','row','LW%','column',[]);
g(1,1).set_order_options('row',[15 10 5 0],...
    'color',{'CrossB_CrossY','CrossB_wtY','wtB_CrossY','wtB_wtY'});

C.Plate1 = Complete(1).Plate7;
g(1,2) = gramm('x',C.Plate1.day,'y',C.Plate1.TKCPerY,...
    'color',C.Plate1.BactYeastPair,...
    'subset',C.Plate1.day~=0 & (C.Plate1.YeastFcn~="None" & C.Plate1.BactFcn~="None"));
g(1,2).facet_grid(C.Plate1.Perc_LW,[]);
g(1,2).stat_summary('geom','line');
g(1,2).geom_point();
g(1,2).axe_property('Yscale','log','Ylim',[1E-8 5E-2]);
g(1,2).set_title('trans Culture');
g(1,2).set_names('x','Day','y',["Transconjugant CFU/","Yeast Count"],...
    'color','Coculture Pair','row','LW%','column',[]);
g(1,2).set_order_options('row',[15 10 5 0],...
    'color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

C.Plate1 = Complete(1).Plate14.FlowData;
g(1,3) = gramm('x',C.Plate1.day,'y',C.Plate1.TKCPerY,...
    'color',C.Plate1.BactYeastPair,...
    'subset',C.Plate1.day~=0 & (C.Plate1.YeastFcn~="None" & C.Plate1.BactFcn~="None"));
g(1,3).facet_grid(C.Plate1.Perc_LW,[]);
g(1,3).stat_summary('geom','line');
g(1,3).geom_point();
g(1,3).axe_property('Yscale','log','Ylim',[1E-8 5E-2]);
g(1,3).set_title('trans Colonies');
g(1,3).set_names('x','Day','y',["Transconjugant CFU/","Yeast Count"],...
    'color','Coculture Pair','row','LW%','column',[]);
g(1,3).set_order_options('row',[100 10 5 0],...
    'color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

g.set_title('');

g.draw();


%% SF4 alt 1 for dissertation

clear g;
figure();

C.Plate1 = Complete(1).Plate3;
g(1,1) = gramm('x',C.Plate1.day,'y',C.Plate1.TKCPerY,...
    'color',C.Plate1.BactYeastPair,...
    'subset',C.Plate1.day~=0 & (C.Plate1.YeastFcn~="None" & C.Plate1.BactFcn~="None"));
g(1,1).facet_grid(C.Plate1.LW_Percent,[]);
g(1,1).stat_summary('geom','line');
g(1,1).geom_point();
g(1,1).axe_property('Yscale','log');
g(1,1).set_title('cis Culture');
g(1,1).set_names('x','Day','y','CFU/Yeast','color','Coculture Pair','row','LW%','column',[]);
g(1,1).set_order_options('row',[15 10 5 0],...
    'color',{'CrossB_CrossY','CrossB_wtY','wtB_CrossY','wtB_wtY'});

C.Plate1 = Complete(1).Plate7;
g(1,2) = gramm('x',C.Plate1.day,'y',C.Plate1.TKCPerY,...
    'color',C.Plate1.BactYeastPair,...
    'subset',C.Plate1.day~=0 & (C.Plate1.YeastFcn~="None" & C.Plate1.BactFcn~="None"));
g(1,2).facet_grid(C.Plate1.Perc_LW,[]);
g(1,2).stat_summary('geom','line');
g(1,2).geom_point();
g(1,2).axe_property('Yscale','log');
g(1,2).set_title('trans Culture');
g(1,2).set_names('x','Day','y','CFU/Yeast','color','Coculture Pair','row','LW%','column',[]);
g(1,2).set_order_options('row',[15 10 5 0],...
    'color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

g.set_title('Transconjugant Fraction of Yeast');

g.draw();

%% SF4 alt 2 for dissertation

clear g;
figure();

C.Plate1 = Complete(1).Plate14.FlowData;
g(1,1) = gramm('x',C.Plate1.day,'y',C.Plate1.TKCPerY,...
    'color',C.Plate1.BactYeastPair,...
    'subset',C.Plate1.day~=0 & (C.Plate1.YeastFcn~="None" & C.Plate1.BactFcn~="None"));
g(1,1).facet_grid(C.Plate1.Perc_LW,[]);
g(1,1).stat_summary('geom','line');
g(1,1).geom_point();
g(1,1).axe_property('Yscale','log');
g(1,1).set_title('Transconjugant Fraction of Yeast, trans Colonies');
g(1,1).set_names('x','Day','y','CFU/Yeast','color','Coculture Pair','row','LW%','column',[]);
g(1,1).set_order_options('row',[100 10 5 0],...
    'color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

g.set_title('Transconjugant Fraction of Yeast, trans Colonies');

g.draw();



%% SF5: Clumps ellipse fits, area v. coincident bact

dataS = Complete(1).Plate16.ImageData;

clear g;
figure();

g(1,1) = gramm('x',dataS.Areas,'y',dataS.Robjs,'color',dataS.BactYeastPair,'subset',...
    (dataS.YeastFcn~='None' & dataS.BactFcn~='None' & dataS.Dilution~=1));
g(1,1).facet_grid(dataS.PercLW,dataS.Mannose);
g(1,1).stat_ellipse();
g(1,1).axe_property('YLim',[0 4],'XLim',[0 600]);

% g(1,1).set_title('Yeast event sizes');
g(1,1).set_names('x','Yeast Event Area (px^2)','y',["Coincident Bact","Per Yeast Event"],...
    'color','Coculture Pair','row','LW%','column','Mannose');
g(1,1).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'},'row',[15 10 5 0],'column',{'N','Y'});
g(1,1).set_point_options('base_size',2);

% g.set_title(["Yeast event area by coincident bact count","(95% of events within shaded regions)"]);

g.draw();



%% SF6: Clumping_Histos

dataS = Complete(1).Plate16.ImageData;

%Index with cols as pairings, rows as LW%. First Mannose-, then +
LWPerc = [15 10 5 0];
Pair = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'};


    for j = 1:4
        p = Pair{j};

        IdxM{1,j} = find(dataS.Day==6 & dataS.BactYeastPair==p & dataS.Mannose=='N' & dataS.Dilution~=1);
        IdxP{1,j} = find(dataS.Day==6 & dataS.BactYeastPair==p & dataS.Mannose=='Y' & dataS.Dilution~=1);

        pos = j;
        
        figure(1)
        h = subplot(1,4,pos);
        hold on
        
        hM = histogram(cell2mat(cellfun(@(x)x(:),dataS.Robjs(IdxM{1,j}),'un',0)),...
            'FaceAlpha',0.4,'BinWidth',1);
        hP = histogram(cell2mat(cellfun(@(x)x(:),dataS.Robjs(IdxP{1,j}),'un',0)),...
            'FaceAlpha',0.4,'BinWidth',1);
        xlabel('Number of coincident bacteria');
        
        h.XLim = [0 10];
        h.YScale = 'log';
        hold off
        
        title(p,'Interpreter','none');

%         sgtitle('Bacteria co-existing with yeast events, day6 (all LW%)');
        legend('Mannose-','Mannose+');

        
        Mmed = mean(cell2mat(cellfun(@(x)x(:),dataS.Areas(IdxM{1,j}),'UniformOutput',false)));
        Mlab = {['Mean=' num2str(Mmed)]};
        Pmed = mean(cell2mat(cellfun(@(x)x(:),dataS.Areas(IdxP{1,j}),'UniformOutput',false)));
        Plab = {['Mean=' num2str(Pmed)]};
        
        figure(2)
        a = subplot(1,4,pos);
        hold on
        
        aM = histogram(cell2mat(cellfun(@(x)x(:),dataS.Areas(IdxM{1,j}),'un',0)),...
            'FaceAlpha',0.4,'BinWidth',10);
        aP = histogram(cell2mat(cellfun(@(x)x(:),dataS.Areas(IdxP{1,j}),'un',0)),...
            'FaceAlpha',0.4,'BinWidth',10);
        xlabel('Yeast event size (px^2)');
        
        M = xline(Mmed, 'Color',[0 0 1],'label',Mlab);
        P = xline(Pmed, 'Color',[1 .5 0],'label',Plab);

        if Mmed>Pmed
            P.LabelHorizontalAlignment = 'left';
        else
            M.LabelHorizontalAlignment = 'left';
        end

        a.XLim = [0 500];
        hold off
        
        title(p,'Interpreter','none');

%         sgtitle('Yeast event areas, day6 (all LW%)');
        legend('Mannose-','Mannose+');

    end


%% SF7: Clumps D:R, TKC

clear g;
figure();

EndPoints.Plate23 = Complete(1).Plate16.TecanData;

maxCh = max(EndPoints.Plate23.Cherry(EndPoints.Plate23.day>0));
maxCit = max(EndPoints.Plate23.Citrine(EndPoints.Plate23.day>0));

EndPoints.Plate23.ChNorm = EndPoints.Plate23.Cherry ./ maxCh;
EndPoints.Plate23.CitNorm = EndPoints.Plate23.Citrine ./ maxCit;
EndPoints.Plate23.BYnormRatio = EndPoints.Plate23.ChNorm ./ EndPoints.Plate23.CitNorm;

EndPoints.Plate23.TKCPerRatio = EndPoints.Plate23.TKC ./ EndPoints.Plate23.BYnormRatio;
EndPoints.Plate23.TKCPerCitrine = EndPoints.Plate23.TKC ./ EndPoints.Plate23.Citrine;

g(1,1) = gramm('x',EndPoints.Plate23.day,'y',EndPoints.Plate23.BYnormRatio,...
    'color',EndPoints.Plate23.BactYeastPair,...
    'subset',(EndPoints.Plate23.YeastFcn~="None" & EndPoints.Plate23.BactFcn~="None"));
g(1,1).facet_grid(EndPoints.Plate23.PercLW,EndPoints.Plate23.Mannose);
g(1,1).stat_summary('type','quartile');
g(1,1).geom_point();
g(1,1).axe_property('Yscale','log');
g(1,1).set_title('D:R ratios');
g(1,1).set_names('x','Day','y',["Normalized bact FU (mCherry) /","Normalized yeast FU (ymCitrine)"],...
    'color','Coculture Pair','row','LW%','column','Mannose');
g(1,1).set_order_options('row',[15 10 5 0],'column',{'N','Y'},...
    'color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

g(1,2) = gramm('x',EndPoints.Plate23.day,'y',EndPoints.Plate23.TKC,...
    'color',EndPoints.Plate23.BactYeastPair,...
    'subset',(EndPoints.Plate23.YeastFcn~="None" & EndPoints.Plate23.BactFcn~="None"));
g(1,2).facet_grid(EndPoints.Plate23.PercLW,EndPoints.Plate23.Mannose);
g(1,2).stat_summary('geom','line');
g(1,2).geom_point('dodge',1);
% g(1,2).axe_property('Yscale','log');
g(1,2).set_title('Raw IDC Counts');
g(1,2).set_names('x','Day','y','Transconjugant CFU','color','Coculture Pair','row','LW%','column','Mannose');
g(1,2).set_order_options('row',[15 10 5 0],'column',{'N','Y'},...
    'color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

% g.set_title('Mannose Plates Cell Ratios and TKC');

g.draw();

%% SF7: defense alt

clear g;
figure();

EndPoints.Plate23 = Complete(1).Plate16.TecanData;

maxCh = max(EndPoints.Plate23.Cherry(EndPoints.Plate23.day>0));
maxCit = max(EndPoints.Plate23.Citrine(EndPoints.Plate23.day>0));

EndPoints.Plate23.ChNorm = EndPoints.Plate23.Cherry ./ maxCh;
EndPoints.Plate23.CitNorm = EndPoints.Plate23.Citrine ./ maxCit;
EndPoints.Plate23.BYnormRatio = EndPoints.Plate23.ChNorm ./ EndPoints.Plate23.CitNorm;

EndPoints.Plate23.TKCPerRatio = EndPoints.Plate23.TKC ./ EndPoints.Plate23.BYnormRatio;
EndPoints.Plate23.TKCPerCitrine = EndPoints.Plate23.TKC ./ EndPoints.Plate23.Citrine;

g(1,1) = gramm('x',EndPoints.Plate23.day,'y',EndPoints.Plate23.BYnormRatio,...
    'color',EndPoints.Plate23.BactYeastPair,...
    'subset',(EndPoints.Plate23.YeastFcn~="None" & EndPoints.Plate23.BactFcn~="None"));
g(1,1).facet_grid(EndPoints.Plate23.Mannose,EndPoints.Plate23.PercLW);
g(1,1).stat_summary('type','quartile');
g(1,1).geom_point();
g(1,1).axe_property('Yscale','log');
g(1,1).set_title('D:R ratio (normalized Cherry/Citrine)');
g(1,1).set_names('x','Day','y',["Normalized mCherry /","Normalized ymCitrine"],'color','Coculture Pair',...
    'column','LW%','row','Mannose');
g(1,1).set_order_options('column',[15 10 5 0],'row',{'N','Y'},...
    'color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'});

g.set_title('Donor-to-recipient ratios');

g.draw();


%% SF12: Colony Radial_Day6

data = Complete(1).Plate14.ImageData;

clear g;
figure();

g(1,1) = gramm('x',data.RadiusMillimetersV,'y',data.RadialAvgR,'color',data.BactYeastPair,'subset',data.Day==6);
g(1,1).facet_grid(data.Perc_LW,[]);
g(1,1).stat_summary('type','std');
g(1,1).axe_property('XLim',[0 3]);
g(1,1).set_title('Average mCherry per Radius');
g(1,1).set_names('x','Radius (mm)','y','Average Red Intensity','color','Coculture Pair','row','LW%');
g(1,1).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'},...
    'row',[100 10 5 0]);
g(1,1).set_stat_options('alpha',0.001);
g(1,1).set_line_options('base_size',2);

g(1,2) = gramm('x',data.RadiusMillimetersV,'y',data.RadialAvgG,'color',data.BactYeastPair,'subset',data.Day==6);
g(1,2).facet_grid(data.Perc_LW,[]);
g(1,2).stat_summary('type','std');
g(1,2).axe_property('XLim',[0 3]);
g(1,2).set_title('Average ymCitrine per Radius');
g(1,2).set_names('x','Radius (mm)','y','Average Yellow Intensity','color','Coculture Pair','row','LW%');
g(1,2).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'},...
    'row',[100 10 5 0]);
g(1,2).set_stat_options('alpha',0.001);
g(1,2).set_line_options('base_size',2);

g(1,3) = gramm('x',data.RadiusMillimetersV,'y',data.RadialICQ,'color',data.BactYeastPair,'subset',data.Day==6);
g(1,3).facet_grid(data.Perc_LW,[]);
g(1,3).stat_summary('type','std');
g(1,3).axe_property('XLim',[0 3]);
g(1,3).set_title('Colocalization per Radius');
g(1,3).set_names('x','Radius (mm)','y','ICQ','color','Coculture Pair','row','LW%');
g(1,3).set_order_options('color',{'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'},...
    'row',[100 10 5 0]);
g(1,3).set_stat_options('alpha',0.001);
g(1,3).set_line_options('base_size',2);

% g.set_title('Radial Signals, Day 6');

g.draw();



%% SF14: Rescue batch culture growths

clear g;
figure();

g(1,1) = gramm('x',Complete(1).All.day,'y',Complete(1).All.Bnormalized,...
    'linestyle',Complete(1).All.YeastFcn,'size',Complete(1).All.YeastFcn,...
    'subset',(Complete(1).All.BactFcn=="crossB" & ...
    Complete(1).All.YeastFcn~="wtY" & Complete(1).All.Assay=="Rescue" &...
    Complete(1).All.day~=0));
g(1,1).facet_grid(Complete(1).All.Perc_AA,[]);
g(1,1).stat_summary('type','std');
g(1,1).set_title('crossB_crossB');
g(1,1).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair','row','UH%','column',[]);
g(1,1).set_order_options('row',[15 10 5 0],...
    'linestyle',{'crossY','None'},'size',{'None','crossY'});
g(1,1).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,1).set_color_options('map',[1 0 0],'legend','expand');
g(1,1).no_legend();

g(1,1).update('y',Complete(1).All.Ynormalized,'linestyle',Complete(1).All.BactFcn,...
    'size',Complete(1).All.BactFcn,'subset',(Complete(1).All.YeastFcn=="crossY" & ...
    Complete(1).All.BactFcn~="wtB" & Complete(1).All.Assay=="Rescue" &...
    Complete(1).All.day~=0));
g(1,1).stat_summary('type','std');
g(1,1).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair');
g(1,1).set_order_options('linestyle',{'crossB','None'},'size',{'None','crossB'},...
    'row',[15 10 5 0])
g(1,1).set_color_options('map',[.65 .5 0],'legend','expand');
g(1,1).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,1).no_legend();

g(1,2) = gramm('x',Complete(1).All.day,'y',Complete(1).All.Bnormalized,...
    'linestyle',Complete(1).All.YeastFcn,'size',Complete(1).All.YeastFcn,...
    'subset',(Complete(1).All.BactFcn=="crossB" &...
    Complete(1).All.YeastFcn~="crossY" & Complete(1).All.Assay=="Rescue" &...
    Complete(1).All.day~=0));
g(1,2).facet_grid(Complete(1).All.Perc_AA,[]);
g(1,2).stat_summary('type','std');
g(1,2).set_title('crossB_wtY');
g(1,2).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair','row','UH%','column',[]);
g(1,2).set_order_options('row',[15 10 5 0],'size',{'None','wtY'},...
    'linestyle',{'wtY','None'});
g(1,2).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,2).set_color_options('map',[1 0 0],'legend','expand');
g(1,2).no_legend();

g(1,2).update('y',Complete(1).All.Ynormalized,'linestyle',Complete(1).All.BactFcn,...
    'size',Complete(1).All.BactFcn,'subset',(Complete(1).All.YeastFcn=="wtY" & ...
    Complete(1).All.BactFcn~="wtB" & Complete(1).All.Assay=="Rescue" &...
    Complete(1).All.day~=0));
g(1,2).stat_summary('type','std');
g(1,2).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair');
g(1,2).set_order_options('row',[15 10 5 0],'linestyle',{'crossB','None'},'size',{'None','crossB'});
g(1,2).set_color_options('map',[.65 .5 0],'legend','expand');
g(1,2).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,2).no_legend();

g(1,3) = gramm('x',Complete(1).All.day,'y',Complete(1).All.Bnormalized,...
    'linestyle',Complete(1).All.YeastFcn,'size',Complete(1).All.YeastFcn,...
    'subset',(Complete(1).All.BactFcn=="wtB") & ...
    Complete(1).All.YeastFcn~="wtY" & Complete(1).All.Assay=="Rescue" &...
    Complete(1).All.day~=0);
g(1,3).facet_grid(Complete(1).All.Perc_AA,[]);
g(1,3).stat_summary('type','std');
g(1,3).set_title('wtB_crossY');
g(1,3).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair','row','UH%','column',[]);
g(1,3).set_order_options('row',[15 10 5 0],'size',{'None','crossY'},...
    'linestyle',{'crossY','None'});
g(1,3).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,3).set_color_options('map',[1 0 0],'legend','expand');
g(1,3).no_legend();

g(1,3).update('y',Complete(1).All.Ynormalized,'linestyle',Complete(1).All.BactFcn,...
    'size',Complete(1).All.BactFcn,'subset',(Complete(1).All.YeastFcn=="crossY" &...
    Complete(1).All.BactFcn~="crossB" & Complete(1).All.Assay=="Rescue" &...
    Complete(1).All.day~=0));
g(1,3).stat_summary('type','std');
g(1,3).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair');
g(1,3).set_order_options('linestyle',{'wtB','None'},'size',{'None','wtB'},'row',[15 10 5 0]);
g(1,3).set_color_options('map',[.65 .5 0],'legend','expand');
g(1,3).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,3).no_legend();

g(1,4) = gramm('x',Complete(1).All.day,'y',Complete(1).All.Bnormalized,'linestyle',Complete(1).All.YeastFcn,...
    'size',Complete(1).All.YeastFcn,'subset',(Complete(1).All.BactFcn=="wtB") & ...
    Complete(1).All.YeastFcn~="crossY" & Complete(1).All.Assay=="Rescue" &...
    Complete(1).All.day~=0);
g(1,4).facet_grid(Complete(1).All.Perc_AA,[]);
g(1,4).stat_summary('type','std');
g(1,4).set_title('wtB_wtY');
g(1,4).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair','row','UH%','column',[]);
g(1,4).set_order_options('row',[15 10 5 0],'size',{'None','wtY'},...
    'linestyle',{'wtY','None'});
g(1,4).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,4).set_color_options('map',[1 0 0],'legend','expand');
g(1,4).no_legend();

g(1,4).update('y',Complete(1).All.Ynormalized,'linestyle',Complete(1).All.BactFcn,...
    'size',Complete(1).All.BactFcn,'subset',(Complete(1).All.YeastFcn=="wtY" & ...
    Complete(1).All.BactFcn~="crossB" & Complete(1).All.Assay=="Rescue" &...
    Complete(1).All.day~=0));
g(1,4).stat_summary('type','std');
g(1,4).set_names('x','Day','y','Normalized cells','linestyle','Coculture Pair');
g(1,4).set_order_options('row',[15 10 5 0],'linestyle',{'wtB','None'},'size',{'None','wtB'});
g(1,4).set_color_options('map',[.65 .5 0],'legend','expand');
g(1,4).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,4).no_legend();


% g.set_title('Growth: Rescue All Pairs');

g.draw();


%% SF18: CRISPR Growth

AllData.Plate1 = Complete(1).Plate15.TecanData;

maxCh = max(AllData.Plate1.Cherry(AllData.Plate1.day>1));
maxCit = max(AllData.Plate1.Citrine(AllData.Plate1.day>1));

clear g;
figure();
g(1,1) = gramm('x',AllData.Plate1.time./3600,'y',AllData.Plate1.Cherry./maxCh,...
    'linestyle',AllData.Plate1.YeastFcn,'size',AllData.Plate1.YeastFcn,...
    'subset',AllData.Plate1.BactFcn=='cutB' & AllData.Plate1.PercU==0);
g(1,1).facet_grid(AllData.Plate1.PercW,[]);
g(1,1).stat_summary('type','std');
g(1,1).set_title('Cutter donor');
g(1,1).axe_property('YLim',[0 1]);
g(1,1).set_names('x','t (hr)','y','Normalized FU','linestyle','Coculture Pair','row','W%');
g(1,1).set_order_options('row',[15 10 5 1],...
    'linestyle',{'crisprY','None'},'size',{'None','crisprY'});
g(1,1).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,1).set_color_options('map',[1 0 0],'legend','expand');
g(1,1).no_legend();

g(1,1).update('y',AllData.Plate1.Citrine./maxCit,'linestyle',AllData.Plate1.BactFcn,...
    'size',AllData.Plate1.BactFcn,'subset',AllData.Plate1.YeastFcn=='crisprY' & ...
    AllData.Plate1.BactFcn~='nocutB' & AllData.Plate1.PercU==0);
g(1,1).stat_summary('type','std');
g(1,1).geom_vline('xintercept',509335/3600);
g(1,1).set_names('x','t (hr)','y','Normalized FU','linestyle','Coculture Pair');
g(1,1).set_order_options('linestyle',{'cutB','None'},'size',{'None','cutB'},...
    'row',[15 10 5 1]);
g(1,1).set_color_options('map',[.7 .6 0],'legend','expand');
g(1,1).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,1).no_legend();

g(1,2) = gramm('x',AllData.Plate1.time./3600,'y',AllData.Plate1.Cherry./maxCh,...
    'linestyle',AllData.Plate1.YeastFcn,'size',AllData.Plate1.YeastFcn,...
    'subset',AllData.Plate1.BactFcn=='nocutB' & AllData.Plate1.PercU==0);
g(1,2).facet_grid(AllData.Plate1.PercW,[]);
g(1,2).stat_summary('type','std');
g(1,2).set_title('No-oriT donor');
g(1,2).axe_property('YLim',[0 1]);
g(1,2).set_names('x','t (hr)','y','Normalized FU','linestyle','Coculture Pair','row','W%');
g(1,2).set_order_options('row',[15 10 5 1],...
    'linestyle',{'crisprY','None'},'size',{'None','crisprY'});
g(1,2).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,2).set_color_options('map',[1 0 0],'legend','expand');
g(1,2).no_legend();

g(1,2).update('y',AllData.Plate1.Citrine./maxCit,'linestyle',AllData.Plate1.BactFcn,...
    'size',AllData.Plate1.BactFcn,'subset',AllData.Plate1.YeastFcn=='crisprY' & ...
    AllData.Plate1.BactFcn~='cutB' & AllData.Plate1.PercU==0);
g(1,2).stat_summary('type','std');
g(1,2).geom_vline('xintercept',509335/3600);
g(1,2).set_names('x','t (hr)','y','Normalized FU','linestyle','Coculture Pair');
g(1,2).set_order_options('linestyle',{'nocutB','None'},'size',{'None','nocutB'},...
    'row',[15 10 5 1]);
g(1,2).set_color_options('map',[.7 .6 0],'legend','expand');
g(1,2).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
g(1,2).no_legend();


% g.set_title('CRISPR Plate Growth (normalized)');

g.draw();


%% SI: CRISPR BFP 

AllData.Plate1 = Complete(1).Plate15.TecanEndpoints;

clear g;
figure();
g(1,1) = gramm('x',AllData.Plate1.time./3600,'y',AllData.Plate1.NormBFP,...
   'linestyle',AllData.Plate1.BactFcn,...
   'subset',AllData.Plate1.YeastFcn=='crisprY' & AllData.Plate1.PercU==100);
g(1,1).facet_grid(AllData.Plate1.PercW,[]);
g(1,1).stat_summary('type','std');
g(1,1).geom_vline('xintercept',564036/3600);
g(1,1).set_title('Cutter donor');
g(1,1).axe_property('YLim',[0 1]);
g(1,1).set_names('x','t (hr)','y','Normalized FU','row','W%','linestyle','Yeast Pairing');
g(1,1).set_order_options('row',[15 10 5 1],'linestyle',{'cutB','nocutB','None'});
g(1,1).set_line_options('base_size',2,'step_size',2,'styles',{'-','-.',':'});
g(1,1).set_color_options('map',[0 0 .7],'legend','expand');

g.set_title('CRISPR Yeast Growth (normalized)');

g.draw();

%% SI19: CRISPR Fold Plasmid loss Yeast


EndPoints.Plate1 = Complete(1).Plate15.TecanEndpoints;

clear g;
figure();

g(1,1) = gramm('x',EndPoints.Plate1.day,'y',EndPoints.Plate1.FoldBFPDecreaseMono,...
    'color',EndPoints.Plate1.BactYeastPair,...
    'subset',(EndPoints.Plate1.YeastFcn~="None" & EndPoints.Plate1.PercU==100));
g(1,1).facet_grid(EndPoints.Plate1.PercW,[]);
g(1,1).stat_summary('geom','bar');
g(1,1).geom_vline('xintercept',6.5);
g(1,1).axe_property('Ylim',[-10 -.9]);
% g(1,1).set_title('Fold decrease in BFP-URA from Monoculture, 100% Ura');
g(1,1).set_names('x','Day','y',["-Monoculture Control BFP","/ BFP"],'color','Coculture Pair','row','W%');
g(1,1).set_order_options('row',[15 10 5 1],...
    'color',{'cutB_crisprY','nocutB_crisprY'});


% g.set_title('Comparing Fold-Loss of BFP-URA from Yeast Cocultures (no URA selection)');

g.draw();


%% Dissertation 3.4.2.I: TKC w/ diff agar %


clear g;
figure();

g(1,1) = gramm('x',Complete(1).Plate6.day,'y',Complete(1).Plate6.TKC,...
    'color',Complete(1).Plate6.BactYeastPair,...
    'subset',(Complete(1).Plate6.YeastFcn~="None" & Complete(1).Plate6.BactFcn~="None" &...
    Complete(1).Plate6.Agar_Percent==0.2));
g(1,1).facet_grid(Complete(1).Plate6.LW_Percent,[]);
g(1,1).stat_summary('geom','line');
g(1,1).geom_point('dodge',1);
g(1,1).axe_property('Yscale','log','YLim',[.1 7000]);
g(1,1).set_title('TKC in motile colonies, 0.2% Agar');
g(1,1).set_names('x','Day','y','CFU','color','Coculture Pair','row','LW%');
g(1,1).set_order_options('row',[10 0],...
    'color',{'CrossB_CrossY','CrossB_wtY','wtB_CrossY','wtB_wtY'});

g(1,2) = gramm('x',Complete(1).Plate6.day,'y',Complete(1).Plate6.TKC,...
    'color',Complete(1).Plate6.BactYeastPair,...
    'subset',(Complete(1).Plate6.YeastFcn~="None" & Complete(1).Plate6.BactFcn~="None" &...
    Complete(1).Plate6.Agar_Percent==2));
g(1,2).facet_grid(Complete(1).Plate6.LW_Percent,[]);
g(1,2).stat_summary('geom','line');
g(1,2).geom_point('dodge',1);
g(1,2).axe_property('Yscale','log','YLim',[.1 7000]);
g(1,2).set_title('TKC in non-motile colonies, 2% Agar');
g(1,2).set_names('x','Day','y','CFU','color','Coculture Pair','row','LW%');
g(1,2).set_order_options('row',[10 0],...
    'color',{'CrossB_CrossY','CrossB_wtY','wtB_CrossY','wtB_wtY'});

g.set_title('Mixed colony TKC at different agar %');

g.draw();


%% END OF FIGURES








% Next up... data consolidation, statistical tests








%% Consolidate culture experiments

Vars = Complete(1).AllCulture.Properties.VariableNames;

% Complete(1).AllCulture = table();

for plate = [3 4 5 7 8 12 13 15]
    plt = ['Plate' num2str(plate)];
    
    temp = table();
    Complete(1).(plt).Plate(:) = plate;
    
    for V = 1:numel(Vars)
        var = Vars{V};
        
        temp.(var) = Complete(1).(plt).(var);
        
    end
    
    Complete(1).AllCulture = vertcat(Complete(1).AllCulture,temp);

end

%% Consolidate colony experiments

Vars = Complete(1).AllColonies.Properties.VariableNames;
% Complete(1).AllCulture = table();

for plate = [6 11 14]
    plt = ['Plate' num2str(plate)];
    
    temp = table();
    Complete(1).(plt).Plate(:) = plate;
    
    for V = 1:numel(Vars)
        var = Vars{V};
        
        temp.(var) = Complete(1).(plt).(var);
        
    end
    
    Complete(1).AllColonies = vertcat(Complete(1).AllColonies,temp);

end

%% Consolidate Rescue Plates

Cols = {'experiment','Plate','day','well','Yeast','Bacteria','YeastFcn','BactFcn','BactYeastPair','Perc_UH',...
    'TKC','Bnormalized','Ynormalized','BYnormRatio','TKCPerRatioNorm','Contaminated'};

Rescues = [];

for p = [4 8]
    plate = ['Plate' num2str(p)];
    
    Rescues0 = table();
    
    for c = 1:numel(Cols)
        col = Cols{c};

        Rescues0.(col)(:) = Complete(1).(plate).(col)(:);
        
    end
    
    Rescues = vertcat(Rescues,Rescues0);
    
end

Complete(1).AllRescues = Rescues;
Complete(2).AllRescues = {'Flow Plates 4,8'};

%% Calculate log versions of data

Complete(1).AllCulture.logBYnormRatio = log(Complete(1).AllCulture.BYnormRatio);
Complete(1).AllCulture.logTKC = log(Complete(1).AllCulture.TKC);
Complete(1).AllCulture.logTKCPerRatioNorm = log(Complete(1).AllCulture.TKCPerRatioNorm);

Complete(1).AllCulture.logBYnormRatio(Complete(1).AllCulture.logBYnormRatio==Inf) = nan;
Complete(1).AllCulture.logTKC(Complete(1).AllCulture.logTKC==Inf) = nan;
Complete(1).AllCulture.logTKCPerRatioNorm(Complete(1).AllCulture.logTKCPerRatioNorm==Inf) = nan;
Complete(1).AllCulture.logBYnormRatio(Complete(1).AllCulture.logBYnormRatio==-Inf) = nan;
Complete(1).AllCulture.logTKC(Complete(1).AllCulture.logTKC==-Inf) = nan;
Complete(1).AllCulture.logTKCPerRatioNorm(Complete(1).AllCulture.logTKCPerRatioNorm==-Inf) = nan;


Complete(1).AllColonies.logBYnormRatio = log(Complete(1).AllColonies.BYnormRatio);
Complete(1).AllColonies.logTKC = log(Complete(1).AllColonies.TKC);
Complete(1).AllColonies.logTKCPerRatioNorm = log(Complete(1).AllColonies.TKCPerRatioNorm);

Complete(1).AllColonies.logBYnormRatio(Complete(1).AllColonies.logBYnormRatio==Inf) = nan;
Complete(1).AllColonies.logTKC(Complete(1).AllColonies.logTKC==Inf) = nan;
Complete(1).AllColonies.logTKCPerRatioNorm(Complete(1).AllColonies.logTKCPerRatioNorm==Inf) = nan;
Complete(1).AllColonies.logBYnormRatio(Complete(1).AllColonies.logBYnormRatio==-Inf) = nan;
Complete(1).AllColonies.logTKC(Complete(1).AllColonies.logTKC==-Inf) = nan;
Complete(1).AllColonies.logTKCPerRatioNorm(Complete(1).AllColonies.logTKCPerRatioNorm==-Inf) = nan;


%% Designate assay types

for i = [3 5 7 12 13]
    Complete(1).AllCulture.Assay(Complete(1).AllCulture.Plate == i) = categorical("Dynamics");
end

for j = [4 8]
    Complete(1).AllCulture.Assay(Complete(1).AllCulture.Plate == j) = categorical("Rescue");
end

%% Combine all used

Complete(1).All = vertcat(Complete(1).AllCulture,Complete(1).AllColonies);

%% Get B, Y normalized

Plates = unique(Complete(1).All.Plate);
Pcol = {'LW_Percent','Perc_UH','LW_Percent','LW_Percent','Perc_LW','Perc_UH','','Perc_LW','Perc_LW','Perc_LW',''};


for p = 1:numel(Plates)
    plate = Plates(p);
    plt = ['Plate' num2str(plate)];    
    
    Wells = unique(Complete(1).All.well(Complete(1).All.Plate==plate));
    pcol = Pcol{p};
    
    for w = 1:numel(Wells)
        well = Wells(w);
        
        Days = unique(Complete(1).All.day(Complete(1).All.Plate==plate & Complete(1).All.well==well));
        
        for d = 1:numel(Days)
            day = Days(d);
            
            if day>0 && plate~=11 && plate~=15 && plate~=14 && plate~=6
                Bread = Complete(1).(plt).Bnormalized(Complete(1).(plt).well==well &...
                    Complete(1).(plt).day==day & ~isnan(Complete(1).(plt).Bnormalized));
                Yread = Complete(1).(plt).Ynormalized(Complete(1).(plt).well==well &...
                    Complete(1).(plt).day==day & ~isnan(Complete(1).(plt).Ynormalized));
                Y100 = Complete(1).(plt).YPer100Mic(Complete(1).(plt).well==well &...
                    Complete(1).(plt).day==day & ~isnan(Complete(1).(plt).YPer100Mic));
                B100 = Complete(1).(plt).BPer100Mic(Complete(1).(plt).well==well &...
                    Complete(1).(plt).day==day & ~isnan(Complete(1).(plt).BPer100Mic));
                
                
                Perc = Complete(1).(plt).(pcol)(Complete(1).(plt).well==well &...
                    Complete(1).(plt).day==day);
                
                if isempty(Bread)
                Complete(1).All.Bnormalized(Complete(1).All.Plate==plate &...
                    Complete(1).All.well==well & Complete(1).All.day==day) = nan;
                else
                Complete(1).All.Bnormalized(Complete(1).All.Plate==plate &...
                    Complete(1).All.well==well & Complete(1).All.day==day) = Bread;
                end

                if isempty(Yread)
                Complete(1).All.Ynormalized(Complete(1).All.Plate==plate &...
                    Complete(1).All.well==well & Complete(1).All.day==day) = nan;
                else
                Complete(1).All.Ynormalized(Complete(1).All.Plate==plate &...
                    Complete(1).All.well==well & Complete(1).All.day==day) = Yread;
                end
                
                if isempty(Y100)
                Complete(1).All.YPer100Mic(Complete(1).All.Plate==plate &...
                    Complete(1).All.well==well & Complete(1).All.day==day) = nan;
                else
                Complete(1).All.YPer100Mic(Complete(1).All.Plate==plate &...
                    Complete(1).All.well==well & Complete(1).All.day==day) = Y100;
                end
                
                if isempty(B100)
                Complete(1).All.BPer100Mic(Complete(1).All.Plate==plate &...
                    Complete(1).All.well==well & Complete(1).All.day==day) = nan;
                else
                Complete(1).All.BPer100Mic(Complete(1).All.Plate==plate &...
                    Complete(1).All.well==well & Complete(1).All.day==day) = B100;
                end

                
                Complete(1).All.Perc_AA(Complete(1).All.Plate==plate &...
                    Complete(1).All.well==well & Complete(1).All.day==day) = Perc;
            end
            
        end
    end
end

%% Fill in gaps


Plates = 11;
FlowPlate = {'A','B'};

for p = 1:numel(Plates)
    plate = Plates(p);
    plt = ['Plate' num2str(plate)];    
        
    Days = unique(Complete(1).All.day(Complete(1).All.Plate==plate));
                
    for d = 1:numel(Days)
        day = Days(d);
            
        Wells = unique(Complete(1).All.well(Complete(1).All.Plate==plate & Complete(1).All.day==day));

        for w = 1:numel(Wells)
            well = Wells(w);
        
            BYs = unique(Complete(1).All.BactYeastPair(Complete(1).All.Plate==plate &...
                Complete(1).All.well==well & Complete(1).All.day==day));

            for pr = 1:numel(BYs)
                by = BYs(pr);
                
                if day>0 
                    Bread = Complete(1).(plt).FlowData.Bnormalized(Complete(1).(plt).FlowData.well==well &...
                        Complete(1).(plt).FlowData.day==day & Complete(1).(plt).FlowData.BactYeastPair==by &...
                        ~isnan(Complete(1).(plt).FlowData.Bnormalized));
                    Yread = Complete(1).(plt).FlowData.Ynormalized(Complete(1).(plt).FlowData.well==well &...
                        Complete(1).(plt).FlowData.day==day & Complete(1).(plt).FlowData.BactYeastPair==by &...
                        ~isnan(Complete(1).(plt).FlowData.Ynormalized));
                    Y100 = Complete(1).(plt).FlowData.YPer100Mic(Complete(1).(plt).FlowData.well==well &...
                        Complete(1).(plt).FlowData.day==day & Complete(1).(plt).FlowData.BactYeastPair==by &...
                        ~isnan(Complete(1).(plt).FlowData.YPer100Mic));
                    
                    
                    
                    Perc = Complete(1).(plt).FlowData.Perc_LW(Complete(1).(plt).FlowData.well==well &...
                        Complete(1).(plt).FlowData.day==day & Complete(1).(plt).FlowData.BactYeastPair==by);
                    Type = Complete(1).(plt).FlowData.Type(Complete(1).(plt).FlowData.well==well &...
                        Complete(1).(plt).FlowData.day==day & Complete(1).(plt).FlowData.BactYeastPair==by);
                    
                    
                    if isempty(Bread)
                    Complete(1).All.Bnormalized(Complete(1).All.Plate==plate &...
                        Complete(1).All.well==well & Complete(1).All.BactYeastPair==by &...
                        Complete(1).All.day==day) = nan;
                    else
                    Complete(1).All.Bnormalized(Complete(1).All.Plate==plate &...
                        Complete(1).All.well==well & Complete(1).All.BactYeastPair==by &...
                        Complete(1).All.day==day) = Bread;
                    end

                    if isempty(Yread)
                    Complete(1).All.Ynormalized(Complete(1).All.Plate==plate &...
                        Complete(1).All.well==well & Complete(1).All.BactYeastPair==by &...
                        Complete(1).All.day==day) = nan;
                    else
                    Complete(1).All.Ynormalized(Complete(1).All.Plate==plate &...
                        Complete(1).All.well==well & Complete(1).All.BactYeastPair==by &...
                        Complete(1).All.day==day) = Yread;
                    end
                    
                    if isempty(Y100)
                    Complete(1).All.YPer100Mic(Complete(1).All.Plate==plate &...
                        Complete(1).All.well==well & Complete(1).All.BactYeastPair==by &...
                        Complete(1).All.day==day) = nan;
                    else
                    Complete(1).All.YPer100Mic(Complete(1).All.Plate==plate &...
                        Complete(1).All.well==well & Complete(1).All.BactYeastPair==by &...
                        Complete(1).All.day==day) = Y100;
                    end

                    
                    Complete(1).All.Perc_AA(Complete(1).All.Plate==plate &...
                        Complete(1).All.well==well & Complete(1).All.BactYeastPair==by &...
                        Complete(1).All.day==day) = Perc;
                    Complete(1).All.Txfer(Complete(1).All.Plate==plate &...
                        Complete(1).All.well==well & Complete(1).All.BactYeastPair==by &...
                        Complete(1).All.day==day) = Type;

                end
            end
        end
    end
end

Complete(1).All.LimAA(Complete(1).All.Plate==11) = categorical("LW");


%% Enter cis trans

for plate = [3:6 8]
    Complete(1).All.Txfer(Complete(1).All.Plate==plate) = categorical("Cis");
end

for plate = [7 12:15]
    Complete(1).All.Txfer(Complete(1).All.Plate==plate) = categorical("Trans");
end
   
%% Error and significance, mannose plate


% CI calculated via stat_summary in F3D plot code, copied here for
% convenience
PlusCI = [-0.0435574698629874,0.0866543113897658;0.0349021713123923,0.140782562250965;-0.0970576488157222,0.292266344848396;0.000587153859681278,0.377396602224529];
MinCI = [1.74707734078518,5.24365045960643;1.67056234899049,10.0743241597373;0.585023557887802,5.95156450721151;-1.18144077707485,26.5655271606254];
LW = [0 5 10 15];

for i = 1:4
    lw = LW(i);
    
    IdxMin = find(Complete(1).Plate16.TecanData.PercLW==lw &...
        Complete(1).Plate16.TecanData.BactYeastPair=='crossB_wtY' &...
        Complete(1).Plate16.TecanData.day==6 &...
        Complete(1).Plate16.TecanData.Mannose=='N');
    
    IdxPlus = find(Complete(1).Plate16.TecanData.PercLW==lw &...
        Complete(1).Plate16.TecanData.BactYeastPair=='crossB_wtY' &...
        Complete(1).Plate16.TecanData.day==6 &...
        Complete(1).Plate16.TecanData.Mannose=='Y');
    
    IdxAll = vertcat(IdxMin,IdxPlus);
    
    MinVals = Complete(1).Plate16.TecanData.NormBYRatio(IdxMin);
    PlusVals = Complete(1).Plate16.TecanData.NormBYRatio(IdxPlus);
    
    [h,p,ci,stats] = ttest2(MinVals,PlusVals);
    
    Complete(1).Plate16.TecanData.CImin(IdxMin) = MinCI(i,1);
    Complete(1).Plate16.TecanData.CImax(IdxMin) = MinCI(i,2);
    Complete(1).Plate16.TecanData.CImin(IdxPlus) = PlusCI(i,1);
    Complete(1).Plate16.TecanData.CImax(IdxPlus) = PlusCI(i,2);
    
    Complete(1).Plate16.TecanData.pval(IdxAll) = p;
    
    for n = 1:numel(IdxAll)
        row = IdxAll(n);
        Complete(1).Plate16.TecanData.ttest_stats{row} = stats;
    end

    clear IdxMin IdxPlus IdxAll MinVals PlusVals
    
    
end


%% Anovas, Li ICQ colonies (1-way, between colonies)

% CI calculated via stat_summary in F5C plot code, copied into new variable
% "CI" for convenience
LW = [100 0];
Types = {'Cis','Trans'};
Days = 1:6;
Pairings = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'};

for l = 1:numel(LW)
    lw = LW(l);
    
    for t = 1:numel(Types)
        type = Types{t};
        
        for d = 1:numel(Days)
            day = Days(d);
            
            Idx = find(Complete(1).Plate11.ImageData.Perc_LW==lw &...
                Complete(1).Plate11.ImageData.Day==day & Complete(1).Plate11.ImageData.BactFcn~='None' &...
                Complete(1).Plate11.ImageData.YeastFcn~='None' & Complete(1).Plate11.ImageData.Type==type);
           
            ICQ = Complete(1).Plate11.ImageData.ICQ(Idx);
            Pairs = Complete(1).Plate11.ImageData.BactYeastPair(Idx);

            [p,table,stats] = anovan(ICQ,{Pairs},'sstype',2,'alpha',0.5);

            F = table{2,6};
            
            Complete(1).Plate11.ImageData.ICQ_AOV1_p(Idx) = p;
            Complete(1).Plate11.ImageData.ICQ_AOV1_F(Idx) = F;
            
            for i = 1:numel(Idx)
                row = Idx(i);
                Complete(1).Plate11.ImageData.ICQ_AOV1_stats{row} = stats;
            end
            
            if t==1
                for p = 1:numel(Pairings)
                    pair = Pairings{p};
                    
                    IdxP = find(Complete(1).Plate11.ImageData.Perc_LW==lw &...
                        Complete(1).Plate11.ImageData.Day==day & ...
                        Complete(1).Plate11.ImageData.BactYeastPair==pair &...
                        Complete(1).Plate11.ImageData.Type==type);

%                     Complete(1).Plate11.ImageData.CImin(IdxP) = CI.yci{p+4*(l-1),1}(d,1);
%                     Complete(1).Plate11.ImageData.CImax(IdxP) = CI.yci{p+4*(l-1),1}(d,2);
                end
            end
            clear Idx ICQ Pairs p table stats
            
        end
    end
end


%% Anovas, Li ICQ colonies (2-way, between colonies and perc lw)

% CI calculated via stat_summary in F5C plot code, copied into new variable
% "CI" for convenience
% LW = [100 0];
Types = {'Cis','Trans'};
Days = 1:6;
Pairings = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'};

    for t = 1:numel(Types)
        type = Types{t};
        
        for d = 1:numel(Days)
            day = Days(d);
            
            Idx = find(Complete(1).Plate11.ImageData.Day==day & Complete(1).Plate11.ImageData.BactFcn~='None' &...
                Complete(1).Plate11.ImageData.YeastFcn~='None' & Complete(1).Plate11.ImageData.Type==type);
           
            ICQ = Complete(1).Plate11.ImageData.ICQ(Idx);
            Pairs = Complete(1).Plate11.ImageData.BactYeastPair(Idx);
            LWs = Complete(1).Plate11.ImageData.Perc_LW(Idx);
            
            [p,table,stats] = anovan(ICQ,{Pairs,LWs},'sstype',2,'alpha',0.5);

                        
            for i = 1:numel(Idx)
                row = Idx(i);
                Complete(1).Plate11.ImageData.ICQ_AOV2_p{row} = p;
                Complete(1).Plate11.ImageData.ICQ_AOV2_stats{row} = stats;
            end
            
%             if t==1
%                 for p = 1:numel(Pairings)
%                     pair = Pairings{p};
%                     
%                     IdxP = find(Complete(1).Plate11.ImageData.Day==day & ...
%                         Complete(1).Plate11.ImageData.BactYeastPair==pair &...
%                         Complete(1).Plate11.ImageData.Type==type);
% 
%                     Complete(1).Plate11.ImageData.CImin(IdxP) = CI.yci{p+4*(l-1),1}(d,1);
%                     Complete(1).Plate11.ImageData.CImax(IdxP) = CI.yci{p+4*(l-1),1}(d,2);
%                 end
%             end
            clear Idx ICQ Pairs LWs p table stats
            
        end
    end

    
%% ttest, ICQ b/w LW% day 6

Pairs = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'};

for d = 1:6
for i = 1:numel(Pairs)
    pair = Pairs{i};
    
    IdxMin = find(Complete(1).Plate11.ImageData.BactYeastPair==pair &...
        Complete(1).Plate11.ImageData.Perc_LW==0 &...
        Complete(1).Plate11.ImageData.Day==d &...
        Complete(1).Plate11.ImageData.Type=='Cis');
    
    IdxPlus = find(Complete(1).Plate11.ImageData.Perc_LW==100 &...
        Complete(1).Plate11.ImageData.BactYeastPair==pair &...
        Complete(1).Plate11.ImageData.Day==d &...
        Complete(1).Plate11.ImageData.Type=='Cis');
    
    IdxAll = vertcat(IdxMin,IdxPlus);
    
    MinVals = Complete(1).Plate11.ImageData.ICQ(IdxMin);
    PlusVals = Complete(1).Plate11.ImageData.ICQ(IdxPlus);
    
    [h,p,ci,stats] = ttest2(MinVals,PlusVals);
    
%     Complete(1).Plate11.ImageData.CImin(IdxMin) = MinCI(i,1);
%     Complete(1).Plate11.ImageData.CImax(IdxMin) = MinCI(i,2);
%     Complete(1).Plate11.ImageData.CImin(IdxPlus) = PlusCI(i,1);
%     Complete(1).Plate11.ImageData.CImax(IdxPlus) = PlusCI(i,2);
%     
    Complete(1).Plate11.ImageData.LW_ttest_pval(IdxAll) = p;
    
    for n = 1:numel(IdxAll)
        row = IdxAll(n);
        Complete(1).Plate11.ImageData.ttest_stats{row} = stats;
    end

    clear IdxMin IdxPlus IdxAll MinVals PlusVals
    
    
end
end
    
Complete(1).Plate11.ImageData = movevars(Complete(1).Plate11.ImageData, 'LW_ttest_pval', 'Before', 'Type');
%% Error and significance, killer plate fold decrease


% CI calculated via stat_summary in F5C plot code, copied into new variable
% "CI" for convenience
W = [15 10 5 1];
Pairs = {'cutB_crispry','nocutB_crispry'};

for i = 1:numel(W)
    w = W(i);
    
    for day = 1:12
        
    IdxMin = find(Complete(1).Plate15.TecanEndpoints.PercW==w &...
        Complete(1).Plate15.TecanEndpoints.BactYeastPair=='nocutB_crisprY' &...
        Complete(1).Plate15.TecanEndpoints.day==day &...
        Complete(1).Plate15.TecanEndpoints.PercU==0);
    
    IdxPlus = find(Complete(1).Plate15.TecanEndpoints.PercW==w &...
        Complete(1).Plate15.TecanEndpoints.BactYeastPair=='cutB_crisprY' &...
        Complete(1).Plate15.TecanEndpoints.day==day &...
        Complete(1).Plate15.TecanEndpoints.PercU==0);
    
    IdxAll = vertcat(IdxMin,IdxPlus);
    
    MinVals = Complete(1).Plate15.TecanEndpoints.FoldDecreaseMono(IdxMin);
    PlusVals = Complete(1).Plate15.TecanEndpoints.FoldDecreaseMono(IdxPlus);
    
    [h,p] = ttest2(MinVals,PlusVals);
    
    Complete(1).Plate15.TecanEndpoints.FDM_pval(IdxAll) = p;

    Complete(1).Plate15.TecanEndpoints.FDM_CImin(IdxPlus) = CI.yci{i*2-1,1}(day,1);
    Complete(1).Plate15.TecanEndpoints.FDM_CImax(IdxPlus) = CI.yci{i*2-1,1}(day,2);
    Complete(1).Plate15.TecanEndpoints.FDM_CImin(IdxMin) = CI.yci{i*2,1}(day,1);
    Complete(1).Plate15.TecanEndpoints.FDM_CImax(IdxMin) = CI.yci{i*2,1}(day,2);
    
    
    clear IdxMin IdxPlus IdxAll MinVals PlusVals
    
    end
end

%% Add alt labels

pairs = {'crossB_crossY','crossB_wtY','wtB_crossY','wtB_wtY'};
altPairs = {'E_{cross} S_{cross}','E_{cross} S','E S_{cross}','E S'};

for p = 1:numel(pairs)
    pair = pairs{p};
    altPair = altPairs{p};
    
    Complete(1).All.AltNom(Complete(1).All.BactYeastPair==pair) = categorical(cellstr(altPair));
    
    
end

