%% Rates and Carrying Capacities: plot data and store fit info

fitdata = [];
fitdata0 = table();
close all
% exp = string(Experiments{8,1});
% date = double(Experiments{8,2});

% Enter in your search conditions here
plt = 'Plate2';
wellList = unique(AllData.(plt).well(AllData.(plt).BactYeastPair=='wtY'),'stable');
TimeBounds = [172*1000 255*1000]; %Useful: run "time intervals" section of TecAnalyzeSparkMulti
sp = categorical("Yeast");
% perc = 'L_Percent'; %This is the variable name for percent values in data
% LimAA = categorical("L");
Exp = categorical(AllData.(plt).Experiment(1)); %
channelList = {'OD600','Citrine'}; 

% Specifiy lower bounds (LBs), upper bounds (UBs) & init conditions (ICs) per channel (in order of P0 K r)
LBs = {[0.01, 0.01, 0.001];[0, 0, 0]};  
UBs = {[1, 1, 2];[10000, 400000000, 200000]}; 
ICs = {[0.5, 2, .5];[100, 10000, 5]}; 


% Cycle through samples
for c = 1:numel(channelList)
    channel = channelList{c};
    LB = LBs{c};
    UB = UBs{c};
    IC = ICs{c};
    figure
    for w = 1:numel(wellList)
        % Extract relevant growth data
        well = wellList(w);
        clear Idx
        Idx = find(AllData.(plt).well==well & AllData.(plt).time>TimeBounds(1) &...
            AllData.(plt).time<TimeBounds(2));
        t = AllData.(plt).time(Idx);
        t0 = min(t);
        t = t - t0;
        t = t./3600;
        y = AllData.(plt).(channel)(Idx);
%         y = log(y);

        % Fit growth curve using Logistic Growth
        eqn = '(K.*P0.*exp(r.*t))./(K+P0.*(exp(r.*t)-1))'; %K = asymptote, r = growth rate, P0 = init population
%         eqn = 'K./(1+((K-P0)./P0).*exp(-r.*t))'; %K = asymptote, r = growth rate, P0 = init population
        fo = fitoptions('Method','NonlinearLeastSquares');
        fn = fittype(eqn,'coefficients',{'P0','K','r'},'independent','t','options',fo);
        [f, gof] = fit(t,y,fn,'StartPoint',IC,'Lower',LB,'Upper',UB);%
        
        % Plot
        hold on
        p(1) = plot(t,y,'LineWidth',1);%,'color',[1 1-w/length(wellList) 0+w/length(wellList)]);
        title(channel)
        p(2) = plot(f,'--');
        set(p(2),'Color',[0.5 0.5 0.5],'LineWidth',0.25)
legend('off')        
hold off

%         Store growth curve parameters
        fitdata0 = [];
        fitdata0.well = well;
        fitdata0.channel = channel;
        fitdata0.fitEqn = string(formula(f));
        fitParams = coeffnames(f)';
        fitParamValues = coeffvalues(f);
        for i = 1:length(fitParams)
            fitdata0.(fitParams{i}) = fitParamValues(i);
        end
        gofParams = fieldnames(gof);
        for i = 1:length(gofParams)
            v = gofParams{i};
            fitdata0.(gofParams{i}) = gof.(v);
        end
        
        fitdata0.Experiment = Exp;
        fitdata0.Perc_AA = 100; %AllData.(plt).(perc)(Idx(1))
        fitdata0.LimAA = categorical("All"); %LimAA
        fitdata0.Mannose = categorical("N");
        if sp=='Bact'
        fitdata0.Strain = AllData.(plt).Bacteria(Idx(1));
        fitdata0.StrainFcn = AllData.(plt).BactFcn(Idx(1));
        elseif sp=='Yeast'
        fitdata0.Strain = AllData.(plt).Yeast(Idx(1));
        fitdata0.StrainFcn = AllData.(plt).YeastFcn(Idx(1));
        end
        fitdata0.Sp = sp;
        fitdata0.FitDay = AllData.(plt).day(Idx(1));

        fitdata = vertcat(fitdata,fitdata0);
    end
end

fitdata = struct2table(fitdata);
fitdata = movevars(fitdata, 'Experiment', 'Before', 'well');
%% Consolidate across experiments

% I usually just run this in command window whenever I feel fitdata looks
% good enough
CompiledFits = vertcat(CompiledFits,fitdata);

%% Throw into master table

% Ditto above
AllFits.RatesCapacities = CompiledFits;

%% Plot anything for checking

clear g;
figure();

g(1,1) = gramm('x',CompiledFits.Perc_AA,'y',CompiledFits.r,...
    'color',CompiledFits.Strain,'subset',...
    (CompiledFits.Strain=='WT M2' & CompiledFits.Perc_AA==100 & CompiledFits.rsquare>0.9));
g(1,1).facet_grid([],CompiledFits.channel);
g(1,1).stat_boxplot();
g(1,1).geom_point('dodge',1);
g(1,1).set_title('Growth rate');
g(1,1).axe_property('XTick',[]);
g(1,1).set_names('x','','y','Fit rate','Color','Strain');
% g(1,1).set_order_options('row',[15 10 5 1],'column',[100 5],...
%     'linestyle',{'crisprY','None'},'size',{'None','crisprY'});
% g(1,1).set_line_options('base_size',0.5,'step_size',2,'styles',{'-',':'});
% g(1,1).set_color_options('map',[1 0 0],'legend','expand');
% g(1,1).no_legend();

g.set_title('Fit terms');

g.draw();

%% Caculate Average r and K terms from 100% AA growths

AllFits.Summary = table();
Summary = [];
Summary0 = table();

Mets = AllFits.RatesCapacities.Properties.VariableNames;

Idx = find(AllFits.RatesCapacities.Experiment=='230321_Fits');

Strains = unique(AllFits.RatesCapacities.Strain(Idx));
% Fcns = unique(AllFits.RatesCapacities.StrainFcn(Idx));
% Spp = unique(AllFits.RatesCapacities.Sp(Idx));

for i = 1:numel(Strains)
    s = Strains(i);
                        
            Chans = unique(AllFits.RatesCapacities.Channel(AllFits.RatesCapacities.Strain==s));
            
            for l = 1:2
                c = Chans(l);
                
                IdxStrain = find(AllFits.RatesCapacities.Strain==s & AllFits.RatesCapacities.Experiment=='230321_Fits');
                
                clear Summary0
                
                Summary0.Strain = s;
                Summary0.StrainFcn = AllFits.RatesCapacities.StrainFcn(IdxStrain(1));
                Summary0.Spp = AllFits.RatesCapacities.Sp(IdxStrain(1));
                Summary0.Channel = c;
                Summary0.Perc_AA = 100;
                Summary0.LimAA = 100;
                Summary0.FitDay = 2;
                Summary0.fitEqn = AllFits.RatesCapacities.fitEqn(1);
           
                IdxMean = find(AllFits.RatesCapacities.Strain==s & AllFits.RatesCapacities.rsquare>0.85 &...
                    AllFits.RatesCapacities.Channel==c & AllFits.RatesCapacities.Experiment=='230321_Fits');
              
                Summary0.r = nanmean(AllFits.RatesCapacities.r(IdxMean));
                Summary0.K = nanmean(AllFits.RatesCapacities.K(IdxMean));
                Summary0.P0 = nanmean(AllFits.RatesCapacities.P0(IdxMean));
            
                Summary = vertcat(Summary,Summary0);
                
            end
end

AllFits.Summary = struct2table(Summary);


%% Solve expanded eqtns with competitive terms (no crossfeeding)

syms y(t) b(t) 
syms Cb Cy

Rb = AllFits.Summary.r(AllFits.Summary.Strain=='WT M2' & AllFits.Summary.Channel=='Cherry'); 
Kb = AllFits.Summary.K(AllFits.Summary.Strain=='WT M2' & AllFits.Summary.Channel=='Cherry'); 

Ry = AllFits.Summary.r(AllFits.Summary.Strain=='yMM1636' & AllFits.Summary.Channel=='Citrine'); 
Ky = AllFits.Summary.K(AllFits.Summary.Strain=='yMM1636' & AllFits.Summary.Channel=='Citrine'); 

% eqns = [diff(y,t) == Ry * y*(1 - (b*Cb)/Kb - y/Ky), diff(b,t) == Rb * b*(1 - b/Kb - (y*Cy)/Ky)];
eqn = diff(y,t) == Ry * y*(1 - (b*Cb)/Kb - y/Ky)
Sy = dsolve(eqn,'Implicit',true)

%% Try using embedded ODE solver and lsqcurvefit

% from this post: https://www.mathworks.com/matlabcentral/answers/43439-monod-kinetics-and-curve-fitting#comment_89455

Rb = AllFits.Summary.r(AllFits.Summary.Strain=='WT M2' & AllFits.Summary.Channel=='Cherry'); 
Kb = AllFits.Summary.K(AllFits.Summary.Strain=='WT M2' & AllFits.Summary.Channel=='Cherry'); 

Ry = AllFits.Summary.r(AllFits.Summary.Strain=='yMM1636' & AllFits.Summary.Channel=='Citrine'); 
Ky = AllFits.Summary.K(AllFits.Summary.Strain=='yMM1636' & AllFits.Summary.Channel=='Citrine'); 

B = AllData.Plate2.Cherry(AllData.Plate2.time>256000 & AllData.Plate2.time<340000);
Y = AllData.Plate2.Citrine(AllData.Plate2.time>256000 & AllData.Plate2.time<340000);

time = AllData.Plate2.time(AllData.Plate2.time>256000 & AllData.Plate2.time<340000);

[Sol,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqcurvefit(@TKC_ODEs1,Rb,Ry,Kb,Ky,B,Y,Time);
