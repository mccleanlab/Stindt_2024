%% Save molarity of each AA limited in experiments

AllFits.KS_Solution = table();

% First make variable that includes 100% molarity values of each AA modified in experiments,
% Based on https://openwetware.org/wiki/McClean:KS_Amino_Acid_Supplement
AllFits.KS_Solution.L(1) = (100/1000)/131.2;
AllFits.KS_Solution.W(1) = (50/1000)/204.2;
AllFits.KS_Solution.H(1) = (20/1000)/209.6;
AllFits.KS_Solution.U(1) = (20/1000)/112.1;

% Also put those molarities in terms of molecules per well, i.e. 200uL
AAs = {'L','W','H','U'};
for i = 1:4
    AA = AAs{i};
    AllFits.KS_Solution.(AA)(2) = AllFits.KS_Solution.(AA)(1) * 6.0221408e23 * 0.5 / 1000;
end
AllFits.KS_Solution.Properties.RowNames = ["Molar","Molecules_200uL"];

%% Then use that to calculate molarity of limited AAs in experiments

% If it's limited, calculate...
AllFits.R_K_k.L_Molar(AllFits.R_K_k.LimAA=="L" | AllFits.R_K_k.LimAA=="LW") = ...
    (AllFits.R_K_k.Perc_AA(AllFits.R_K_k.LimAA=="L"| AllFits.R_K_k.LimAA=="LW").*AllFits.KS_Solution.L(1))./100;
% Or else set it to 100%
AllFits.R_K_k.L_Molar(AllFits.R_K_k.LimAA~="LW" & AllFits.R_K_k.LimAA~="L") = AllFits.KS_Solution.L(1); 

%Repeat for other AAs
AllFits.R_K_k.W_Molar(AllFits.R_K_k.LimAA=="W" | AllFits.R_K_k.LimAA=="LW") = ...
(AllFits.R_K_k.Perc_AA(AllFits.R_K_k.LimAA=="W" | AllFits.R_K_k.LimAA=="LW").*AllFits.KS_Solution.W(1))./100;
AllFits.R_K_k.W_Molar(AllFits.R_K_k.LimAA~="W" & AllFits.R_K_k.LimAA~="LW") = AllFits.KS_Solution.W(1);

AllFits.R_K_k.H_Molar(AllFits.R_K_k.LimAA=="H" | AllFits.R_K_k.LimAA=="UH") = ...
(AllFits.R_K_k.Perc_AA(AllFits.R_K_k.LimAA=="H" | AllFits.R_K_k.LimAA=="UH").*AllFits.KS_Solution.H(1))./100;
AllFits.R_K_k.H_Molar(AllFits.R_K_k.LimAA~="H" & AllFits.R_K_k.LimAA~="UH") = AllFits.KS_Solution.H(1);

AllFits.R_K_k.U_Molar(AllFits.R_K_k.LimAA=="UH") = (AllFits.R_K_k.Perc_AA(AllFits.R_K_k.LimAA=="UH").*AllFits.KS_Solution.U(1))./100;
AllFits.R_K_k.U_Molar(AllFits.R_K_k.LimAA~="UH") = AllFits.KS_Solution.U(1);

%% Generate Rmax per each strain, using 100% AA data

strains = unique(AllFits.R_K_k.Strain);

for i = 1:numel(strains)

strain = strains(i);

sp = unique(AllFits.R_K_k.Sp(AllFits.R_K_k.Strain==strain));

    if sp=="Yeast"
    chans = {'OD600','Citrine'};
    elseif sp=="Bact"
    chans = {'OD600','Cherry'};
    end

    for j = 1:2

    chan = chans{j};
    
    Rmax = mean(AllFits.R_K_k.r(AllFits.R_K_k.Strain==strain & AllFits.R_K_k.Channel==chan &...
        AllFits.R_K_k.Perc_AA==100 & AllFits.R_K_k.rsquare >= 0.85));
    AllFits.R_K_k.Rmax(AllFits.R_K_k.Strain==strain & AllFits.R_K_k.Channel==chan) = Rmax;

    end

end

%% Go ahead and just calculate num molecules from these molar values

mol = {'L_Molar','W_Molar','H_Molar','U_Molar'};
num = {'L_Num','W_Num','H_Num','U_Num'};

for i = 1:4
m = mol{i};
n = num{i};

AllFits.R_K_k.(n) = AllFits.R_K_k.(m) .* 6.0221408e23 .* 0.5 ./ 1000;
end

%% Use Rmax values to find Monod constant k (Molar)

% Based on R = Rmax (A / (A + k)), rearranged as k = (A * Rmax / R) - A,
% where A is limited amino acid

% "crossB" are leu auxotrophs...
idx = find(AllFits.R_K_k.StrainFcn=='crossB');
AllFits.R_K_k.kL(idx) = (AllFits.R_K_k.L_Molar(idx) .* AllFits.R_K_k.Rmax(idx) ./ AllFits.R_K_k.r(idx)) - AllFits.R_K_k.L_Molar(idx);

% "crossY" are trp auxotrophs...
clear idx
idx = find(AllFits.R_K_k.StrainFcn=='crossY');
AllFits.R_K_k.kW(idx) = (AllFits.R_K_k.W_Molar(idx) .* AllFits.R_K_k.Rmax(idx) ./ AllFits.R_K_k.r(idx)) - AllFits.R_K_k.W_Molar(idx);

% All yeast are either auxotrophic for U or H, too, depending on strain
idxU = find(AllFits.R_K_k.Strain=='yMM1585');
AllFits.R_K_k.kU(idxU) = (AllFits.R_K_k.U_Molar(idxU) .* AllFits.R_K_k.Rmax(idxU) ./ AllFits.R_K_k.r(idxU)) - AllFits.R_K_k.U_Molar(idxU);

idxH = find(AllFits.R_K_k.Strain=='yMM1636' | AllFits.R_K_k.Strain=='yMM1720');
AllFits.R_K_k.kH(idxU) = (AllFits.R_K_k.H_Molar(idxU) .* AllFits.R_K_k.Rmax(idxU) ./ AllFits.R_K_k.r(idxU)) - AllFits.R_K_k.H_Molar(idxU);

% Set negative values to 0
AllFits.R_K_k.kL(AllFits.R_K_k.kL<0) = 0;
AllFits.R_K_k.kW(AllFits.R_K_k.kL<0) = 0;
AllFits.R_K_k.kH(AllFits.R_K_k.kL<0) = 0;
AllFits.R_K_k.kU(AllFits.R_K_k.kL<0) = 0;

% Set irrelevant values to nan
AllFits.R_K_k.kL(AllFits.R_K_k.L_Molar==AllFits.KS_Solution.L) = nan;
AllFits.R_K_k.kW(AllFits.R_K_k.W_Molar==AllFits.KS_Solution.W) = nan;
AllFits.R_K_k.kU(AllFits.R_K_k.U_Molar==AllFits.KS_Solution.U) = nan;
AllFits.R_K_k.kH(AllFits.R_K_k.H_Molar==AllFits.KS_Solution.H) = nan;

%% Use Rmax values to find Monod constant k (Numbers, to make sure they're the same)

% Based on R = Rmax (A / (A + k)), rearranged as k = (A * Rmax / R) - A,
% where A is limited amino acid

% "crossB" are leu auxotrophs...
idx = find(AllFits.R_K_k.StrainFcn=='crossB');
AllFits.R_K_k.kL_Num(idx) = (AllFits.R_K_k.L_Num(idx) .* AllFits.R_K_k.Rmax(idx) ./ AllFits.R_K_k.r(idx)) - AllFits.R_K_k.L_Num(idx);

% "crossY" are trp auxotrophs...
clear idx
idx = find(AllFits.R_K_k.StrainFcn=='crossY');
AllFits.R_K_k.kW_Num(idx) = (AllFits.R_K_k.W_Num(idx) .* AllFits.R_K_k.Rmax(idx) ./ AllFits.R_K_k.r(idx)) - AllFits.R_K_k.W_Num(idx);

% All yeast are either auxotrophic for U or H, too, depending on strain
idxU = find(AllFits.R_K_k.Strain=='yMM1585');
AllFits.R_K_k.kU_Num(idxU) = (AllFits.R_K_k.U_Num(idxU) .* AllFits.R_K_k.Rmax(idxU) ./ AllFits.R_K_k.r(idxU)) - AllFits.R_K_k.U_Num(idxU);

idxH = find(AllFits.R_K_k.Strain=='yMM1636' | AllFits.R_K_k.Strain=='yMM1720');
AllFits.R_K_k.kH_Num(idxU) = (AllFits.R_K_k.H_Num(idxU) .* AllFits.R_K_k.Rmax(idxU) ./ AllFits.R_K_k.r(idxU)) - AllFits.R_K_k.H_Num(idxU);

% Set negative values to 0
AllFits.R_K_k.kL_Num(AllFits.R_K_k.kL_Num<0) = 0;
AllFits.R_K_k.kW_Num(AllFits.R_K_k.kL_Num<0) = 0;
AllFits.R_K_k.kH_Num(AllFits.R_K_k.kL_Num<0) = 0;
AllFits.R_K_k.kU_Num(AllFits.R_K_k.kL_Num<0) = 0;

% Set irrelevant values to nan
AllFits.R_K_k.kL_Num(AllFits.R_K_k.L_Molar==AllFits.KS_Solution.L(1)) = nan;
AllFits.R_K_k.kW_Num(AllFits.R_K_k.W_Molar==AllFits.KS_Solution.W(1)) = nan;
AllFits.R_K_k.kU_Num(AllFits.R_K_k.U_Molar==AllFits.KS_Solution.U(1)) = nan;
AllFits.R_K_k.kH_Num(AllFits.R_K_k.H_Molar==AllFits.KS_Solution.H(1)) = nan;


%% Take a look at those k values

ks = {'kL','kW','kU','kH'};
kNs = {'kL_Num','kW_Num','kU_Num','kH_Num'};

for i = 1:4

k = ks{i};
kN = kNs{i};

subplot(2,4,i)

histogram(AllFits.R_K_k.(k))
meanAll = nanmean(AllFits.R_K_k.(k));
meanNonZero = nanmean(AllFits.R_K_k.(k)(AllFits.R_K_k.(k)>0));
xline(meanAll,'Color','b','label',['Mean of all = ' num2str(meanAll)]);
xline(meanNonZero,'Color','g','label',['Mean of nonzero values = ' num2str(meanNonZero)]);

tit = ['Computed values for ' k];
title(tit,'Interpreter','none');

subplot(2,4,i+4)

histogram(AllFits.R_K_k.(kN))
meanAll = nanmean(AllFits.R_K_k.(kN));
meanNonZero = nanmean(AllFits.R_K_k.(kN)(AllFits.R_K_k.(kN)>0));
xline(meanAll,'Color','b','label',['Mean of all = ' num2str(meanAll)]);
xline(meanNonZero,'Color','g','label',['Mean of nonzero values = ' num2str(meanNonZero)]);

tit = ['Computed values for ' kN];
title(tit,'Interpreter','none');


end

%% Also look at R and Rmax values

strains = unique(AllFits.R_K_k.Strain(AllFits.R_K_k.Perc_AA==100));

for i = 1:numel(strains)

strain = strains(i);
edges = linspace(0,1,10);
subplot(2,numel(strains),i)

h1 = histogram(AllFits.R_K_k.r(AllFits.R_K_k.Strain==strain & AllFits.R_K_k.Perc_AA==100),'BinEdges',edges);
% grid on
% xlim([0, 1]);
meanMax = nanmean(AllFits.R_K_k.r(AllFits.R_K_k.Strain==strain & AllFits.R_K_k.Perc_AA==100));
xline(meanMax,'Color','b','label',['Mean Rmax = ' num2str(meanMax)]);

tit = ['Rmax values for ' strain];
title(tit,'Interpreter','none');

subplot(2,numel(strains),i+numel(strains))

h2 = histogram(AllFits.R_K_k.r(AllFits.R_K_k.Strain==strain & AllFits.R_K_k.Perc_AA~=100),'BinEdges',edges);
% xlim([0, 1]);
meanRs = nanmean(AllFits.R_K_k.r(AllFits.R_K_k.Strain==strain & AllFits.R_K_k.Perc_AA~=100));
xline(meanRs,'Color','b','label',['Mean R (AA<100%) = ' num2str(meanRs)]);

tit = ["R values", "(AA<100%) for " strain];
title(tit,'Interpreter','none');

end

%% Save out k's

AllFits.Saved_Params.k.kL(1) = nanmean(AllFits.R_K_k.kL(AllFits.R_K_k.kL>0));
AllFits.Saved_Params.k.kW(1) = nanmean(AllFits.R_K_k.kW(AllFits.R_K_k.kW>0));
AllFits.Saved_Params.k.kU(1) = nanmean(AllFits.R_K_k.kU(AllFits.R_K_k.kU>0));
AllFits.Saved_Params.k.kH(1) = nanmean(AllFits.R_K_k.kH(AllFits.R_K_k.kH>0));

AllFits.Saved_Params.k.kL(2) = nanmean(AllFits.R_K_k.kL_Num(AllFits.R_K_k.kL>0));
AllFits.Saved_Params.k.kW(2) = nanmean(AllFits.R_K_k.kW_Num(AllFits.R_K_k.kW>0));
AllFits.Saved_Params.k.kU(2) = nanmean(AllFits.R_K_k.kU_Num(AllFits.R_K_k.kU>0));
AllFits.Saved_Params.k.kH(2) = nanmean(AllFits.R_K_k.kH_Num(AllFits.R_K_k.kH>0));

AllFits.Saved_Params.k.Properties.RowNames = ["Molar","Molecules_200uL"];