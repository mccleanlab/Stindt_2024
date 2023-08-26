%% Notes on Amino Acid stuff, copied from Mathematica script:

% For production (benefit) terms, the benefit is really just per-cell
% secretion of limited amino acid, when you look at eqtns:
%
% Zhao 2011 has Trp secretion at 0.022 gTrp/gCDW*hr (that's cell dry
% weight/hr), which, for an individual CDW of 3e-13 g for E. coli (ECMDB
% stats: ecmdb.ca/e_coli_stats), comes out to 1.95e7 molecules/cell*hr
% (=10^7.29 molecules/200uL) Muller 2014 estimates Leu secretion to be
% 7e5molecules/cell*s = 2.52e9 molecules/cell*hr (=10^9.40 molecules/200uL)
% 
% Monod terms ("k") split into kX, where X is amino acid specific, based on
% the relevant strain that needs it (bact require L, yeast require W, U,
% and/or H). Leu degradation (k-1 = 0.03) taken from Kobayashi 2010, from
% Fig 3 (at 30°C, 10^3/T(Kelvin) = 3.3). Assuming the same for Trp for now, since I can't find obvious
% candidate figures. Both need more research.
% 
% Finally, this model includes formation- and growth of transconjugants, T.
% The formation of these is based on reaction-diffusion coincidence of
% bacteria B and yeast Y, with a constant transfer rate, \[Gamma], that is
% thus only dependent on the concentrations of B and Y. T grows with the
% same parameters as Y (i.e. there's assumed no fitness cost). Because of
% the inclusion of T, all terms that are Y-dependent are converted to (Y +
% T)-dependence, e.g. both strains produce limited amino acid and are
% summed in production term; the strains reach carrying capacity together.
% Parameter value for \[Gamma] based on Volkova et. al.

% First just copy KS_Solution into Saved_Params as "g", based on equations
AllFits.Saved_Params.g = AllFits.KS_Solution;

% Save degradation terms
AllFits.Saved_Params.d = table();
AllFits.Saved_Params.d.L = 0.03;
AllFits.Saved_Params.d.W = 0.03;
AllFits.Saved_Params.d.H = 0.03;
AllFits.Saved_Params.d.U = 0.03;
AllFits.Saved_Params.d.Properties.RowNames = "1/s";

% Crossfeeder info (secretion rates)
AllFits.Saved_Params.a = table();
AllFits.Saved_Params.a.L(2) = 1.95e7;
AllFits.Saved_Params.a.W(2) = 2.52e9;
AllFits.Saved_Params.a.L(1) = (AllFits.Saved_Params.a.L(2) * 1E6)/(200 * 6.0221408e23);
AllFits.Saved_Params.a.W(1) = (AllFits.Saved_Params.a.W(2) * 1E6)/(200 * 6.0221408e23);
AllFits.Saved_Params.a.Properties.RowNames = ["Molar/Cell*hr","Molecules/Cell*hr"];