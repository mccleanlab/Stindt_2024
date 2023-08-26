function dydt = Rescue_ODEs(t, y, params)

% What's different from Clump_ODEs?? 
%
% 1) Transconjugants are
% not subject to amino acid terms, so any yeast that recieve plasmid are
% able to grow normally regardless of global AA. 
%
% 2) Transconjugant growth rate modified to be somewhere in between
% clumped-yeast rate and free-yeast rate, since it's intuitively reasonable
% that rescued yeast wouldn't stay permanently clumped

% Basal growth rates
Rb = params(1);
Ry = params(2);

% Basal carrying capacities
Kb = params(3);
Ky = params(4);

% Incomplete niche overlap terms
Cb = params(5);
Cy = params(6);

% AA terms are designated by WHICH STRAIN NEEDS THEM, which varies per
% experiment (for Gy, could be W, U, or H)

% Global AA added to media
Gb = params(7);
Gy = params(8);

% AA secreted from opposite cell
Ab = params(9);
Ay = params(10);

% Monod term for needed AA
kb = params(11);
ky = params(12);

% Degradation of AA (just using one term)
D = params(13);

% TKC term (rate given non-clumped coincidence)
Gam = params(14);

% Basal death rates
Delb = params(15);
Dely = params(16);

% Multiplier for secretion of opposing cell when together in proximity
% (clumped)
Pb = params(17);
Py = params(18);

% Clump terms denoted as "L" to avoid confusion w/ niche overlap "C"

% Max of each cell type in clumps
KLb = params(3);
KLy = params(4);

% Rates of clumped cells
Rlb = params(19);
Rly = params(20);

% Rate of clumping per cell coincidence in well-mixed culture
Rl = params(21);

% TKC term for clumped cells, in 2 flavors: GamL is dependent on collisions
% of clumped B & Y, while GamC is dependent on number of clumps
GamL = params(22);
% GamC = params(22);


% Conversion factors for transforming cell collisions from respective
% fluorescence values (Cher & Cit) into cells, then back into Citrine for
% TKC tracking. From AllFits.Saved_Params.Fluor_Conversion.

Cher = .4437;
Cit = .0710;



dydt = zeros(6,1);

B = y(1,:);
Y = y(2,:);
T = y(3,:);
L = y(4,:); % Total clumps 
Lb = y(5,:); % Clumped bact
Ly = y(6,:); % Clumped yeast

% Bact, yeast
dydt(1) = (Rb .* B) .* ((Ab .* Y + Gb) ./ (Ab .* Y + Gb + kb .* D)) .* ...
    (1 - ((B + Lb) ./ Kb) - (Cy .* (T+Y+Ly) ./ Ky)) - (Delb .* B ./ Kb) -...
    Rl .* (B .* Y) ./ (B + Y) - (Rl .* (L .* B) ./ (KLb + Kb)); 
dydt(2) = (Ry .* Y) .* (Ay .* B + Gy) ./ (Ay .* B + Gy + ky .* D) .* ...
    (1 - ((Y + Ly+T) ./ Ky) - (Cb .* (B + Lb) ./ Kb)) - (Dely .* Y ./ Ky) -...
    Rl .* (B .* Y) ./ (B + Y) - (Rl .* (L .* Y) ./ (KLy + Ky));

% Transconjugants: Choose whether to use GamL or GamC (see above)
% dydt(3) = ((Rly .* T) .* (Ay .* B + Gy) ./ (Ay .* B + Gy + ky .* D) .* ...
%     (1 - (Y ./ Ky) - (Cb .* B ./ Kb))) - (Dely .* Y ./ Ky) + Gam .* (B .* Y) ./ (B + Y) +...
%     GamC .* L;
dydt(3) = (((Rly+Ry/4) .* T) .* ...
    (1 - ((Y+Ly+T) ./ Ky) - (Cb .* B ./ Kb))) - (Dely .* (T) ./ Ky) +...
    Cit * Gam .* ((B./Cher) .* (Y./Cit)) ./ (B./Cher + Y./Cit) +...
    Cit * GamL .* ((Lb./Cher) .* (Ly./Cit)) ./ (Lb./Cher + Ly./Cit);

% Clumps: tot, bact, yeast
dydt(4) = Rl .* (B .* Y) ./ (B + Y);
dydt(5) = (Rlb .* Lb) .* ((Pb * Ab .* Ly + Gb) ./ (Pb * Ab .* Ly + Gb + kb .* D)) .* ...
    (1 - ((Lb) ./ KLb) - (Cy .* (Ly) ./ KLy)) + (Rl.* (L .* B) ./ (L+B)) - (Delb .* Lb ./ KLb); 
dydt(6) = (Rly .* Ly) .* (Py * Ay .* Lb + Gb) ./ (Py * Ay .* Lb + Gb + kb .* D) .* ...
    (1 - ((Ly) ./ KLy) - (Cb .* (Lb) ./ KLb)) + (Rl .* (L .* Y) ./ (L + Y)) - (Dely .* Ly ./ KLy); 


end