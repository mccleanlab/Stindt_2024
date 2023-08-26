function dydt = Conjugation_ODEs_wDeath(t, y, params)

Rb = params(1);
Ry = params(2);
Kb = params(3);
Ky = params(4);
Cb = params(5);
Cy = params(6);

% AA terms are designated by WHICH STRAIN NEEDS THEM, which varies per
% experiment (for Gy, could be W, U, or H)
Gb = params(7);
Gy = params(8);
Ab = params(9);
Ay = params(10);
kb = params(11);
ky = params(12);

D = params(13);
Gam = params(14);
Delb = params(15);
Dely = params(16);

% Conversion factors for transforming cell collisions from respective
% fluorescence values (Cher & Cit) into cells, then back into Citrine for
% TKC tracking. From AllFits.Saved_Params.Fluor_Conversion.

Cher = .4437;
Cit = .0710;

dydt = zeros(3,1);

B = y(1,:);
Y = y(2,:);
T = y(3,:);

dydt(1) = ((Rb .* B) .* (Ab .* Y + Gb) ./ (Ab .* Y + Gb + kb .* D) .* ...
    (1 - (B ./ Kb) - (Cy .* Y ./ Ky))) - (Delb .* B ./ Kb); 
dydt(2) = ((Ry .* Y) .* (Ay .* B + Gy) ./ (Ay .* B + Gy + ky .* D) .* ...
    (1 - (Y ./ Ky) - (Cb .* B ./ Kb))) - (Dely .* Y ./ Ky);
dydt(3) = ((Ry .* T) .* (Ay .* B + Gy) ./ (Ay .* B + Gy + ky .* D) .* ...
    (1 - ((Y) ./ Ky) - (Cb .* B ./ Kb))) - (Dely .* T ./ Ky) +...
    Cit * Gam .* ((B./Cher) .* ((Y-T)./Cit)) ./ (B./Cher + Y./Cit);

end
