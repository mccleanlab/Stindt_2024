function dydt = Competition_ODEs(t, y, params)

Rb = params(1);
Ry = params(2);
Kb = params(3);
Ky = params(4);
Cb = params(5);
Cy = params(6);

dydt = zeros(2,1);
B = y(1,:);
Y = y(2,:);

dydt(1) = (Rb .* B) .* (1 - (B ./ Kb) - (Cy .* Y ./ Ky));
dydt(2) = (Ry .* Y) .* (1 - (Y ./ Ky) - (Cb .* B ./ Kb));

end