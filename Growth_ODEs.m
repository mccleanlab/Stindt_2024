function [dydt dbdt] = Growth_ODEs(t,y,b)

% Global Parameters
d3 = 0.08; % Turn off if modeling mRNA decay
k4 = 15;
% d4 = 0.003;
d4 = 0.001;
k5 = 0.06;

% Promoter specific parameters
k1 = promoter_params(1);
d1 = promoter_params(2);
k2 = promoter_params(3);
K = promoter_params(4);
n = promoter_params(5);
d2 = promoter_params(6);
k3 = promoter_params(7);

% ODE indices
dydt = zeros(6,1);
P_unbound = y(1);
P_bound = y(2);
P_active = y(3);
mRNA = y(4);
YFP = y(5);
mYFP = y(6);

% Basic equations
K = K_scale*K;
dydt(1) = d1*P_bound - ((k1*Msn2^n)/(K^n + Msn2^n))*P_unbound;
dydt(2) = ((k1*Msn2^n)/(K^n + Msn2^n))*P_unbound + d2*P_active - d1*P_bound - k2*Msn2*P_bound;
dydt(3) = k2*Msn2*P_bound - d2*P_active;
dydt(4) = k3*P_active - d3*mRNA;
dydt(5) = k4*mRNA - d4*YFP - k5*YFP;
dydt(6) = k5*YFP - d4*mYFP;

