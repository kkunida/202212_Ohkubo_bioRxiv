function [dVdt,dV_Glcdt,dV_Glndt,dXvdt,dXtdt,mu] = multiFeedGrowthStateFcn(V,Xv,Xt,Metab,Q_Glc,Q_Gln)
%% parameters
Kdeath_Amm = 1.76;
K_Glc = 0.75;
K_Gln = 0.075;
Ki_Amm = 28.48;
Ki_Lac = 171.76;
n = 2;
mu_max = 0.058;
mu_death_max = 0.03; % 0.06 in [Kontoravdi+,2010]
Klysis = 0.03;

%% inputs
Q_in = Q_Glc + Q_Gln; % total feed rate

%% state variables
Glc = Metab(1);
Gln = Metab(2);
Lac = Metab(3);
Amm = Metab(4);

%% growth and death rate
f_lim = (Glc/(K_Glc + Glc))*(Gln/(K_Gln + Gln));
f_inh = (Ki_Lac/(Ki_Lac + Lac))*(Ki_Amm/(Ki_Amm + Amm));
mu = mu_max*f_lim*f_inh;
mu_death = mu_death_max/(1 + (Kdeath_Amm/Amm)^n);

%% ODEs
dVdt = Q_in;
dV_Glcdt = Q_Glc;
dV_Glndt = Q_Gln;
dXvdt = (mu - mu_death)*Xv - (Q_in/V)*Xv;
dXtdt = mu*Xv - (Q_in/V)*Xt - Klysis*(Xt - Xv);

end