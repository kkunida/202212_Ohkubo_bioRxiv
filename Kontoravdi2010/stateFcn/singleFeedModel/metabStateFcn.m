function dMetabdt = metabStateFcn(V,Xv,Metab,Q_in,mu)
%% parameters
Kdeg_Gln = 9.6E-3;
mnt_Glc = 4.9E-14;
Y_Amm_Gln = 0.45;
Y_Lac_Glc = 2;
Y_X_Glc = 2.6E+8;
Y_X_Gln = 8E+8;
alpha1 = 3.4E-13;
alpha2 = 4;

Metab_in = [500; 100; 0; 0;];

%% state
% Glc = Metab(1);
Gln = Metab(2);
% Lac = Metab(3);
% Amm = Metab(4);

%% production rate
mnt_Gln = (alpha1*Gln)/(alpha2 + Gln);

q_Glc = -mu/Y_X_Glc - mnt_Glc;
q_Gln = -mu/Y_X_Gln - mnt_Gln;
q_Lac = -Y_Lac_Glc*q_Glc;
q_Amm = -Y_Amm_Gln*q_Gln;

q_Metab = [q_Glc;q_Gln;q_Lac;q_Amm];
Fdeg = Kdeg_Gln*Gln*[0;-1;0;1]; % Gln is degraded to Amm

%% ODEs
dMetabdt = q_Metab.*Xv - (Q_in/V).*(Metab - Metab_in) + Fdeg; 

end