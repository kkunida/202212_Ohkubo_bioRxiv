function [dmAb_extdt,dmRNAdt,dAssydt,dmAb_goldt] = mAbStateFcn(V,Xv,mAb_ext,mRNA,Assy,mAb_gol,Q_in,mu)
%% parameter
Kdecay = 0.1;
K_Assy = 1E-6;
K_ER2gol = 0.69;
K_gol2ext = 0.14;
copyNum_Hc = 100;
copyNum_Lc = 100;
kTS_Hc = 3000;
kTS_Lc = 4500;
kTL_Hc = 17;
kTL_Lc = 11.5;
gamma1 = 0.1;
gamma2 = 2;
epsilon1 = 0.995;
epsilon2 = 1;
mass_mAb = 2.4908084E-16; % [mg/molecule] = 150 kDa

copyNum = [copyNum_Hc; copyNum_Lc];
kTS = [kTS_Hc; kTS_Lc];
kTL = [kTL_Hc; kTL_Lc];

%% state
Hc = Assy(1);
Lc = Assy(2);
Hc2 = Assy(3);
Hc2Lc = Assy(4);
mAb_ER = Assy(5);

%% reaction rate
r = zeros(3,1);
r(1) = (1/3)*K_Assy*Hc^2;
r(2) = 2*K_Assy*Hc2*Lc;
r(3) = K_Assy*Hc2Lc*Lc;

F_TS = kTS.*copyNum;
F_TL = kTL.*mRNA;

F_in_Assy = [F_TL;0;0;0];
F_ER2gol = K_ER2gol*mAb_ER;
F_out_Assy = [0;0;0;0;F_ER2gol];

F_gol2ext = K_gol2ext*mAb_gol;
q_mAb = mass_mAb*epsilon2*F_gol2ext;

%% stoichiometric matrix
S_Assy = [-2  0  0;     % Hc
           0 -1 -1;     % Lc
           1 -1  0;     % Hc2
           0  1 -1;     % Hc2Lc    
           0  0  1;     % Hc2Lc2 (mAb_ER)
          ];

%% ODEs
dmRNAdt = F_TS - Kdecay.*mRNA;
dAssydt = S_Assy*r + F_in_Assy - F_out_Assy;
dAssydt(2) = dAssydt(2) - 0.1*Lc; % modified from [Kontoravdi+,2010]
dmAb_goldt = epsilon1*F_ER2gol - F_gol2ext;
dmAb_extdt = (gamma2 - gamma1*mu)*q_mAb*Xv - (Q_in/V)*mAb_ext;

end