function dxdt = multiFeedRStateFcn(x,u)
%% states
V = x(1);
Xv = x(2);
Xt = x(3);
Metab = x(4:7);
mAb_ext = x(8);
% V_Glc = x(17);
% V_Gln = x(18);

%% inputs
Q_Glc = u(1);
Q_Gln = u(2);
Q_in = Q_Glc + Q_Gln;

%% ODEs
[dVdt,dV_Glcdt,dV_Glndt,dXvdt,dXtdt,mu] = multiFeedGrowthStateFcn(V,Xv,Xt,Metab,Q_Glc,Q_Gln);
dMetabdt = multiFeedMetabStateFcn(V,Xv,Metab,Q_Glc,Q_Gln,mu);

%% mAb production submodel
% parameters
Y_mAb_X = 3.0583E-12;
mnt_mAb = 1.1989E-08;

q_mAb = Y_mAb_X*mu + mnt_mAb;
dmAb_extdt = q_mAb*Xv - (Q_in/V)*mAb_ext;

%% outputs
dxdt = zeros(10,1);
dxdt(1) = dVdt;
dxdt(2) = dXvdt;
dxdt(3) = dXtdt;
dxdt(4:7) = dMetabdt;
dxdt(8) = dmAb_extdt;
dxdt(9) = dV_Glcdt;
dxdt(10) = dV_Glndt;

% lower limit of states
dxdt(x < 0 & dxdt < 0) = 0;

end