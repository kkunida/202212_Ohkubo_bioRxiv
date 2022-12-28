function dxdt = multiFeedStateFcnError(x,u)
%% states
V = x(1);
Xv = x(2);
Xt = x(3);
Metab = x(4:7);
mAb_ext = x(8);
mRNA = x(9:10);
Assy = x(11:15);
mAb_gol = x(16);
% V_Glc = x(17);
% V_Gln = x(18);

%% inputs
Q_Glc = u(1);
Q_Gln = u(2);
Q_in = Q_Glc + Q_Gln;

%% ODEs
[dVdt,dV_Glcdt,dV_Glndt,dXvdt,dXtdt,mu] = multiFeedGrowthStateFcnError(V,Xv,Xt,Metab,Q_Glc,Q_Gln);
dMetabdt = multiFeedMetabStateFcnError(V,Xv,Metab,Q_Glc,Q_Gln,mu);
[dmAb_extdt,dmRNAdt,dAssydt,dmAb_goldt] = mAbStateFcn(V,Xv,mAb_ext,mRNA,Assy,mAb_gol,Q_in,mu);

%% outputs
dxdt = zeros(18,1);
dxdt(1) = dVdt;
dxdt(2) = dXvdt;
dxdt(3) = dXtdt;
dxdt(4:7) = dMetabdt;
dxdt(8) = dmAb_extdt;
dxdt(9:10) = dmRNAdt;
dxdt(11:15) = dAssydt;
dxdt(16) = dmAb_goldt;
dxdt(17) = dV_Glcdt;
dxdt(18) = dV_Glndt;

% lower limit of states
dxdt(x < 0 & dxdt < 0) = 0;

end