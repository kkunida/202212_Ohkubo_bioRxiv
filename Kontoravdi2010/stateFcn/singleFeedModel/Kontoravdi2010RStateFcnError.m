function dxdt = Kontoravdi2010RStateFcnError(x,Q_in)
%% states
V = x(1);
Xv = x(2);
Xt = x(3);
Metab = x(4:7);
mAb_ext = x(8);

%% ODEs
[dVdt,dXvdt,dXtdt,mu] = growthStateFcnError(V,Xv,Xt,Metab,Q_in);
dMetabdt = metabStateFcnError(V,Xv,Metab,Q_in,mu);

%% mAb production submodel
% parameters
Y_mAb_X = 3.0583E-12;
mnt_mAb = 1.1989E-08;

q_mAb = Y_mAb_X*mu + mnt_mAb;
dmAb_extdt = q_mAb*Xv - (Q_in/V)*mAb_ext;

%% outputs
dxdt = zeros(8,1);
dxdt(1) = dVdt;
dxdt(2) = dXvdt;
dxdt(3) = dXtdt;
dxdt(4:7) = dMetabdt;
dxdt(8) = dmAb_extdt;

% lower limit of states
dxdt(x < 0 & dxdt < 0) = 0;

end