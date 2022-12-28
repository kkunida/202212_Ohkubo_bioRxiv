function dxdt = Kontoravdi2010StateFcn(x,Q_in)
%% states
V = x(1);
Xv = x(2);
Xt = x(3);
Metab = x(4:7);
mAb_ext = x(8);
mRNA = x(9:10);
Assy = x(11:15);
mAb_gol = x(16);

%% ODEs
[dVdt,dXvdt,dXtdt,mu] = growthStateFcn(V,Xv,Xt,Metab,Q_in);
dMetabdt = metabStateFcn(V,Xv,Metab,Q_in,mu);
[dmAb_extdt,dmRNAdt,dAssydt,dmAb_goldt] = mAbStateFcn(V,Xv,mAb_ext,mRNA,Assy,mAb_gol,Q_in,mu);

%% outputs
dxdt = zeros(16,1);
dxdt(1) = dVdt;
dxdt(2) = dXvdt;
dxdt(3) = dXtdt;
dxdt(4:7) = dMetabdt;
dxdt(8) = dmAb_extdt;
dxdt(9:10) = dmRNAdt;
dxdt(11:15) = dAssydt;
dxdt(16) = dmAb_goldt;
end