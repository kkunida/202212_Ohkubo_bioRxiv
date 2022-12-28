function dxdt = deTremblay1992StateFcnError(x,u)
%% Parameters
mu_max = (1.09/24);       % [/h] % P-M mismatch
mu_death_max3 = 0.69/(24^3); % [/h^3]
Kdeath_Lac = (0.01/24)*2;   % [/(mM*h)] % P-M mismatch
Kdeath_Amm = (0.06/24)*2;   % [/(mM*h)] % P-M mismatch
Kdeath_Gln = 0.02;      % [mM]
Y_X_Glc = 1.09E+8;      % [cell/mmol]
Y_X_Gln = 3.8E+8;       % [cell/mmol]
Y_Lac_Glc = 1.8*2;        % [mmol/mmol] % P-M mismatch
Y_Amm_Gln = 0.85*2;       % [mmol/mmol] % P-M mismatch
mnt_Glc = 0.17E-8/24;   % [mmol/(cell*h)]
Km_Glc = 19.0;          % [mM]
K_Glc = 1.0;            % [mM]
K_Gln = 0.3;            % [mM]
K_mu = 0.02/24;         % [/h]
alpha0 = 2.57E-8/24;    % [mg/(cell*h)]
beta = 0.35E-8/24;      % [mg/(cell*h)]

Klysis = 0.03;          % [/h] from [Kontoravdi+,2010,Comput.Chem.Eng.]

%% Input variables and feed concentration
Q_in = u;
Metab_in = [25; 4; 0; 0;];

% % for multiple feed case
% if size(u,1) == 2 && Q_in > 1E-12
%     Q_Glc = u(1);
%     Q_Gln = u(2);
%     Metab_in(1) = 50*(Q_Glc/Q_in);    % Glc feed conc: 50 mM
%     Metab_in(2) = 8*(Q_Gln/Q_in);     % Gln feed conc: 8 mM
% end

%% State variables
V = x(1);
Xv = x(2);
Xt = x(3);
Metab = x(4:7);
mAb = x(8);

% for multiple feed case
% V_Glc = x(9);
% V_Gln = x(10);

Glc = Metab(1);
Gln = Metab(2);
Lac = Metab(3);
Amm = Metab(4);

%% production rate
mu = mu_max*(Glc/(K_Glc + Glc))*(Gln/(K_Gln + Gln));
mu_death = mu_death_max3*(Kdeath_Gln/(Kdeath_Gln + Gln))/((mu_max - Kdeath_Lac*Lac)*(mu_max - Kdeath_Amm*Amm));
q_Glc = -mu/Y_X_Glc - mnt_Glc*(Glc/(Km_Glc + Glc));
q_Gln = -mu/Y_X_Gln;
q_Lac = -Y_Lac_Glc*q_Glc;
q_Amm = -Y_Amm_Gln*q_Gln;
q_Metab = [q_Glc;q_Gln;q_Lac;q_Amm];
q_mAb = alpha0*(mu/(K_mu + mu)) + beta;

%% ODEs
dVdt = Q_in;
dXvdt = (mu - mu_death)*Xv - (Q_in/V)*Xv;
dXtdt = mu*Xv - (Q_in/V)*Xt - Klysis*(Xt - Xv);
dMetabdt = q_Metab.*Xv - (Q_in/V).*(Metab - Metab_in);
dmAbdt = q_mAb*Xv - (Q_in/V).*mAb;
dxdt = [dVdt; dXvdt; dXtdt; dMetabdt; dmAbdt];

% % for multiple feed case
% if size(u,1) == 2
%     dV_Glc = Q_Glc;
%     dV_Gln = Q_Gln;
%     dxdt = [dxdt; dV_Glc; dV_Gln];
% end

% mu_net = mu - mu_death;

% lower limit of states
dxdt(x < 0 & dxdt < 0) = 0;
end