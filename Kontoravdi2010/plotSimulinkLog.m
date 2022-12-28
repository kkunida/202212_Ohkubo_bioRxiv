function mAb_massFinal = plotSimulinkLog(simout)
stateName = {'V','Xv','Xt','Glc','Gln','Lac','Amm','mAb_ext', ...
                'mRNA_Hc', 'mRNA_Lc', 'Hc', 'Lc', 'Hc2', 'Hc2Lc', ...
                'mAb_ER', 'mAb_gol' ...
            };
stateUnit = {'L', 'cell/L', 'cell/L', 'mM', 'mM', 'mM', 'mM', 'mg/L', ...
                'molecule', 'molecule', 'molecule', 'molecule', ...
                'molecule', 'molecule', 'molecule', 'molecule' ...
            };

T = simout.yout{1}.Values.Time;
X = simout.yout{1}.Values.Data;
T_u = simout.yout{2}.Values.Time;
U = simout.yout{2}.Values.Data;

mAb_mass = X(:,1).*X(:,8);
mAb_massFinal = mAb_mass(end);

fprintf( ...
            [ ...
                'Final liquid volume: %g L\n' ...
                'Final mAb conc: %g mg/L\n' ...
                'Final mAb amount: %g mg\n' ...
            ], ...
            X(end,1), ...
            X(end,8), ...
            mAb_massFinal ...
        );

set(groot,'defaulttextinterpreter','none');
figure
for stateIdx = 1:8
    subplot(5,2,stateIdx)
    plot(T,X(:,stateIdx))

    title(stateName{stateIdx})
    ylabel(['[' stateUnit{stateIdx} ']'])
    % ylim([0 inf])
    ylim padded
    xlim([0 168])
    xticks(24*(1:7))
end


subplot(5,2,9)
plot(T,mAb_mass)
% legend({'ode15s','NLMPC'},'location','northwest')
ylabel('Mass [mg]')
xlim([0 168])
xticks(24*(1:7))
title('mAb mass (x1*x8)')

subplot(5,2,10)
stairs(T_u,U(:,:))
xlabel('Time [h]')
xlim([0 168])
xticks(24*(1:7))
ylabel('[L/h]')
title('Feed rate (u1)')
sgtitle('Optimized input and states')
end