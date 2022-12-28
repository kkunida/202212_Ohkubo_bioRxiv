function [t_ode15s,X_ode15s] = simTrueMFModel(Info,Ts,x0)

Nstep = size(Info.Xopt,1) - 1;

t_ode15s = 0;
X_ode15s = x0'; 
t0 = 0;
for i = 1:Nstep
    % The manipulated variables are constant only during a single step.
    % ODEFUN needs to be recreated everytime a new prediction step starts 
    % because MV values change.
    u = Info.MVopt(i,:)';
    ODEFUN = @(t,x) multiFeedStateFcn(x,u);

    TSPAN = [t0, t0 + Ts]; % Ts = 0.01
    X_stepend = X_ode15s(end,:)';
    [TOUT,XOUT] = ode15s(ODEFUN,TSPAN,X_stepend);
    t_ode15s = [t_ode15s; TOUT(2:end)];
    X_ode15s = [X_ode15s; XOUT(2:end,:)];
    t0 = t0 + Ts;
end

mAbFinal = X_ode15s(end,8)*X_ode15s(end,1);

fprintf( ...
            [ ...
                'Final liquid volume: %g L\n' ...
                'Final mAb conc: %g mg/L\n' ...
                'Final mAb amount: %g mg\n' ...
            ], ...
            X_ode15s(end,1), ...
            X_ode15s(end,8), ...
            mAbFinal ...
        );

%%
% Plot the states predicted by the MPC with the reduced model and Euler method 
% Plot actual states calculated with the original model and ode15s.

stateName = {'V','Xv','Xt','Glc','Gln','Lac','Amm','mAb_ext', ...
                'mRNA_Hc', 'mRNA_Lc', 'Hc', 'Lc', 'Hc2', 'Hc2Lc', ...
                'mAb_ER', 'mAb_gol', 'V_Glc', 'V_Gln' ...
            };
stateUnit = {'L', 'cell/L', 'cell/L', 'mM', 'mM', 'mM', 'mM', 'mg/L', ...
                'molecule', 'molecule', 'molecule', 'molecule', ...
                'molecule', 'molecule', 'molecule', 'molecule', 'L', 'L' ...
            };

set(groot,'defaulttextinterpreter','none');
set(groot,'defaultlegendinterpreter','none');
figure
for stateIdx = 1:8
    subplot(5,2,stateIdx)
    plot(t_ode15s,X_ode15s(:,stateIdx), ...
        (0:Nstep)*Ts, Info.Xopt(:,stateIdx),'o' ...
    )

    title(stateName{stateIdx})
    ylabel(['[' stateUnit{stateIdx} ']'])
    % ylim([0 inf])
    ylim padded
    xlim([0 168])
    xticks(24*(1:7))
    % legend({'ode15s','NLMPC (ode15s)'},'location','northwest')
end

subplot(5,2,9)
plot(t_ode15s,X_ode15s(:,1).*X_ode15s(:,8), ...
        (0:Nstep)*Ts, Info.Xopt(:,1).*Info.Xopt(:,8),'o' ...
    )
% legend({'ode15s','NLMPC'},'location','northwest')
ylabel('Mass [mg]')
xlim([0 168])
xticks(24*(1:7))
title('mAb mass (x1*x8)')

subplot(5,2,10)
stairs((0:Nstep)*Ts,Info.MVopt(:,1))
hold on
stairs((0:Nstep)*Ts,Info.MVopt(:,2))
hold off
xlabel('Time [h]')
xlim([0 168])
xticks(24*(1:7))
ylabel('[L/h]')
title('Feed rate')
legend({'Q_Glc','Q_Gln'})
sgtitle('Optimized input and states')

end