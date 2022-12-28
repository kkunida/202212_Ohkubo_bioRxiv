function [] = plotConstFeed(simout)
% for deTremblay1992 model
Tf = 240; % [h] process duration

T_U = simout.yout{1}.Values.Time;
U = simout.yout{1}.Values.Data;

T_true = simout.yout{2}.Values.Time;
X_true = simout.yout{2}.Values.Data;

T_true_error = simout.yout{3}.Values.Time;
X_true_error = simout.yout{3}.Values.Data;

%%
% Plot the states predicted by the MPC with the reduced model and Euler method 
% Plot actual states calculated with the original model and ode15s.

stateName = {'Culture volume','Viable cell density','Total cell density','Glucose','Glutamine','Lactate','Ammonia','mAb', ...
                'mRNA_Hc', 'mRNA_Lc', 'Hc', 'Lc', 'Hc2', 'Hc2Lc', ...
                'mAb_ER', 'mAb_gol', 'V_Glc', 'V_Gln' ...
            };
stateUnit = {'L', 'cell/L', 'cell/L', 'mM', 'mM', 'mM', 'mM', 'mg/L', ...
                'molecule', 'molecule', 'molecule', 'molecule', ...
                'molecule', 'molecule', 'molecule', 'molecule', 'L', 'L' ...
            };
% inputName = {'Feed rate'};

loweryLim = zeros(18,1);
loweryLim(1) = 0.8;

set(groot,'defaulttextinterpreter','none');
set(groot,'defaultlegendinterpreter','none');
set(groot,'DefaultAxesFontName','Arial');
set(groot,'DefaultAxesLineWidth',0.5);
set(groot,'DefaultLineLineWidth',1);
set(groot,'DefaultStairLineWidth',1);

h = figure;
h.Position = [0 0 900 600];

tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

% feed rate
nexttile
plot(T_U,U)
yline(U(1,1))
ax = gca;
ax.PlotBoxAspectRatio = [1.2 0.8 1];
xlabel('Time [h]')
xlim([0 Tf])
ylim([0 1E-2])
xticks(48*(1:Tf/48))
xlabel('Time [h]')
ylabel('[L/h]')
title('Feed rate')

% state variables
pos = 3;
for stateIdx = [1 2 4 5 6 7 8]
    nexttile(pos)
    pos = pos + 1;

    plot(T_true,X_true(:,stateIdx), ...
            T_true_error,X_true_error(:,stateIdx),'--' ...
        )
    ax = gca;
    ax.PlotBoxAspectRatio = [1.2 0.8 1];
    title(stateName{stateIdx})
    xlabel('Time [h]')
    ylabel(['[' stateUnit{stateIdx} ']'])
    ylim([loweryLim(stateIdx) inf])
    xlim([0 Tf])
    xticks(48*(1:Tf/48))
end
legend({'Plant' 'Model with PMM'},'Location','northwest');

end