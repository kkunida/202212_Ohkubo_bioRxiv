function [] = plotSingleFeed(Info_true,Info_error,simout,Ts,x0)
% for Kontoravdi2010 model

[t_true,X_true] = simTrueSFModel(Info_true,Ts,x0);
[t_error,X_error] = simTrueSFModel(Info_error,Ts,x0);
t_MPC = simout.yout{1}.Values.Time;
X_MPC = simout.yout{1}.Values.Data;
t_u_MPC = simout.yout{2}.Values.Time;
U_MPC = simout.yout{2}.Values.Data;


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
inputName = {'Feed rate'};

loweryLim = zeros(18,1);
loweryLim(1) = 0.2;

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
inputNum = size(Info_true.MVopt,2);
for inputIdx = 1:inputNum
    nexttile
    Nstep = size(Info_true.MVopt,1) - 1;
    stairs((0:Nstep)*Ts,Info_true.MVopt(:,inputIdx))
    hold on
    Nstep = size(Info_error.MVopt,1) - 1;
    stairs((0:Nstep)*Ts,Info_error.MVopt(:,inputIdx))
    hold on
    stairs(t_u_MPC,U_MPC(:,inputIdx))
    hold off
    ax = gca;
    ax.PlotBoxAspectRatio = [1.2 0.8 1];
    xlabel('Time [h]')
    xlim([0 168])
    ylim([0 20E-4])
    xticks(24*(1:7))
    xlabel('Time [h]')
    ylabel('[L/h]')
    title(inputName{inputIdx})
    if inputIdx == 1
        legend({'Off-line, full model, no PMM' 'Off-line, reduced model, PMM' 'On-line, reduced model, PMM'},'Location','northeast');
    end
end

% state variables
pos = 3;
for stateIdx = [1 2 4 5 6 7 8]
    nexttile(pos)
    pos = pos + 1;

    plot(t_true,X_true(:,stateIdx), ...
            t_error,X_error(:,stateIdx), ...
            t_MPC, X_MPC(:,stateIdx) ...
        )
    ax = gca;
    ax.PlotBoxAspectRatio = [1.2 0.8 1];
    title(stateName{stateIdx})
    xlabel('Time [h]')
    ylabel(['[' stateUnit{stateIdx} ']'])
    ylim([loweryLim(stateIdx) inf])
    xlim([0 168])
    xticks(24*(1:7))
end

end