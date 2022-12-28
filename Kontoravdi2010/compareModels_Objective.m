function vals = compareModels_Objective(p,Simulator,Exp)
%%
% Define a signal tracking requirement to compute how well the model output
% matches the experiment data. Configure the tracking requirement so that
% it returns the tracking error residuals (rather than the
% sum-squared-error) and does not normalize the errors.
%
r = sdo.requirements.SignalTracking;
r.Type      = '==';
r.Method    = 'Residuals';
r.Normalize = 'off';
% the maximum absolute value of the measured signal is used for normalization

%%
% Update the experiments with the estimated parameter values.
%
Exp  = setEstimatedValues(Exp,p);

%%
% Simulate the model and compare model outputs with measured experiment
% data.
%
Error = [];
for expIdx = 1:numel(Exp)
    
    Simulator = createSimulator(Exp(expIdx),Simulator);
    Simulator = sim(Simulator);

    SimLog  = find(Simulator.LoggedData,get_param('compareModels','SignalLoggingName'));
    mAbDifferenceLog = find(SimLog,'mAbDifference');
   
    mAbDifferenceError = evalRequirement(r,mAbDifferenceLog.Values,Exp(expIdx).OutputData.Values);
    
    Error = [Error; mAbDifferenceError(:)];
end

%%
% Return the residual errors to the optimization solver.
%
vals.F = Error;
end