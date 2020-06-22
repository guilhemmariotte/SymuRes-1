%% Scenario 2: SCref (no control), OD demand with path flow calculation
%--------------------------------------------------------------------------

Simulation.Duration = 3.45*3600; % Simulation duration [s]

% Assignment.Periods = [0 Simulation.Duration];
Assignment.Periods = [0:(20*60):Simulation.Duration Simulation.Duration];
Assignment.PredefRoute = 0;
Assignment.Convergence = 1;
Assignment.model = 1; % DUE
Assignment.Behavior = 1; % rational
% Assignment.model = 104; % Wardrop with penalty
Assignment.NumShortestPath = 3;



% Reservoir weights for assignment model 103
% (we want to favor paths that include R1 and/or R2)
Assignment.ReservoirWeight = zeros(1,NumRes);
Assignment.ReservoirWeight(1) = 2;
Assignment.ReservoirWeight(2) = 2;

% Penalized reservoir paths and penalized coeff for assignment model 104
% (we want to ban paths with transfers between R5 and R1, and between R4
% and R1)
Assignment.PenalizedResPath = {[5 1],[5 1 3],[5 1 4],[1 5],[1 5 2],[1 4],[4 1]};
% Assignment.PenalizedResPath = {};
Assignment.PenalizedAssignCoeff = 0.05;


% PID control
%--------------------------------------------------------------------------
ictrl = 1;
Simulation.Control(ictrl).Type = 'PID'; % controler type
Simulation.Control(ictrl).Active = 0; % boolean, activate or not the controler
Simulation.Control(ictrl).ResID = 1; % list of reservoirs to protect
Simulation.Control(ictrl).MacroNodesID = 1; % list of macro nodes to apply inflow regulation
Simulation.Control(ictrl).SetPointAcc = 0; % control set point: accumulation to track [veh]
Simulation.Control(ictrl).PrevError = 0; % previous step error, for the derivative part [veh]
Simulation.Control(ictrl).IntError = 0; % integrated error, for the integral part [veh]
Simulation.Control(ictrl).MinInflow = 0.1; % min allowed inflow [veh/s]
Simulation.Control(ictrl).MaxInflow = 10; % max allowed inflow [veh/s]
Simulation.Control(ictrl).Param = [0.3 0 0]; % [Kp Ki Kd] proportional, integral and derivative gains


