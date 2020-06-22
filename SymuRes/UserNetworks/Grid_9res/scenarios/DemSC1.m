%% Scenario 1: SCref (no control), direct route definition from route set
%--------------------------------------------------------------------------

% All macro OD pairs
%--------------------------------------------------------------------------
for od = 1:NumODmacro
    i = 1;
    ODmacro(od).Demand(i).Purpose = 'cartrip';
    ODmacro(od).Demand(i).Time = 0; % [s]
    ODmacro(od).Demand(i).Data = 0; % [veh/s]
end

Simulation.Duration = 3.45*3600; % Simulation duration [s]

Assignment.Periods = [0 Simulation.Duration];
Assignment.PredefRoute = 1;
Assignment.Convergence = 0;

load(['UserNetworks/' Simulation.Network '/networkdata/Grid_9res_routeset_SCref.mat']) % Route

% MacroNode(49).Capacity.Time = 0;
% MacroNode(49).Capacity.Data = 0.45;

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


