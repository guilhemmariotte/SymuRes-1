%% Scenario 3: one PID control, direct route definition from route set
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

Assignment.Periods = [0:(2*60):Simulation.Duration Simulation.Duration];
Assignment.PredefRoute = 1;
Assignment.Convergence = 0;

load(['UserNetworks/' Simulation.Network '/networkdata/Grid_9res_routeset_SCref.mat']) % Route


% PID control for entering R5
%--------------------------------------------------------------------------
ictrl = 1;
Simulation.Control(ictrl).Type = 'PID'; % controler type
Simulation.Control(ictrl).Active = 1; % boolean, activate or not the controler
Simulation.Control(ictrl).ResID = 5; % list of reservoirs to protect

% Get the list of border nodes at the perimeter
Temp_nodelist = [];
for r = Simulation.Control(ictrl).ResID
    for inode = Reservoir(r).MacroNodesID
        if strcmp(MacroNode(inode).Type,'border') && MacroNode(inode).ResID(2) == r && ...
                ~ismember(MacroNode(inode).ResID(1),Simulation.Control(ictrl).ResID)
            Temp_nodelist = [Temp_nodelist inode];
        end
    end
end
Simulation.Control(ictrl).MacroNodesID = Temp_nodelist; % list of macro nodes to apply inflow regulation

% Set the control set point: sum of critical accumulations
Temp_n0 = 0;
for r = Simulation.Control(ictrl).ResID
    Temp_n0 = Temp_n0 + Reservoir(r).CritAcc;
end
Simulation.Control(ictrl).SetPointAcc = Temp_n0; % control set point: accumulation to track [veh]
Simulation.Control(ictrl).PrevError = 0; % previous step error, for the derivative part [veh]
Simulation.Control(ictrl).IntError = 0; % integrated error, for the integral part [veh]

% Inflow bounds
Simulation.Control(ictrl).MinInflow = 0.1; % min allowed inflow [veh/s]
Simulation.Control(ictrl).MaxInflow = 10; % max allowed inflow [veh/s]

% Controler parameters
Simulation.Control(ictrl).Param = [0.008 0 0]; % [Kp Ki Kd] proportional, integral and derivative gains

