%% Scenario 4: multi PID control, OD demand with path flow calculation
%--------------------------------------------------------------------------

Simulation.Duration = 3.45*3600; % Simulation duration [s]

Assignment.Periods = [0:(5*60):Simulation.Duration Simulation.Duration];
Assignment.PredefRoute = 0;
Assignment.Convergence = 1;
Assignment.model = 1; % DUE
Assignment.Behavior = 1; % rational
% Assignment.model = 104; % Wardrop with penalty
Assignment.NumShortestPath = 3;


% PID control for entering R5
%--------------------------------------------------------------------------
Temp_BorderRes = [4 2 6 8];
Temp_Kp = [0.002 0.002 0.002 0.002];
Temp_Ki = [0 0 0 0];
Temp_Kd = [0 0 0 0];
for ictrl = 1:length(Temp_BorderRes)
    
    % PID control for the border with r2
    r2 = Temp_BorderRes(ictrl); % reservoir sending flow to the controled area
    Simulation.Control(ictrl).Type = 'PID'; % controler type
    Simulation.Control(ictrl).Active = 1; % boolean, activate or not the controler
    Simulation.Control(ictrl).ResID = 5; % list of reservoirs to protect
    
    % Get the list of border nodes at the perimeter
    Temp_nodelist = [];
    for r = Simulation.Control(ictrl).ResID
        for inode = Reservoir(r).MacroNodesID
            if strcmp(MacroNode(inode).Type,'border') && isequal(MacroNode(inode).ResID,[r2 r])
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
    Simulation.Control(ictrl).Param = [Temp_Kp(ictrl) Temp_Ki(ictrl) Temp_Kd(ictrl)]; % [Kp Ki Kd] proportional, integral and derivative gains
end
