%% DEMAND DEFINITION
%--------------------------------------------------------------------------

% All macro OD pairs
%--------------------------------------------------------------------------
% for od = 1:NumODmacro
%     i = 1;
%     ODmacro(od).Demand(i).Purpose = 'cartrip';
%     ODmacro(od).Demand(i).Time = 0; % [s]
%     ODmacro(od).Demand(i).Data = 0; % [veh/s]
% end

addpath(['UserNetworks/' Simulation.Network '/scenarios/'])


%% Scenario 1: SCref (no control), direct route definition from route set
%--------------------------------------------------------------------------

if strcmp(Simulation.Name(1:4),'SC10')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'decrdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC1
end

if strcmp(Simulation.Name(1:4),'SC11')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC1
end

if strcmp(Simulation.Name(1:4),'SC12')
    Simulation.MergeModel = 'endogenous'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC1
end

if strcmp(Simulation.Name(1:4),'SC13')
    Simulation.MergeModel = 'demfifo'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC1
end


%% Scenario 2: SCref (no control), OD demand with path flow calculation
%--------------------------------------------------------------------------

if strcmp(Simulation.Name(1:4),'SC20')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'decrdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC2
end

if strcmp(Simulation.Name(1:4),'SC21')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC2
end


%% Scenario 3: one PID control, direct route definition from route set
%--------------------------------------------------------------------------

if strcmp(Simulation.Name(1:4),'SC30')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'decrdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC3
end

if strcmp(Simulation.Name(1:4),'SC31')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC3
end


%% Scenario 4: multi PID control, direct route definition from route set
%--------------------------------------------------------------------------

if strcmp(Simulation.Name(1:4),'SC40')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'decrdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC4
end

if strcmp(Simulation.Name(1:4),'SC41')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC4
end


%% Scenario 5: multi PID control, OD demand with path flow calculation
%--------------------------------------------------------------------------

if strcmp(Simulation.Name(1:4),'SC50')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'decrdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC5
end

if strcmp(Simulation.Name(1:4),'SC51')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC5
end



rmpath(['UserNetworks/' Simulation.Network '/scenarios/'])

