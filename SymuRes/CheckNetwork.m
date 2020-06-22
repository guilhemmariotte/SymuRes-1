%% Check a network before simulation (Reservoir and Route)
%--------------------------------------------------------------------------

clear all
clc

addpath('Utilityfunctions/','FDfunctions/')
addpath('MFDsolver/','Assignment/','UserNetworks/','PostProc/','Route/','Convergence/')


%% Simulation definition
%--------------------------------------------------------------------------

% Choice of a network defined by user
Simulation.Network = 'City_2res';
% Simulation name
Simulation.Name = 'SC11';


%% Simulation initialization
%--------------------------------------------------------------------------
tic;
disp ' '
disp '********************************************'
disp 'SymuRes V1.1 - network and demand definition'

addpath(['UserNetworks/' Simulation.Network '/'])

SimulSettings

Assignment.CurrentPeriodID = 1;
Assignment.CurrentTime = Assignment.Periods(Assignment.CurrentPeriodID);
Assignment.CurIteration = 1;

RouteCalc

AssignCalc

rmpath(['UserNetworks/' Simulation.Network '/'])

toc;


%% Plot network
%--------------------------------------------------------------------------

% Plot reservoirs
figure
ResList = 1:NumRes;
opts.fontname = 'Arial';
opts.fontsize = 16;
opts.linewidth = 2;
plotResBallConfig(Reservoir,ResList,1,0.5,[1.5 1.5],opts)

figure
plotResNetConfig([],Reservoir,ResList,opts)

% Plot macro nodes
figure
opts.plotlegend = 0;
opts.plotnumnodes = 1;
MacroNodesList = 1:NumMacroNodes;
plotMacroNodes([],Reservoir,ResList,0,MacroNode,MacroNodesList,0,[],[],0,opts)

% Plot routes
figure
RoutesList = 1:NumRoutes;
opts.plotlegend = 0;
opts.nodepath = 1;
plotRoutes([],Reservoir,ResList,0,MacroNode,MacroNodesList,0,Route,RoutesList,1,0,opts)

figure
opts.plotlegend = 0;
opts.plotnumnodes = 0;
plotMacroNodes([],Reservoir,ResList,1,MacroNode,MacroNodesList,1,Route,RoutesList,0,opts)

figure
RoutesList = 1:NumRoutes;
opts.legloc = 'East';
plotResRouteDem([],Reservoir,ResList,Route,RoutesList,'demand','trafficolor',opts)


