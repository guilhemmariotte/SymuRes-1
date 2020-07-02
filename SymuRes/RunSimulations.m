%% Run a list of simulations
%--------------------------------------------------------------------------

clear all
clc

addpath('Utilityfunctions/','FDfunctions/')
addpath('MFDsolver/','Assignment/','UserNetworks/','PostProc/','Route/','Convergence/')


%% Simulation definition
%--------------------------------------------------------------------------

% SINGLERES: single reservoir with one trip
% Network = 'SingleRes';
% NameSuffix = '_2';
% SolverList = [1 2];
% ScenarioList = {'SC11','SC12','SC13','SC21','SC22','SC23'};
% % ScenarioList = {'SC13','SC23'};

% SINGLERES: single reservoir with two trips
% Network = 'SingleRes';
% NameSuffix = '_2';
% SolverList = [1 2];
% % ScenarioList = {'SC31','SC32','SC33','SC34','SC35','SC41','SC42','SC43','SC44','SC45',...
% %     'SC51','SC52','SC53','SC54','SC55'};
% % ScenarioList = {'SC31','SC32','SC33','SC41','SC42','SC43','SC51','SC52','SC53'};
% ScenarioList = {'SC51','SC52','SC53'};
% % ScenarioList = {'SC34','SC44','SC54'};
% % ScenarioList = {'SC35','SC45','SC55'};

% TWORES: two reservoirs with two trips
% Network = 'TwoRes';
% NameSuffix = '';
% SolverList = [1 2];
% ScenarioList = {'SC10','SC11','SC12','SC13','SC14','SC15','SC20','SC21','SC22','SC23','SC24','SC25',...
%     'SC30','SC31','SC32','SC33','SC34','SC35','SC40','SC41','SC42','SC43','SC44','SC45'};

% TWORES: two reservoirs with two trips, single exit
% Network = 'TwoRes';
% NameSuffix = '_2';
% SolverList = [1 2];
% % ScenarioList = {'SC52','SC53','SC54'};
% % ScenarioList = {'SC15','SC25','SC35','SC45'};
% ScenarioList = {'SC42'};


% BRAESS: Braess network
% Network = 'Braess';
% NameSuffix = '_2';
% SolverList = [1 2];
% % ScenarioList = {'SC11','SC12','SC13','SC21','SC31','SC41','SC51','SC61'};
% % ScenarioList = {'SC11','SC12','SC13','SC14','SC21'};
% ScenarioList = {'SC11','SC31'};


% BRAESS: Braess network
Network = 'Braess_2modes';
NameSuffix = '';
SolverList = [1 2];
% ScenarioList = {'SC11','SC12','SC13','SC21','SC31','SC41','SC51','SC61'};
% ScenarioList = {'SC11','SC12','SC13','SC14','SC21'};
ScenarioList = {'SC11','SC21','SC31'};


% SINGLERESSYMUVIA: Single reservoir comparison with Symuvia
% Network = 'SingleResSymuvia';
% NameSuffix = '';
% SolverList = [1];
% ScenarioList = {'SC361','SC362','SC363','SC364','SC365','SC366'};
% % ScenarioList = {'SC591','SC592','SC593','SC594','SC595','SC596'};
% % ScenarioList = {'SC991','SC992','SC993','SC994','SC995','SC996'};


% CITY_2RES: City CBD-periph network
% Network = 'City_2res';
% NameSuffix = '';
% SolverList = [2];
% % ScenarioList = {'SC11','SC31','SC21','SC41'};
% ScenarioList = {'SC11','SC31'};

% GRID_9RES: Grid 9 res
% Network = 'Grid_9res';
% NameSuffix = '_1';
% SolverList = [1 2];
% % ScenarioList = {'SC10','SC20','SC30','SC40','SC50'};
% ScenarioList = {'SC11','SC21','SC31','SC41','SC51'};
% % ScenarioList = {'SC41','SC51'};


% GRID_25RES: Grid 9 res
% Network = 'Grid_25res';
% NameSuffix = '';
% SolverList = [1];
% % ScenarioList = {'SC11','SC21','SC31','SC41','SC51'};
% ScenarioList = {'SC31'};



for isolver = 1:length(SolverList)
    for iscenario = 1:length(ScenarioList)
        
        % Choice of a network defined by user
        Simulation.Network = Network;
        
        % Choice of the solver
        % 1: accbased / 2: tripbased
        Simulation.Solver = SolverList(isolver);
        
        % Simulation name
        Simulation.Name = [ScenarioList{iscenario} NameSuffix];
        
        % Run the simulation
        Main_Control2
        %Main_DTA
        
    end
end




