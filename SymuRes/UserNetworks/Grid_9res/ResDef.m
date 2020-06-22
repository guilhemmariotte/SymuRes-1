%% RESERVOIR DEFINITION
%--------------------------------------------------------------------------

%%% RESERVOIR
% Structure which defines all information related to the reservoirs

% Reservoir.Centroid ; virtual positioning [x y] for plotting purpose
% Reservoir.AdjacentRes ; vector, list of adjacent reservoir IDs
% Reservoir.NetLength ; total network length in the reservoir [m]
% Reservoir.FreeflowSpeed ; free-flow speed in the reservoir [m/s]
% Reservoir.MaxProd ; reservoir maximum (critical) production [veh.m/s]
% Reservoir.MaxAcc ; reservoir maximum (jam) accumulation [veh]
% Reservoir.CritAcc ; reservoir critical accumulation [veh]

% Load network structures
load(['UserNetworks/' Simulation.Network '/networkdata/Grid_9res_reservoirs.mat']) % Reservoir, MacroNode
load(['UserNetworks/' Simulation.Network '/networkdata/Grid_9res_ODmacro_SCref.mat']) % ODmacro
load(['UserNetworks/' Simulation.Network '/networkdata/Grid_9res_MFDfit_SCref2.mat']) % MFDpoints

NumRes = length(Reservoir);
NumMacroNodes = length(MacroNode);
NumODmacro = length(ODmacro);

% MFD function
% param = [nj nc Pc]
MFDfct = @paraboFD;

% Entry supply function
% param = [nj nc Pc a1*nc a2*nc b*Pc], with 0 < a1 < 1 < a2, and 1 < b
Entryfct = @paraboEntryFD;
% Entryfct = @(n,param) (n <= param(4)).*param(6) + ...
%     (param(4) < n).*(n <= param(5)).*(param(6)+(n-param(4))./(param(5)-param(4)).*(MFDfct(param(5),param(1:3))-param(6))) + ...
%     (param(5) < n).*MFDfct(n,param(1:3));

Exitfct = @paraboExitFD;

% Reservoir definition

CFit = ones(1,NumRes);
% CFit(4) = 1.2;
% CFit(8) = 1.3;
% CFit(4) = 1.1;
% CFit(8) = 1.1;

% All reservoirs
%--------------------------------------------------------------------------
for r = 1:NumRes
    Reservoir(r).MaxAcc = MFDpoints(r).njam; % [veh]
    Reservoir(r).MaxProd = CFit(r)*MFDpoints(r).Pc; % [veh.m/s]
    Reservoir(r).CritAcc = MFDpoints(r).nc; % [veh]
    Reservoir(r).FreeflowSpeed = MFDpoints(r).u; % [m/s]
    Reservoir(r).MFDfctParam = [Reservoir(r).MaxAcc Reservoir(r).CritAcc Reservoir(r).MaxProd];
    Reservoir(r).EntryfctParam = [Reservoir(r).MaxAcc Reservoir(r).CritAcc Reservoir(r).MaxProd ...
        1.5*Reservoir(r).CritAcc 2*Reservoir(r).CritAcc 1.4*Reservoir(r).MaxProd];
    Reservoir(r).ExitfctParam = Reservoir(r).MFDfctParam;
end

clear MFDpoints






