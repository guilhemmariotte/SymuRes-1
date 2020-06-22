%% BUILD THE ODMACRO STRUCTURE
%--------------------------------------------------------------------------
% Run this script to create and save the ODmacro structure
% Build the regional OD matrix (all macro node origins/external entries to
% all macro node destinations/external exits) with the set of possible
% macro routes and their trip lengths associated to each OD pair.
% A demand profile can be also associated to each OD. Otherwise, the demand
% is to be defined in the simulation scenario (DemDef file).
%
% Based on the MacroChoice file built with the trip set file for Symuvia
% from Oceane Mascart
%
% June 2020 - Guilhem Mariotte

clear all
clc

% Network studied and data files
NetworkName = 'Grid_9res';
FileName = 'Grid_9res';
TripLengthFile = 'SCref';


% Load network structures
load(['UserNetworks/' NetworkName '/networkdata/' FileName '_reservoirs.mat']) % Reservoir, MacroNode
load(['UserNetworks/' NetworkName '/networkdata/' FileName '_network.mat']) % Link
load(['UserNetworks/' NetworkName '/networkdata/' FileName '_triplengths_' TripLengthFile '.mat']) % MacroChoice

NumRes = length(Reservoir);
NumMacroNodes = length(MacroNode);
NumLinks = length(Link);
NumPossibleRoutes = length(MacroChoice);


%% Calculate the flow per route and modify the reservoir path
%--------------------------------------------------------------------------

TimeStep = 600; % time step for the demand profiles definition [s]
SimulationDuration = 3.45*3600; % [s]

TripLengthAdd = 2*294; % [m] default trip length to consider in the additional intermediate res

Nrows = max([Reservoir(1:NumRes).RowID]);
Ncols = max([Reservoir(1:NumRes).ColID]);
GridRes = zeros(Nrows,Ncols); % reservoir ID by grid row and column ID
for r = 1:NumRes
    i = Reservoir(r).RowID;
    j = Reservoir(r).ColID;
    GridRes(i,j) = r;
end

for i = 1:NumPossibleRoutes
    
    % Modification of the reservoir path in case of non adjacent reservoirs
    Temp_respath = MacroChoice(i).MacroPath;
    Temp_newrespath = Temp_respath(1);
    Temp_newmeanlengths = MacroChoice(i).Mean(1);
    Temp_newstdlengths = MacroChoice(i).Std(1);
    if length(Temp_respath) > 1
        for i_r = 2:length(Temp_respath)
            r = Temp_respath(i_r); % current res in the path
            r2 = Temp_respath(i_r-1); % previous res in the path
            if ~ismember(r,Reservoir(r2).AdjacentRes)
                % if two successive res not adjacent
                if Reservoir(r).RowID == Reservoir(r2).RowID+1 && Reservoir(r).ColID == Reservoir(r2).ColID+1
                    r3 = GridRes(Reservoir(r2).RowID,Reservoir(r2).ColID+1);
                elseif Reservoir(r).RowID == Reservoir(r2).RowID-1 && Reservoir(r).ColID == Reservoir(r2).ColID-1
                    r3 = GridRes(Reservoir(r2).RowID,Reservoir(r2).ColID-1);
                elseif Reservoir(r).RowID == Reservoir(r2).RowID+1 && Reservoir(r).ColID == Reservoir(r2).ColID-1
                    r3 = GridRes(Reservoir(r2).RowID+1,Reservoir(r2).ColID);
                elseif Reservoir(r).RowID == Reservoir(r2).RowID-1 && Reservoir(r).ColID == Reservoir(r2).ColID+1
                    r3 = GridRes(Reservoir(r2).RowID-1,Reservoir(r2).ColID);
                else
                    r3 = [];
                    warning(['Did not find an intermediate reservoir for MacroChoice route ' int2str(i)])
                end
                Temp_newrespath = [Temp_newrespath r3]; % add an intermediate reservoir
                Temp_newmeanlengths = [Temp_newmeanlengths TripLengthAdd];
                Temp_newstdlengths = [Temp_newstdlengths 0];
            end
            Temp_newrespath = [Temp_newrespath r];
            Temp_newmeanlengths = [Temp_newmeanlengths MacroChoice(i).Mean(i_r)];
            Temp_newstdlengths = [Temp_newstdlengths MacroChoice(i).Std(i_r)];
        end
    end
    MacroChoice(i).MacroPath0 = Temp_respath; % old reservoir path
    MacroChoice(i).MacroPath = Temp_newrespath; % new reservoir path with (possibly) additional res
    MacroChoice(i).Std0 = MacroChoice(i).Std;
    MacroChoice(i).Std = Temp_newstdlengths; % new reservoir trip length std
    MacroChoice(i).Mean0 = MacroChoice(i).Mean;
    MacroChoice(i).Mean = Temp_newmeanlengths; % new reservoir trip length mean
    
    
    % Flow estimation (demand for the route)
    Temp_times0 = sort(MacroChoice(i).CreationTimes);
    Temp_N0 = 1:length(MacroChoice(i).CreationTimes);
    for i_t = 2:length(Temp_times0)
        if Temp_times0(i_t) == Temp_times0(i_t-1)
            Temp_times0(i_t) = Temp_times0(i_t) + 0.01; % slight shift to avoid numeric troubles
        end
    end
    [Temp_times, Temp_N] = pwfct(Temp_times0,Temp_N0,TimeStep,0,SimulationDuration);
    Temp_dem = deriv(Temp_times,Temp_N); % demand profile as the derivative of the cumul veh count
    if isnan(sum(Temp_dem))
        warning(['Demand for MacroChoice route ' int2str(i) ' has NaNs'])
    end
    Temp_dem = (Temp_dem >= 0).*Temp_dem;
    MacroChoice(i).Demand.Time = Temp_times; % [s]
    MacroChoice(i).Demand.Data = Temp_dem; % [veh/s]
end

DemandTime = MacroChoice(1).Demand.Time;
NumTimes = length(DemandTime);


%% Build the ODmacro structure
%--------------------------------------------------------------------------

% Find origins and destinations (entry or exit through a macro node)
Temp_orires = [];
Temp_orinodes = [];
Temp_destres = [];
Temp_destnodes = [];
for r = 1:NumRes
    for i = Reservoir(r).MacroNodesID
        if strcmp(MacroNode(i).Type,'origin') || strcmp(MacroNode(i).Type,'externalentry') % origins
            Temp_orires = [Temp_orires r];
            Temp_orinodes = [Temp_orinodes i];
        end
        if strcmp(MacroNode(i).Type,'destination') || strcmp(MacroNode(i).Type,'externalexit') % destinations
            Temp_destres = [Temp_destres r];
            Temp_destnodes = [Temp_destnodes i];
        end
    end
end
NumOrigins = length(Temp_orires);
NumDestinations = length(Temp_destres);

NumODmacro = NumOrigins*NumDestinations;
ODmacro = struct('OriginID',cell(1,NumODmacro));

Temp_MacroDemUsed = zeros(1,NumPossibleRoutes); % boolean, 1 if the route demand has been used in ODmacro
Route = struct('ODmacroID',cell(1,NumPossibleRoutes));
iroute2 = 1;

od = 1;
for o = 1:NumOrigins % loop on all possible origins
    for d = 1:NumDestinations % loop on all possible destinations
        ro = Temp_orires(o);
        rd = Temp_destres(d);
        no = Temp_orinodes(o);
        nd = Temp_destnodes(d);
        
        ODmacro(od).OriginID = o;
        ODmacro(od).DestinationID = d;
        ODmacro(od).ResOriginID = ro;
        ODmacro(od).ResDestinationID = rd;
        ODmacro(od).NodeOriginID = no;
        ODmacro(od).NodeDestinationID = nd;
        
        % Possible routes with their lengths
        Temp_data = zeros(1,NumTimes);
        iroute = 1;
        for i = 1:NumPossibleRoutes
            if MacroChoice(i).OriginRes == ro && MacroChoice(i).DestinationRes == rd
                Temp_validroute = 1; % valid route by default
                Temp_respath = MacroChoice(i).MacroPath;
                Temp_nodepath = [];
                if length(Temp_respath) > 1
                    for i_r = 1:(length(Temp_respath)-1)
                        inode = 1;
                        while inode <= NumMacroNodes && ~isequal([Temp_respath(i_r) Temp_respath(i_r+1)],MacroNode(inode).ResID)
                            inode = inode + 1;
                        end
                        if inode > NumMacroNodes
                            warning(['Possible route ' int2str(i) ' cannot get a MacroNode path'])
                            Temp_validroute = 0; % not a valid route
                        else
                            Temp_nodepath = [Temp_nodepath inode];
                        end
                    end
                end
                
                if Temp_validroute == 1
                    
                    % Record the possible route only if valid
                    ODmacro(od).PossibleRoute(iroute).ResPath = Temp_respath;
                    ODmacro(od).PossibleRoute(iroute).NodePath = [no Temp_nodepath nd];
                    ODmacro(od).PossibleRoute(iroute).TripLengths = MacroChoice(i).Mean;
                    ODmacro(od).PossibleRoute(iroute).NumMicroTrips = MacroChoice(i).Ntrips;
                    
                    if Temp_MacroDemUsed(i) == 0 && ...
                            strcmp(MacroNode(no).Type,'externalentry') && strcmp(MacroNode(nd).Type,'externalexit')
                        % Record the demand and the route
                        Temp_data = Temp_data + MacroChoice(i).Demand.Data;
                        Temp_MacroDemUsed(i) = 1;
                        Route(iroute2).ODmacroID = od;
                        Route(iroute2).ResPath = ODmacro(od).PossibleRoute(iroute).ResPath;
                        Route(iroute2).NodePath = ODmacro(od).PossibleRoute(iroute).NodePath;
                        Route(iroute2).TripLengths = ODmacro(od).PossibleRoute(iroute).TripLengths;
                        j = 1;
                        Route(iroute2).Demand0(j).Purpose = 'cartrip';
                        Route(iroute2).Demand0(j).Time = DemandTime; % [s]
                        Route(iroute2).Demand0(j).Data = MacroChoice(i).Demand.Data; % [veh/s]
                        iroute2 = iroute2 + 1;
                    end
                    
                    iroute = iroute + 1;
                end
            end
        end
        ODmacro(od).NumPossibleRoutes = iroute - 1;
        
        % Demand
        j = 1;
        ODmacro(od).Demand(j).Purpose = 'cartrip';
        ODmacro(od).Demand(j).Time = DemandTime; % [s]
        ODmacro(od).Demand(j).Data = Temp_data; % [veh/s]
        
        od = od + 1;
    end
end


NumRoutes = iroute2 - 1;

% Check total demand
Temp_demtot = zeros(1,NumTimes);
for i = 1:NumPossibleRoutes
    Temp_demtot = Temp_demtot + MacroChoice(i).Demand.Data;
end
Temp_demtot2 = zeros(1,NumTimes);
for od = 1:NumODmacro
    Temp_demtot2 = Temp_demtot2 + ODmacro(od).Demand(1).Data;
end
Temp_demtot3 = zeros(1,NumTimes);
for iroute = 1:NumRoutes
    Temp_demtot3 = Temp_demtot3 + Route(iroute).Demand0(1).Data;
end
Temp_demtotcomp = [Temp_demtot; Temp_demtot2; Temp_demtot3];


%% Save the ODmacro structure
%--------------------------------------------------------------------------

save(['UserNetworks/' NetworkName '/networkdata/' FileName '_ODmacro_' TripLengthFile '.mat'],'ODmacro')
save(['UserNetworks/' NetworkName '/networkdata/' FileName '_routeset_' TripLengthFile '.mat'],'Route')

