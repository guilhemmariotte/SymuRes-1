%% BUILD THE MACRO TRIP SET
%--------------------------------------------------------------------------
% Run this script to create and save the MacroChoice structure
% This structure contains a set of possible macro routes with their trip
% length information and number of corresponding micro trips on the real
% network
%
% Use the trip set file for Symuvia from Oceane Mascart
%
% June 2020 - Guilhem Mariotte

clear all
clc

% Network studied and data files
NetworkName = 'Grid_9res';
FileName = 'Grid_9res';
DataName = 'tripset';
DataName2 = 'SCref';

% Load network structures
load(['UserNetworks/' NetworkName '/networkdata/' FileName '_reservoirs.mat']) % Reservoir, MacroNode
load(['UserNetworks/' NetworkName '/networkdata/' FileName '_network.mat']) % Link

NumRes = length(Reservoir);
NumMacroNodes = length(MacroNode);
NumLinks = length(Link);


%% Get the vehicle trip set
%--------------------------------------------------------------------------

txtfile = fopen(['UserNetworks/' NetworkName '/networkdata/' DataName '.csv']);
dataset0 = textscan(txtfile,'%s %s %s %s %s',1,'Delimiter',';');
dataset = textscan(txtfile,'%s %s %f %s %s','Delimiter',';');
fclose(txtfile);

Ndata = length(dataset{1});
Vehicle = struct('LinkPath',cell(1,Ndata));

for idata = 1:Ndata
    
    % Get the link path
    Temp_path = dataset{4}(idata);
    Temp_path = strsplit(Temp_path{1},' ');
    Temp_linkpath = []; % sequence of link IDs
    for i = 1:length(Temp_path)
        j = 1;
        while j <= NumLinks && ~strcmp(Link(j).ID,Temp_path{i})
            j = j + 1;
        end
        if j > NumLinks % no ID match found
            warning(['Link with ID ' Temp_path{i} ' in vehicle path ' int2str(idata) ' is not recognized'])
        else
            Temp_linkpath = [Temp_linkpath j];
        end
    end
    Vehicle(idata).LinkPath = Temp_linkpath;
    
    % Get the reservoir path
    if ~isempty(Vehicle(idata).LinkPath)
        j = Vehicle(idata).LinkPath(1);
        Temp_respath = Link(j).ResID; % sequence of reservoir IDs
        Temp_triplengths = zeros(1,10*NumRes); % default reservoir path length
        Temp_triplengths(1) = Link(j).Length;
        if length(Vehicle(idata).LinkPath) > 1
            i_r = 1;
            for i = 2:length(Vehicle(idata).LinkPath)
                j = Vehicle(idata).LinkPath(i);
                j2 = Vehicle(idata).LinkPath(i-1);
                if Link(j).ResID ~= Link(j2).ResID
                    Temp_respath = [Temp_respath Link(j).ResID];
                    i_r = i_r + 1;
                end
                Temp_triplengths(i_r) = Temp_triplengths(i_r) + Link(j).Length;
            end
        end
        Vehicle(idata).ResPath = Temp_respath;
        Vehicle(idata).TripLengths = Temp_triplengths(1:i_r);
    end
    
    % Creation time
    Vehicle(idata).CreationTime = dataset{3}(idata); % [s]
    
    % Number of micro routes corresponding to the macro route
    Vehicle(idata).NumMicroTrips = 1;
    
end


%% Save the vehicle trip set
%--------------------------------------------------------------------------

save(['UserNetworks/' NetworkName '/networkdata/' NetworkName '_tripset_' DataName2 '.mat'],'Vehicle')


%% Get the macro routes from the trip set
%--------------------------------------------------------------------------

VehList = 1:length(Vehicle);
NumVeh = length(VehList); % number of vehicles

% Gather vehicles by similar Node paths
VehSimilarity = zeros(NumVeh,NumVeh);
for i = 1:NumVeh % loop on all vehicles
    for j = 1:NumVeh % loop on all vehicles
        if isequal(Vehicle(VehList(i)).ResPath,Vehicle(VehList(j)).ResPath)
            VehSimilarity(i,j) = 1;
        end
    end
end
Temp_vehindexes = gatherelements(VehSimilarity);

% Build the MacroChoice structure
NumMacroChoice = length(Temp_vehindexes);
MacroChoice = struct('OriginRes',cell(1,NumMacroChoice));
for ivehset = 1:NumMacroChoice
    if ~isempty(Temp_vehindexes{ivehset})
        iveh = VehList(Temp_vehindexes{ivehset}(1));
        MacroChoice(ivehset).OriginRes = Vehicle(iveh).ResPath(1);
        MacroChoice(ivehset).DestinationRes = Vehicle(iveh).ResPath(end);
        MacroChoice(ivehset).MacroPath = Vehicle(iveh).ResPath;
        %MacroChoice(ivehset).NodePath = Vehicle(iveh).NodePath;
        MacroChoice(ivehset).Ntrips = 0;
        
        % Number of trips
        Temp_nveh = length(Temp_vehindexes{ivehset});
        Temp_ntrips = 0;
        Temp_lengths = [];
        Temp_times = [];
        i = 1;
        for iveh = VehList(Temp_vehindexes{ivehset})
            Temp_ntrips = Temp_ntrips + Vehicle(iveh).NumMicroTrips;
            Temp_lengths = [Temp_lengths  Vehicle(iveh).TripLengths'];
            Temp_times = [Temp_times  Vehicle(iveh).CreationTime];
            i = i + 1;
        end
        MacroChoice(ivehset).Ntrips = Temp_ntrips;
        
        % Trip lengths
        if Temp_ntrips > 1
            MacroChoice(ivehset).Mean = mean(Temp_lengths,2)';
            MacroChoice(ivehset).Std = std(Temp_lengths,0,2)';
        else
            MacroChoice(ivehset).Mean = Temp_lengths';
            MacroChoice(ivehset).Std = Temp_lengths';
        end
        MacroChoice(ivehset).Length = Temp_lengths;
        
        % Creation times
        MacroChoice(ivehset).CreationTimes = Temp_times;
    end
end

% Clean the MacroChoice set
MacroChoice0 = MacroChoice;
NumMacroChoice0 = length(MacroChoice0);
NumTotTrips0 = 0; % total number of trips in the dataset
NumTotTrips = 0; % total number of trips after cleaning the MacroChoice set
NbTripMin = 5; % min number of trips to define a macro route
ichoice = 1;
for ichoice0 = 1:NumMacroChoice0
    NumTotTrips0 = NumTotTrips0 + MacroChoice0(ichoice0).Ntrips;
    if MacroChoice0(ichoice0).Ntrips > NbTripMin && ... % keep macro routes represented by a min number of trips
            isequal(MacroChoice0(ichoice0).MacroPath,unique(MacroChoice0(ichoice0).MacroPath,'stable')) % and without reservoir repetition (e.g. 1,2,1,2,...)
        MacroChoice(ichoice) = MacroChoice0(ichoice0);
        NumTotTrips = NumTotTrips + MacroChoice(ichoice).Ntrips;
        ichoice = ichoice + 1;
    end
end
MacroChoice = MacroChoice(1:(ichoice-1));
disp(['Trip coverage: ' num2str(NumTotTrips/NumTotTrips0)])


%% Save the macro trip set
%--------------------------------------------------------------------------

save(['UserNetworks/' NetworkName '/networkdata/' NetworkName '_triplengths_' DataName2 '.mat'],'MacroChoice')



