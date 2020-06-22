%% BUILD THE RESERVOIR STRUCTURE
%--------------------------------------------------------------------------
% Run this script to create and save the Reservoir structure
%
% Use the automatic build script from the Grid folder
% Use a csv file from Oceane Mascart
%
% June 2020 - Guilhem Mariotte

clear all
clc

addpath('Utilityfunctions/')

% Network studied and data files
NetworkName = 'Grid_9res';
FileName = 'zones';


%% Build the Reservoir structure
%--------------------------------------------------------------------------

% Grid model
% Number of grid rows (> 0)
Nrows = 3;
% Number of grid columns (> 0)
Ncols = 3;

% Number of reservoirs
NumRes = Nrows*Ncols;

Reservoir = struct('Centroid',cell(1,NumRes));

% Spacing between reservoirs
ResSpac = 2;

% Reservoir geometry
r = 1;
for i = 1:Nrows
    for j = 1:Ncols
        xr = (j-1)*ResSpac;
        yr = 2*ResSpac - (i-1)*ResSpac;
        Reservoir(r).ID = ['Zone_' int2str(i) '_' int2str(j)];
        Reservoir(r).Centroid = [xr yr]; % virtual positioning [x y] for plotting purpose
        Reservoir(r).RowID = i;
        Reservoir(r).ColID = j;
        Reservoir(r).BorderPoints = [xr + ResSpac/2*[-1 1 1 -1 -1]; yr + ResSpac/2*[-1 -1 1 1 -1]];        
        r = r + 1;
    end
end

% Adjacent reservoirs
for r = 1:NumRes
    Reservoir(r).AdjacentRes = []; % adjacent reservoirs
    for r2 = 1:NumRes
        Temp_dist = sqrt((Reservoir(r).Centroid(1) - Reservoir(r2).Centroid(1))^2 + (Reservoir(r).Centroid(2) - Reservoir(r2).Centroid(2))^2);
        if Temp_dist < 1.001*ResSpac && r2 ~= r
            Reservoir(r).AdjacentRes = [Reservoir(r).AdjacentRes r2];
        end
    end
end
if NumRes == 1
    Reservoir(1).AdjacentRes = 1;
end


%% Build the Link structure
%--------------------------------------------------------------------------

txtfile = fopen(['UserNetworks/' NetworkName '/networkdata/' FileName '.csv']);
dataset0 = textscan(txtfile,'%s %s %s %s %s %s %s %s %s',1,'Delimiter',';');
dataset = textscan(txtfile,'%s %s %s %s %s %s %s %s %s','Delimiter',';');
fclose(txtfile);

NumLinks = 0;
for r = 1:NumRes
    Ndata = length(dataset{r});
    NumLinks = NumLinks + Ndata;
end
Link = struct('ID',cell(1,NumLinks));

ilink = 1;
for r = 1:NumRes
    Ndata = length(dataset{r});
    Reservoir(r).LinksID = ilink + (0:(Ndata-1));
    for i = 1:Ndata
        Link(ilink).ID = dataset{r}(i);
        Link(ilink).Length = 294; % [m]
        Link(ilink).NumLanes = 2;
        Link(ilink).ResID = r;
        ilink = ilink + 1;
    end
end

% Track some bugs in the network
for j = 1:NumLinks
    if isempty(Link(j).ResID)
        warning(['Link ' int2str(j) ' is not included in any reservoir'])
    end
end

% Reservoir total length
for r = 1:NumRes
    Temp_Ltot = 0;
    for j = Reservoir(r).LinksID
        Temp_Ltot = Temp_Ltot + double(Link(j).NumLanes)*Link(j).Length;
    end
    Reservoir(r).NetLength = Temp_Ltot;
end

% Estimate reservoir max accumulation (based on standard vehicle length)
Temp_kj = 0.17; % lane jam density [veh/m]
for r = 1:NumRes
    Reservoir(r).MaxAcc = floor(Reservoir(r).NetLength*Temp_kj);
end

clear Temp_* i j k l r


%% Build the MacroNode structure
%--------------------------------------------------------------------------

% Origin/entry and destination/exit nodes
i = 1;
for r = 1:NumRes
    MacroNode(i).Type = 'origin';
    MacroNode(i).ResID = r;
    MacroNode(i).Coord = Reservoir(r).Centroid + [0 0.4*ResSpac/2];
    MacroNode(i).Capacity.Time = 0; % [s]
    MacroNode(i).Capacity.Data = 100; % [veh/s]
    i = i + 1;
    
    MacroNode(i).Type = 'destination';
    MacroNode(i).ResID = r;
    MacroNode(i).Coord = Reservoir(r).Centroid + [0 -0.4*ResSpac/2];
    MacroNode(i).Capacity.Time = 0; % [s]
    MacroNode(i).Capacity.Data = 100; % [veh/s]
    i = i + 1;
    
    if Reservoir(r).RowID == 1 || Reservoir(r).RowID == Nrows
        if Reservoir(r).RowID == 1
            Temp_dir = 1;
        elseif Reservoir(r).RowID == Nrows
            Temp_dir = -1;
        end
        MacroNode(i).Type = 'externalentry';
        MacroNode(i).ResID = r;
        MacroNode(i).Coord = Reservoir(r).Centroid + [-0.2*ResSpac/2 Temp_dir*ResSpac/2];
        MacroNode(i).Capacity.Time = 0; % [s]
        MacroNode(i).Capacity.Data = 100; % [veh/s]
        i = i + 1;
        MacroNode(i).Type = 'externalexit';
        MacroNode(i).ResID = r;
        MacroNode(i).Coord = Reservoir(r).Centroid + [0.2*ResSpac/2 Temp_dir*ResSpac/2];
        MacroNode(i).Capacity.Time = 0; % [s]
        MacroNode(i).Capacity.Data = 100; % [veh/s]
        i = i + 1;
    end
    
    if Reservoir(r).ColID == 1 || Reservoir(r).ColID == Ncols
        if Reservoir(r).ColID == 1
            Temp_dir = -1;
        elseif Reservoir(r).ColID == Ncols
            Temp_dir = 1;
        end
        MacroNode(i).Type = 'externalentry';
        MacroNode(i).ResID = r;
        MacroNode(i).Coord = Reservoir(r).Centroid + [Temp_dir*ResSpac/2 0.2*ResSpac/2];
        MacroNode(i).Capacity.Time = 0; % [s]
        MacroNode(i).Capacity.Data = 100; % [veh/s]
        i = i + 1;
        MacroNode(i).Type = 'externalexit';
        MacroNode(i).ResID = r;
        MacroNode(i).Coord = Reservoir(r).Centroid + [Temp_dir*ResSpac/2 -0.2*ResSpac/2];
        MacroNode(i).Capacity.Time = 0; % [s]
        MacroNode(i).Capacity.Data = 100; % [veh/s]
        i = i + 1;
    end
end

% Border nodes
for r = 1:(NumRes-1)
    for r2 = (r+1):NumRes
        if ~isempty(find(Reservoir(r).AdjacentRes == r2,1))
            xr = Reservoir(r).Centroid(1);
            yr = Reservoir(r).Centroid(2);
            xr2 = Reservoir(r2).Centroid(1);
            yr2 = Reservoir(r2).Centroid(2);
            
            % From r to r2
            MacroNode(i).Type = 'border';
            MacroNode(i).ResID = [r r2];
            MacroNode(i).Coord = [(xr + xr2)/2 + 0.2*(yr - yr2) (yr + yr2)/2 + 0.2*(xr2 - xr)];
            MacroNode(i).Capacity.Time = 0; % [s]
            MacroNode(i).Capacity.Data = 100; % [veh/s]
            i = i + 1;
            
            % From r2 to r
            MacroNode(i).Type = 'border';
            MacroNode(i).ResID = [r2 r];
            MacroNode(i).Coord = [(xr + xr2)/2 - 0.2*(yr - yr2) (yr + yr2)/2 - 0.2*(xr2 - xr)];
            MacroNode(i).Capacity.Time = 0; % [s]
            MacroNode(i).Capacity.Data = 100; % [veh/s]
            i = i + 1;
        end
    end
end

NumMacroNodes = i - 1;

% Append to Reservoir structure
for r = 1:NumRes
    Reservoir(r).MacroNodesID = [];
    for i = 1:NumMacroNodes
        if ismember(r,MacroNode(i).ResID)
            Reservoir(r).MacroNodesID = [Reservoir(r).MacroNodesID i];
        end
    end
    Reservoir(r).MacroNodesID = unique(Reservoir(r).MacroNodesID);
end



%% Plot the reservoir network
%--------------------------------------------------------------------------
% Check the reservoir definition

% Plot reservoirs
figure
ResList = 1:NumRes;
opts.fontname = 'Arial';
opts.fontsize = 16;
opts.linewidth = 2;
plotResBallConfig(Reservoir,ResList,1,0.5,[1.5 1.5],opts)

% Plot macro nodes
figure
opts.plotlegend = 0;
opts.plotnumnodes = 1;
MacroNodesList = 1:NumMacroNodes;
plotMacroNodes([],Reservoir,ResList,0,MacroNode,MacroNodesList,0,[],[],0,opts)



%% Save the Reservoir and network structures
%--------------------------------------------------------------------------

save(['UserNetworks/' NetworkName '/networkdata/' NetworkName '_reservoirs.mat'],'Reservoir','MacroNode')
save(['UserNetworks/' NetworkName '/networkdata/' NetworkName '_network.mat'],'Link')


           