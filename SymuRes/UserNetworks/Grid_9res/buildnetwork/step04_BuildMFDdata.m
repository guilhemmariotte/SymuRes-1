%% BUILD THE MFD ESTIMATION
%--------------------------------------------------------------------------
% Run this script to estimate the MFD of reservoirs
%
% Use a csv file from Oceane Mascart
%
% June 2020 - Guilhem Mariotte

clear all
clc

addpath('Utilityfunctions/','FDfunctions/')

% Network studied and data files
NetworkName = 'Grid_9res';
FileName = 'Grid_9res';
DataName = 'Donnees_par_zone';

NumRes = 9;


%% Import data
%--------------------------------------------------------------------------

MFDData = cell(1,NumRes);
AggregPeriod = 3*60; % aggregation period of the data [s]

txtfile = fopen(['UserNetworks/' NetworkName '/networkdata/' DataName '.csv']);
data0 = textscan(txtfile,'%s %s %s %s %s',1,'Delimiter',';');
data = textscan(txtfile,'%f %s %f %f %f','Delimiter',';');
fclose(txtfile);
Ndata = length(data{1})/(NumRes+1);

for r = 1:NumRes
    MFDData{r}.t0 = data{1}((r+1):(NumRes+1):end).*AggregPeriod;
    MFDData{r}.V0 = data{3}((r+1):(NumRes+1):end);
    MFDData{r}.P = data{4}((r+1):(NumRes+1):end);
    MFDData{r}.n = data{5}((r+1):(NumRes+1):end);
    
    for i = 2:length(MFDData{r}.t0)
        % Check time steps
        timediff = MFDData{r}.t0(i) - MFDData{r}.t0(i-1);
        if AggregPeriod < timediff
            disp ' '
            disp(datestr(unix2datetime(MFDData{r}.t0(i-1))))
            disp(datestr(unix2datetime(MFDData{r}.t0(i))))
        end
        if AggregPeriod > timediff
            disp ' '
            disp(datestr(unix2datetime(MFDData{r}.t0(i-1))))
            disp(datestr(unix2datetime(MFDData{r}.t0(i))))
        end
        % Check and correct Inf or NaN values
        if abs(MFDData{r}.n(i)) == Inf || isnan(MFDData{r}.n(i))
            if i < length(MFDData{r}.t0)
                MFDData{r}.n(i) = mean([MFDData{r}.n(i-1) MFDData{r}.n(i+1)]);
            else
                MFDData{r}.n(i) = MFDData{r}.n(i-1);
            end
        end
    end
    
    MFDData{r}.V = MFDData{r}.P./MFDData{r}.n; % mean speed [m/s]
    MFDData{r}.t = (0:AggregPeriod:((length(MFDData{r}.n)-1)*AggregPeriod))'; % local time [s]
end

% Check time
for r = 1:length(MFDData)
    disp '---------------------'
    disp(int2str(r))
    for i = 2:length(MFDData{r}.t0)
        timediff = MFDData{r}.t0(i) - MFDData{r}.t0(i-1);
        if AggregPeriod < timediff && timediff < 3*AggregPeriod
            disp ' '
            disp(datestr(unix2datetime(MFDData{r}.t0(i-1))))
            disp(datestr(unix2datetime(MFDData{r}.t0(i))))
        end
        if 5*AggregPeriod < timediff
            disp ' '
            disp(datestr(unix2datetime(MFDData{r}.t0(i-1))))
            disp(datestr(unix2datetime(MFDData{r}.t0(i))))
        end
    end
end


%% MFD curve fitting
%--------------------------------------------------------------------------
% You may run this section many times until you find a suitable fit

MFDpoints = struct('Pc',cell(1,NumRes));

istart = 1;
iend = Ndata/2;

% MFD function model
MFDfct = @paraboFD;
simfct = @(n,param) MFDest(n,param,MFDfct);

% Fit the data to the model using a Particle Swarm Optimization algorithm
for r = 1:NumRes
    ndata = MFDData{r}.n(istart:iend)';
    Pdata = MFDData{r}.P(istart:iend)';
    %nj = 2*max(ndata);
    nj = 2700;
    nc = nj/2;
    %Pc = max(Pdata);
    Pc = 4000;
    paramini = [nj nc Pc];
    %paramdim = [0.1*paramini' 2*paramini'];
    %paramdim = [0.1*paramini' 5*paramini'];
    paramdim = [2500 3000; 0.1*nc 5*nc; 3500 4500];
    cvgcrit = 0.01; % convergence criterium, min RMSE
    maxiter = 50; % max nb of iterations
    Npart = 200; % nb of particles
    param = ParticleSwarmOpt(ndata,Pdata,simfct,paramini,paramdim,cvgcrit,maxiter,Npart);
    
    MFDpoints(r).Pc = param(3);
    MFDpoints(r).nc = param(2);
    MFDpoints(r).njam = param(1);
    MFDpoints(r).u = 2*param(3)/param(2);
    %MFDpoints(r).njam_theoretical = Reservoir(r).MaxAcc;
end

% Plot the data and fitted curve
ResList = 1:length(MFDData);

Nfig = length(ResList); % number of subfigures
Ncol = 3; % number of columns in the figure
Nrow = (floor(Nfig/Ncol) < Nfig/Ncol).*(floor(Nfig/Ncol) + 1) + (floor(Nfig/Ncol) == Nfig/Ncol).*floor(Nfig/Ncol);

% Production-MFD
figure
color0 = 0.5*[1 1 1];
color1 = 0*[1 1 1];
LW = 2;
MS = 4;
for ifig = 1:Nfig
    subplot(Nrow,Ncol,ifig)
    hold on
    n = MFDData{ResList(ifig)}.n(istart:iend)';
    P = MFDData{ResList(ifig)}.P(istart:iend)';
    plot(n,P,'o','color',color0,'markerfacecolor',color0,'markersize',MS)
    nj = MFDpoints(ResList(ifig)).njam;
    nc = MFDpoints(ResList(ifig)).nc;
    Pc = MFDpoints(ResList(ifig)).Pc;
    n_ = linspace(0,nj,200);
    plot(n_,MFDfct(n_,[nj nc Pc]),'-','color',color1,'linewidth',LW)
    hold off
    grid on
    title(['\itR_{\rm' int2str(ResList(ifig)) '}'])
    if ifig + Ncol > Nfig
        xlabel('accumulation \rm[veh]')
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        ylabel('production \rm[veh.m/s]')
    end
end
set(gcf,'Position',[10 10 1000 700])

% Speed-MFD
figure
for ifig = 1:Nfig
    subplot(Nrow,Ncol,ifig)
    hold on
    hold on
    n = MFDData{ResList(ifig)}.n(istart:iend)';
    V = MFDData{ResList(ifig)}.V(istart:iend)';
    plot(n,V,'o','color',color0,'markerfacecolor',color0,'markersize',MS)
    nj = MFDpoints(ResList(ifig)).njam;
    nc = MFDpoints(ResList(ifig)).nc;
    Pc = MFDpoints(ResList(ifig)).Pc;
    n_ = linspace(0,nj,200);
    plot(n_,MFDfct(n_,[nj nc Pc])./n_,'-','color',color1,'linewidth',LW)
    hold off
    grid on
    title(['\itR_{\rm' int2str(ResList(ifig)) '}'])
    if ifig + Ncol > Nfig
        xlabel('accumulation \rm[veh]')
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        ylabel('speed \rm[m/s]')
    end
end
set(gcf,'Position',[10 10 1000 700])


%% Build the empirical data set organized in days
%--------------------------------------------------------------------------

NumDays = 1;

Day = struct('Reservoir',cell(1,NumDays));

for d = 1:NumDays
    %Day(d).Reservoir = Reservoir;
    Day(d).AggregPeriod = AggregPeriod;
    Day(d).IsWeekend = 0;
    
    Nt = Ndata;
    istart = 1;
    iend = Ndata;
    
    [Temp_daynum, Temp_day] = weekday(unix2datetime(MFDData{1}.t0(istart)));
    if strcmp(Temp_day,'Sat') || strcmp(Temp_day,'Sun')
        Day(d).IsWeekend = 1;
    end
    for r = 1:NumRes
        Day(d).Reservoir(r).MaxAcc = MFDpoints(r).njam; % [veh]
        Day(d).Reservoir(r).MaxProd = MFDpoints(r).Pc; % [veh.m/s]
        Day(d).Reservoir(r).CritAcc = MFDpoints(r).nc; % [veh]
        Day(d).Reservoir(r).FreeflowSpeed = MFDpoints(r).u; % [m/s]
        
        Day(d).Reservoir(r).RealTime = MFDData{r}.t0(istart:iend)';
        Day(d).Reservoir(r).LocalTime = MFDData{r}.t(1:Nt)';
        Day(d).Reservoir(r).Acc = MFDData{r}.n(istart:iend)';
        Day(d).Reservoir(r).MeanSpeed = MFDData{r}.V(istart:iend)';
        Day(d).Reservoir(r).Prod = MFDData{r}.P(istart:iend)';
        
        Day(d).Reservoir(r).Inflow = 0;
        Day(d).Reservoir(r).Outflow = 0;
        Day(d).Reservoir(r).AvgTripLength = 0;
    end
end


%% Plot the data
%--------------------------------------------------------------------------

NumTimes = length(Day(1).Reservoir(1).RealTime);
DaysList = 1:length(Day);
Nplot = length(DaysList);

fontname = 'Arial';
LW = 2; % line width
FS = 16; % font size
FS1 = 18; % title font size
MS = 10; % marker size

figure
hold on
for iplot = 1:Nplot
    if Day(DaysList(iplot)).IsWeekend == 1
        color1 = 0.5*[1 1 1];
    else
        color1 = [0 0 0];
    end
    Temp_t = unix2datetime(Day(DaysList(iplot)).Reservoir(1).RealTime);
    Temp_n = zeros(1,NumTimes);
    for r = 1:NumRes
        Temp_n = Temp_n + Day(DaysList(iplot)).Reservoir(r).Acc;
    end
    plot(Temp_t,Temp_n,'-','color',color1,'linewidth',LW)
end
hold off
datetick('x',1,'keepticks')
xlabel('time','FontName',fontname,'FontSize',FS)
ylabel('accumulation \itn \rm[veh]','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1200 300])

% set(gcf,'PaperSize',[15 4],'PaperPosition',[0 0 15 4])
% print('-dpdf',['img/LyonAccChronicle.pdf'])
% beep


%% Save the MFD fit and data points
%--------------------------------------------------------------------------

DataName = 'SCref2';

save(['UserNetworks/' NetworkName '/networkdata/' FileName '_traffic_' DataName '.mat'],'Day','MFDData')
save(['UserNetworks/' NetworkName '/networkdata/' FileName '_MFDfit_' DataName '.mat'],'MFDpoints')





