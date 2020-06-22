%% PLOT RESULTS
%--------------------------------------------------------------------------
% DO NOT LAUNCH THIS SCRIPT AT ONCE! Many plots and videos would appear at
% the same time. Go to each section to see what it does.

clear all
clc

addpath('Utilityfunctions/','FDfunctions/')
addpath('MFDsolver/','Assignment/','UserNetworks/','PostProc/')


% GRAPHIC DESIGN
%---------------
% cmap_perso: blue; green; red; yellow; violet; light blue; light violet; light yellow
cmap_perso = [51 51 255; 0 204 51; 204 0 0; 204 153 0; 153 0 102; 51 153 153; 204 102 204; 204 204 102]/255;
cmap_perso2 = [51 51 255; 0 204 51; 204 0 0; 204 153 0; 153 0 102;    51 153 153; 153 204 102; 153 102 102; 204 204 102; 204 102 204]/255;
% cmap_ifsttar: dark blue; blue; blue-green; green
cmap_ifsttar = [0 83 151; 0 118 189; 0 166 164; 87 171 39]/255;
marker_perso = ['o' '+' '*' '^' 'x' 'd' 'v' 's' '<' '>'];
line_perso = {'-', '--', ':', '-.','-', '--', ':', '-.','-', '--', ':', '-.','-', '--', ':', '-.'};
fontname = 'Times New Roman';


%% Load simulation files and launch post-processing
%--------------------------------------------------------------------------

% Result from simulation 1
%--------------------------------------------------------------------------
iresu = 1;
Results(iresu).Network = 'Grid_9res'; % Choice of a network defined by user
Results(iresu).Solver = 1; % Choice of the solver. 1: accbased / 2: tripbased
Results(iresu).Name = 'SC11'; % Simulation name
Results(iresu).Name2 = 'acc-based'; % Name to print on the graph legends

if Results(iresu).Solver == 1
    load(['UserNetworks/' Results(iresu).Network '/outputs/Outputs_' Results(iresu).Name '_accbased.mat'])
    PostProc_accbased
elseif Results(iresu).Solver == 2
    load(['UserNetworks/' Results(iresu).Network '/outputs/Outputs_' Results(iresu).Name '_tripbased.mat'])
    PostProc_tripbased
end
Results(iresu).Reservoir = Reservoir;
Results(iresu).Route = Route;
Results(iresu).Assignment = Assignment;
Results(iresu).SimulTime = Simulation.Time;
Results(iresu).ODmacro = ODmacro;


% Result from Symuvia
%--------------------------------------------------------------------------
load(['UserNetworks/' Results(1).Network '/networkdata/Grid_9res_traffic_SCref.mat'])
iresu = 2;
for iday = 1:length(Day)
    Results(iresu).Network = 'Grid_9res'; % Choice of a network defined by user
    Results(iresu).Solver = 1; % Choice of the solver. 1: accbased / 2: tripbased
    Results(iresu).Name = 'Symuvia'; % Simulation name
    Results(iresu).Name2 = 'microsim'; % Name to print on the graph legends
    
    Reservoir = Day(iday).Reservoir;
    istart = 1;
    iend = length(Reservoir(1).Acc);
    for r = 1:length(Reservoir)
        Temp_t = Reservoir(r).LocalTime(1:(iend-istart+1));
        Temp_n = Reservoir(r).Acc(istart:iend);
        Temp_v = Reservoir(r).MeanSpeed(istart:iend);
        Reservoir(r).Acc = resamp(Simulation.Time,Temp_t,Temp_n);
        Reservoir(r).MeanSpeed = resamp(Simulation.Time,Temp_t,Temp_v);
    end
    
    Results(iresu).Reservoir = Reservoir;
    Results(iresu).SimulTime = Simulation.Time; % 2011 data
    
    iresu = iresu + 1;
end


% Result from simulation 2
%--------------------------------------------------------------------------
iresu = 3;
Results(iresu).Network = 'Grid_9res'; % Choice of a network defined by user
Results(iresu).Solver = 2; % Choice of the solver. 1: accbased / 2: tripbased
Results(iresu).Name = 'SC11'; % Simulation name
Results(iresu).Name2 = 'trip-based'; % Name to print on the graph legends

if Results(iresu).Solver == 1
    load(['UserNetworks/' Results(iresu).Network '/outputs/Outputs_' Results(iresu).Name '_accbased.mat'])
    PostProc_accbased
elseif Results(iresu).Solver == 2
    load(['UserNetworks/' Results(iresu).Network '/outputs/Outputs_' Results(iresu).Name '_tripbased.mat'])
    PostProc_tripbased
end
Results(iresu).Reservoir = Reservoir;
Results(iresu).Route = Route;
Results(iresu).Assignment = Assignment;
Results(iresu).SimulTime = Simulation.Time;
Results(iresu).ODmacro = ODmacro;



%% Total accumulation and mean speed in reservoirs
%--------------------------------------------------------------------------

FS = 9; % font size
FS1 = 12; % title font size
figindex = 'abcdefghijklmnopqrstuvwxyz';
PrintFig = 0; % 1: print the figure (save as pdf in the img folder) / 0: do not print
filename = [Results(1).Name '_' Results(2).Name]; % name for printing

% Plot options
ResList = 1:NumRes; % list of reservoirs
ResuList = [1 3 2]; % list of results
TimeRange = [0 Simulation.Duration]; % [s]
AccRange = [0 2000]; % [veh]
FlowRange = [0 3]; % [veh/s]
Nfig = length(ResList); % number of subfigures
Nresu = length(ResuList); % number of results

% Color and line options
% sline = arrayextension(line_perso,Nresu,'column'); % line styles for results
sline = arrayextension({'-'},Nresu,'column'); % line styles for results
% LW = 1*ones(1,Nresu); % line widths for results
LW = linspace(3,1,Nresu); % line widths for results
clight = linspace(0,(Nresu-1)/Nresu,Nresu); % color brightness (0 to 1) for results
% cmap = arrayextension(cmap_perso2,Nresu,'row'); % colors for results
cmap = [linspace(0.3,0.7,(Nresu-1)) 0]'*[1 1 1]; % colors for results

% Figure and subfigure options
Ncol = 3; % number of columns in the figure
Nrow = (floor(Nfig/Ncol) < Nfig/Ncol).*(floor(Nfig/Ncol) + 1) + (floor(Nfig/Ncol) == Nfig/Ncol).*floor(Nfig/Ncol);
[colindex, rowindex] = ind2sub([Ncol Nrow],1:Nfig);
marginleft = 0.1; % figure relative left margin
marginbottom = 0.1; % figure relative bottom margin
intermarginright = 0.15; % subplot relative right margin
intermarginbottom = 0.05; % subplot relative bottom margin
intermargintop = 0.2; % subplot relative top margin
figwidth = 6.5; % whole figure width [inches]
figheight = Nrow*2; % whole figure height [inches]

xt = 7:11;
xtl = cell(1,length(xt));
for i = 1:length(xt)
    xtl{i} = [int2str(xt(i)) ':00'];
end

% Accumulation in reservoirs
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nresu,1)));
hp1 = zeros(1,Nresu);
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    hold on
    i = 1;
    for iresu = 1:Nresu
        Temp_t = Results(ResuList(iresu)).SimulTime;
        Temp_n = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Acc;
        hp1(i) = plot(Temp_t,Temp_n,...
            'linestyle',sline{iresu},'color',cmap(iresu,:),'linewidth',LW(iresu));
        strleg{i} = Results(ResuList(iresu)).Name2;
        i = i + 1;
    end
    hold off
    %axis([TimeRange AccRange])
    xlim(TimeRange)
    title(['(' figindex(ifig) ') - \itR_{\rm' int2str(ResList(ifig)) '}'],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    %xticks(xt)
    %xticklabels(xtl)
    if ifig + Ncol > Nfig
        xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        ylabel('accumulation \itn^r \rm[veh]','FontName',fontname,'FontSize',FS)
    end
    if ifig == 1
        hleg = legend(hp1,strleg);
        legpos = get(hleg,'Position');
        %set(hleg,'Position',[0.6 -0.02 legpos(3) legpos(4)],'FontName',fontname,'FontSize',0.8*FS)
        set(hleg,'Position',[0.6 0.001 legpos(3) legpos(4)],'FontName',fontname,'FontSize',FS) % for 2 col in fig
        %set(hleg,'Position',[0.7 0.05 legpos(3) legpos(4)],'FontName',fontname,'FontSize',0.8*FS) % for 3 col in fig
    end
    xj = marginleft + (colindex(ifig) - 1)*(1 - marginleft)/Ncol; % left pos
    yi = marginbottom + (Nrow - rowindex(ifig) + intermarginbottom)*(1 - marginbottom)/Nrow; % bottom pos
    w = (1 - intermarginright)*(1 - marginleft)/Ncol; % width
    h = (1 - intermargintop - intermarginbottom)*(1 - marginbottom)/Nrow; % height
    set(gca,'Position',[xj yi w h],'FontName',fontname,'FontSize',FS)
end
if PrintFig == 1
    set(gcf,'PaperUnits','inches')
    set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
    print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_acc_' filename '.pdf'])
else
    set(gcf,'Position',[10 10 1000 700])
end


% Mean speed in reservoirs
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nresu,1)));
hp1 = zeros(1,Nresu);
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    hold on
    i = 1;
    for iresu = 1:Nresu
        Temp_t = Results(ResuList(iresu)).SimulTime;
        Temp_v = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MeanSpeed;
        hp1(i) = plot(Temp_t,Temp_v,...
            'linestyle',sline{iresu},'color',cmap(iresu,:),'linewidth',LW(iresu));
        strleg{i} = Results(ResuList(iresu)).Name2;
        i = i + 1;
    end
    hold off
    %axis([TimeRange SpeedRange])
    xlim(TimeRange)
    title(['(' figindex(ifig) ') - \itR_{\rm' int2str(ResList(ifig)) '}'],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    %xticks(xt)
    %xticklabels(xtl)
    if ifig + Ncol > Nfig
        xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        ylabel('mean speed \itv^r \rm[m/s]','FontName',fontname,'FontSize',FS)
    end
    if ifig == 1
        hleg = legend(hp1,strleg);
        legpos = get(hleg,'Position');
        %set(hleg,'Position',[0.6 -0.02 legpos(3) legpos(4)],'FontName',fontname,'FontSize',0.8*FS)
        set(hleg,'Position',[0.6 0.001 legpos(3) legpos(4)],'FontName',fontname,'FontSize',FS) % for 2 col in fig
        %set(hleg,'Position',[0.7 0.05 legpos(3) legpos(4)],'FontName',fontname,'FontSize',0.8*FS) % for 3 col in fig
    end
    xj = marginleft + (colindex(ifig) - 1)*(1 - marginleft)/Ncol; % left pos
    yi = marginbottom + (Nrow - rowindex(ifig) + intermarginbottom)*(1 - marginbottom)/Nrow; % bottom pos
    w = (1 - intermarginright)*(1 - marginleft)/Ncol; % width
    h = (1 - intermargintop - intermarginbottom)*(1 - marginbottom)/Nrow; % height
    set(gca,'Position',[xj yi w h],'FontName',fontname,'FontSize',FS)
end

if PrintFig == 1
    set(gcf,'PaperUnits','inches')
    set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
    print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_speed_' filename '.pdf'])
else
    set(gcf,'Position',[10 10 1000 700])
end



%% Accumulation, inflow and outflow
%--------------------------------------------------------------------------

FS = 9; % font size
FS1 = 11; % title font size
figindex = 'abcdefghijklmnopqrstuvwxyz';
PrintFig = 0; % 1: print the figure (save as pdf in the img folder) / 0: do not print
filename = Results(1).Name; % name for printing

% Plot options
ResList = 1:NumRes; % list of reservoirs
RoutesList = []; % list of routes, put to [] for not plotting route states
ResuList = [1 3]; % list of results
TimeRange = [0 Simulation.Duration]; % [s]
AccRange = [0 1500]; % [veh]
FlowRange = [0 3]; % [veh/s]
PlotResTotalVal = 1; % 1: plot reservoir total states / 0: plot route states only
PlotResInternalDyn = 0; % 1: add plots of internal state (debug) / 0: do not add anything (better when plotting several results)
Nfig = length(ResList); % number of subfigures
Nplot = length(RoutesList); % number of plots per subfigure
Nresu = length(ResuList); % number of results

% Color and line options
sline = arrayextension(line_perso,Nresu,'column'); % line styles for results
LW = 1*ones(1,Nresu); % line widths for results
clight = linspace(0,(Nresu-1)/Nresu,Nresu); % color brightness (0 to 1) for results
cmap = arrayextension(cmap_perso2,Nplot,'row'); % colors for routes
color0 = [0 0 0]; % color for total values (sum over the routes)

% Figure and subfigure options
Ncol = 3; % number of columns in the figure
Nrow = (floor(Nfig/Ncol) < Nfig/Ncol).*(floor(Nfig/Ncol) + 1) + (floor(Nfig/Ncol) == Nfig/Ncol).*floor(Nfig/Ncol);
[colindex, rowindex] = ind2sub([Ncol Nrow],1:Nfig);
marginleft = 0.1; % figure relative left margin
marginbottom = 0.2; % figure relative bottom margin
intermarginright = 0.15; % subplot relative right margin
intermarginbottom = 0.05; % subplot relative bottom margin
intermargintop = 0.2; % subplot relative top margin
figwidth = 6.5; % whole figure width [inches]
figheight = Nrow*2.5; % whole figure height [inches]

% Accumulation
%--------------------------------------------------------------------------
figure
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    hp01 = zeros(1,Nresu);
    hp1 = zeros(Nresu,Nplot);
    strleg01 = cellstr(int2str(zeros(Nresu,1)));
    strleg1 = cellstr(int2str(zeros(Nplot,Nresu)));
    hold on
    for iresu = 1:Nresu
        Temp_t = Results(ResuList(iresu)).SimulTime;
        Temp_ntot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Acc;
        if PlotResTotalVal == 1
            hp01(iresu) = plot(Temp_t,Temp_ntot,...
                'linestyle',sline{iresu},'color',lightencolor(color0,clight(iresu)),'linewidth',LW(iresu));
            strleg01{iresu} = [Results(ResuList(iresu)).Name2 ' Total'];
        end
        for iplot = 1:Nplot
            i_r = Results(ResuList(iresu)).Route(RoutesList(iplot)).ResRouteIndex(ResList(ifig));
            if i_r > 0
                Temp_n = Results(ResuList(iresu)).Reservoir(ResList(ifig)).AccPerRoute(i_r,:);
            else
                Temp_n = zeros(1,length(Temp_t));
            end
            hp1(iresu,iplot) = plot(Temp_t,Temp_n,...
                'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
            strleg1{iplot,iresu} = [Results(ResuList(iresu)).Name2 ' Route ' int2str(RoutesList(iplot))];
        end
    end
    hold off
    
    axis([TimeRange AccRange])
    title(['(' figindex(ifig) ') - \itR_{\rm' int2str(ResList(ifig)) '}'],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    if ifig + Ncol > Nfig
        xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        if isempty(RoutesList)
            ylabel('accumulation \itn^r(t) \rm[veh]','FontName',fontname,'FontSize',FS)
        else
            ylabel('accumulation \itn_p^r(t) \rm[veh]','FontName',fontname,'FontSize',FS)
        end
    end
    xj = marginleft + (colindex(ifig) - 1)*(1 - marginleft)/Ncol; % left pos
    yi = marginbottom + (Nrow - rowindex(ifig) + intermarginbottom)*(1 - marginbottom)/Nrow; % bottom pos
    w = (1 - intermarginright)*(1 - marginleft)/Ncol; % width
    h = (1 - intermargintop - intermarginbottom)*(1 - marginbottom)/Nrow; % height
    set(gca,'Position',[xj yi w h],'FontName',fontname,'FontSize',FS)
    if ifig == 1
        hp = [];
        strleg = [];
        for iresu = 1:Nresu
            if PlotResTotalVal == 1
                if ~isempty(hp1)
                    hp = [hp hp01(iresu) hp1(iresu,:)];
                    strleg = [strleg; strleg01(iresu); strleg1(:,iresu)];
                else
                    hp = [hp hp01(iresu)];
                    strleg = [strleg; strleg01(iresu)];
                end
            else
                if ~isempty(hp1)
                    hp = [hp hp1(iresu,:)];
                    strleg = [strleg; strleg1(:,iresu)];
                end
            end
        end
        hleg = gridLegend(hp,Nresu,strleg,'Location','North','FontName',fontname,'FontSize',FS);
        legpos = get(hleg,'Position');
        set(hleg,'Position',[xj 0.01 legpos(3) legpos(4)])
    end
end

if PrintFig == 1
    set(gcf,'PaperUnits','inches')
    set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
    print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_acc_' filename '.pdf'])
else
    set(gcf,'Position',[10 10 1000 700])
end


% Inflow
%--------------------------------------------------------------------------
figure
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    hp01 = zeros(1,Nresu);
    hp1 = zeros(Nresu,Nplot);
    strleg01 = cellstr(int2str(zeros(Nresu,1)));
    strleg1 = cellstr(int2str(zeros(Nplot,Nresu)));
    hold on
    for iresu = 1:Nresu
        Temp_t = Results(ResuList(iresu)).SimulTime;
        Temp_qintot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Inflow;
        Temp_V = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MeanSpeed;
        Temp_ntot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Acc;
        Temp_Lavg = Results(ResuList(iresu)).Reservoir(ResList(ifig)).AvgTripLength;
        Temp_Pc = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MaxProd;
        if PlotResTotalVal == 1
            hp01(iresu) = plot(Temp_t,Temp_qintot,...
                'linestyle',sline{iresu},'color',lightencolor(color0,clight(iresu)),'linewidth',LW(iresu));
            strleg01{iresu} = [Results(ResuList(iresu)).Name2 ' Total'];
            if PlotResInternalDyn == 1
                plot(Temp_t,Temp_ntot./Temp_Lavg.*Temp_V,...
                    'linestyle','--','color',lightencolor(color0,0.4),'linewidth',3);
                plot(Temp_t,Temp_Pc./Temp_Lavg,...
                    'linestyle',':','color',lightencolor(color0,0.6),'linewidth',5);
            end
        end
        for iplot = 1:Nplot
            i_r = Results(ResuList(iresu)).Route(RoutesList(iplot)).ResRouteIndex(ResList(ifig));
            if i_r > 0
                Temp_qin = Results(ResuList(iresu)).Reservoir(ResList(ifig)).InflowPerRoute(i_r,:);
                Temp_n = Results(ResuList(iresu)).Reservoir(ResList(ifig)).AccPerRoute(i_r,:);
                Temp_Ltrip = Results(ResuList(iresu)).Reservoir(ResList(ifig)).TripLengthPerRoute(i_r);
            else
                Temp_qin = zeros(1,length(Temp_t));
                Temp_n = zeros(1,length(Temp_t));
                Temp_Ltrip = zeros(1,length(Temp_t));
            end
            hp1(iresu,iplot) = plot(Temp_t,Temp_qin,...
                'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
            strleg1{iplot,iresu} = [Results(ResuList(iresu)).Name2 ' Route ' int2str(RoutesList(iplot))];
            if PlotResInternalDyn == 1
                plot(Temp_t,Temp_n./Temp_Ltrip.*Temp_V,...
                    'linestyle','--','color',lightencolor(cmap(iplot,:),0.4),'linewidth',3);
                plot(Temp_t,Temp_n./Temp_ntot.*Temp_Pc./Temp_Ltrip,...
                    'linestyle',':','color',lightencolor(cmap(iplot,:),0.6),'linewidth',5);
            end
        end
    end
    hold off
    
    axis([TimeRange FlowRange])
    title(['(' figindex(ifig) ') - \itR_{\rm' int2str(ResList(ifig)) '}'],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    if ifig + Ncol > Nfig
        xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        if isempty(RoutesList)
            ylabel('inflow \itq_{\rmin}^r(t) \rm[veh/s]','FontName',fontname,'FontSize',FS)
        else
            ylabel('inflow \itq_{\rmin,\itp}^r(t) \rm[veh/s]','FontName',fontname,'FontSize',FS)
        end
    end
    xj = marginleft + (colindex(ifig) - 1)*(1 - marginleft)/Ncol; % left pos
    yi = marginbottom + (Nrow - rowindex(ifig) + intermarginbottom)*(1 - marginbottom)/Nrow; % bottom pos
    w = (1 - intermarginright)*(1 - marginleft)/Ncol; % width
    h = (1 - intermargintop - intermarginbottom)*(1 - marginbottom)/Nrow; % height
    set(gca,'Position',[xj yi w h],'FontName',fontname,'FontSize',FS)
    if ifig == 1
        hp = [];
        strleg = [];
        for iresu = 1:Nresu
            if PlotResTotalVal == 1
                if ~isempty(hp1)
                    hp = [hp hp01(iresu) hp1(iresu,:)];
                    strleg = [strleg; strleg01(iresu); strleg1(:,iresu)];
                else
                    hp = [hp hp01(iresu)];
                    strleg = [strleg; strleg01(iresu)];
                end
            else
                if ~isempty(hp1)
                    hp = [hp hp1(iresu,:)];
                    strleg = [strleg; strleg1(:,iresu)];
                end
            end
        end
        hleg = gridLegend(hp,Nresu,strleg,'Location','North','FontName',fontname,'FontSize',FS);
        legpos = get(hleg,'Position');
        set(hleg,'Position',[xj 0.01 legpos(3) legpos(4)])
    end
end

if PrintFig == 1
    set(gcf,'PaperUnits','inches')
    set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
    print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_inflow_' filename '.pdf'])
else
    set(gcf,'Position',[10 10 1000 700])
end


% Outflow
%--------------------------------------------------------------------------
figure
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    hp01 = zeros(1,Nresu);
    hp1 = zeros(Nresu,Nplot);
    strleg01 = cellstr(int2str(zeros(Nresu,1)));
    strleg1 = cellstr(int2str(zeros(Nplot,Nresu)));
    hold on
    for iresu = 1:Nresu
        Temp_t = Results(ResuList(iresu)).SimulTime;
        Temp_qouttot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Outflow;
        Temp_V = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MeanSpeed;
        Temp_ntot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Acc;
        Temp_Lavg = Results(ResuList(iresu)).Reservoir(ResList(ifig)).AvgTripLength;
        Temp_Pc = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MaxProd;
        if PlotResTotalVal == 1
            hp01(iresu) = plot(Temp_t,Temp_qouttot,...
                'linestyle',sline{iresu},'color',lightencolor(color0,clight(iresu)),'linewidth',LW(iresu));
            strleg01{iresu} = [Results(ResuList(iresu)).Name2 ' Total'];
            if PlotResInternalDyn == 1
                plot(Temp_t,Temp_ntot./Temp_Lavg.*Temp_V,...
                    'linestyle','--','color',lightencolor(color0,0.4),'linewidth',3);
                plot(Temp_t,Temp_Pc./Temp_Lavg,...
                    'linestyle',':','color',lightencolor(color0,0.6),'linewidth',5);
            end
        end
        for iplot = 1:Nplot
            i_r = Results(ResuList(iresu)).Route(RoutesList(iplot)).ResRouteIndex(ResList(ifig));
            if i_r > 0
                Temp_qout = Results(ResuList(iresu)).Reservoir(ResList(ifig)).OutflowPerRoute(i_r,:);
                Temp_n = Results(ResuList(iresu)).Reservoir(ResList(ifig)).AccPerRoute(i_r,:);
                Temp_Ltrip = Results(ResuList(iresu)).Reservoir(ResList(ifig)).TripLengthPerRoute(i_r);
            else
                Temp_qout = zeros(1,length(Temp_t));
                Temp_n = zeros(1,length(Temp_t));
                Temp_Ltrip = zeros(1,length(Temp_t));
            end
            hp1(iresu,iplot) = plot(Temp_t,Temp_qout,...
                'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
            strleg1{iplot,iresu} = [Results(ResuList(iresu)).Name2 ' Route ' int2str(RoutesList(iplot))];
            if PlotResInternalDyn == 1
                plot(Temp_t,Temp_n./Temp_Ltrip.*Temp_V,...
                    'linestyle','--','color',lightencolor(cmap(iplot,:),0.4),'linewidth',3);
                plot(Temp_t,Temp_n./Temp_ntot.*Temp_Pc./Temp_Ltrip,...
                    'linestyle',':','color',lightencolor(cmap(iplot,:),0.6),'linewidth',5);
            end
        end
    end
    hold off
    
    axis([TimeRange FlowRange])
    title(['(' figindex(ifig) ') - \itR_{\rm' int2str(ResList(ifig)) '}'],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    if ifig + Ncol > Nfig
        xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        if isempty(RoutesList)
            ylabel('outflow \itq_{\rmout}^r(t) \rm[veh/s]','FontName',fontname,'FontSize',FS)
        else
            ylabel('outflow \itq_{\rmout,\itp}^r(t) \rm[veh/s]','FontName',fontname,'FontSize',FS)
        end
    end
    xj = marginleft + (colindex(ifig) - 1)*(1 - marginleft)/Ncol; % left pos
    yi = marginbottom + (Nrow - rowindex(ifig) + intermarginbottom)*(1 - marginbottom)/Nrow; % bottom pos
    w = (1 - intermarginright)*(1 - marginleft)/Ncol; % width
    h = (1 - intermargintop - intermarginbottom)*(1 - marginbottom)/Nrow; % height
    set(gca,'Position',[xj yi w h],'FontName',fontname,'FontSize',FS)
    if ifig == 1
        hp = [];
        strleg = [];
        for iresu = 1:Nresu
            if PlotResTotalVal == 1
                if ~isempty(hp1)
                    hp = [hp hp01(iresu) hp1(iresu,:)];
                    strleg = [strleg; strleg01(iresu); strleg1(:,iresu)];
                else
                    hp = [hp hp01(iresu)];
                    strleg = [strleg; strleg01(iresu)];
                end
            else
                if ~isempty(hp1)
                    hp = [hp hp1(iresu,:)];
                    strleg = [strleg; strleg1(:,iresu)];
                end
            end
        end
        hleg = gridLegend(hp,Nresu,strleg,'Location','North','FontName',fontname,'FontSize',FS);
        legpos = get(hleg,'Position');
        set(hleg,'Position',[xj 0.01 legpos(3) legpos(4)])
    end
end

if PrintFig == 1
    set(gcf,'PaperUnits','inches')
    set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
    print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_outflow_' filename '.pdf'])
else
    set(gcf,'Position',[10 10 1000 700])
end


%% Total accumulation and inflow/outflow in reservoirs
%--------------------------------------------------------------------------

FS = 16; % font size

% Plot options
ResList = 1:NumRes; % list of reservoirs
ResuList = [1 2]; % list of results
TimeRange = [0 Simulation.Duration]; % [s]
AccRange = [0 1000]; % [veh]
FlowRange = [0 3]; % [veh/s]
filename = 'test'; % name for printing
Nplot = length(ResList); % number of plots per figure
Nresu = length(ResuList); % number of results

% Color and line options
sline = arrayextension(line_perso,Nresu,'column'); % line styles for results
LW = 2*ones(1,Nresu); % line widths for results
clight = linspace(0,(Nresu-1)/Nresu,Nresu); % color brightness (0 to 1) for results
cmap = arrayextension(cmap_perso2,Nplot,'row'); % colors for reservoirs

% Figure and subfigure options
figwidth = 6.5; % whole figure width [inches]
figheight = 2.5; % whole figure height [inches]

% Accumulation
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nplot*Nresu,1)));
hp1 = zeros(1,Nplot*Nresu);
hold on
i = 1;
for iresu = 1:Nresu
    Temp_t = Results(ResuList(iresu)).SimulTime;
    for iplot = 1:Nplot
        Temp_n = Results(ResuList(iresu)).Reservoir(ResList(iplot)).Acc;
        hp1(i) = plot(Temp_t,Temp_n,...
            'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
        strleg{i} = [Results(ResuList(iresu)).Name2 ' \itR_{\rm' int2str(ResList(iplot)) '}'];
        i = i + 1;
    end
end
hold off
axis([TimeRange AccRange])
xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
ylabel('accumulation \itn^r(t) \rm[veh]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_acctot_' filename '.pdf'])

% Inflow
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nplot*Nresu,1)));
hp1 = zeros(1,Nplot*Nresu);
hold on
i = 1;
for iresu = 1:Nresu
    Temp_t = Results(ResuList(iresu)).SimulTime;
    for iplot = 1:Nplot
        Temp_qin = Results(ResuList(iresu)).Reservoir(ResList(iplot)).Inflow;
        hp1(i) = plot(Temp_t,Temp_qin,...
            'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
        strleg{i} = [Results(ResuList(iresu)).Name2 ' \itR_{\rm' int2str(ResList(iplot)) '}'];
        i = i + 1;
    end
end
hold off
axis([TimeRange FlowRange])
xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
ylabel('inflow \itq_{\rmin}^r(t) \rm[veh/s]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_inflowtot_' filename '.pdf'])

% Outflow
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nplot*Nresu,1)));
hp1 = zeros(1,Nplot*Nresu);
hold on
i = 1;
for iresu = 1:Nresu
    Temp_t = Results(ResuList(iresu)).SimulTime;
    for iplot = 1:Nplot
        Temp_qout = Results(ResuList(iresu)).Reservoir(ResList(iplot)).Outflow;
        hp1(i) = plot(Temp_t,Temp_qout,...
            'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
        strleg{i} = [Results(ResuList(iresu)).Name2 ' \itR_{\rm' int2str(ResList(iplot)) '}'];
        i = i + 1;
    end
end
hold off
axis([TimeRange FlowRange])
xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
ylabel('outflow \itq_{\rmout}^r(t) \rm[veh/s]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_outflowtot_' filename '.pdf'])


%% Production, entry and exit (MFD plane)
%--------------------------------------------------------------------------

FS = 16; % font size
FS1 = 16; % title font size
figindex = 'abcdefghijklmnopqrstuvwxyz';

% Plot options
ResList = 1:NumRes; % list of reservoirs
RoutesList = []; % list of routes, put to [] for not plotting route states
ResuList = [1]; % list of results
AccRange = [0 3000]; % [veh]
ProdRange = [0 6000]; % [veh.m/s]
PlotResTotalVal = 1; % 1: plot reservoir total states / 0: plot route states only
PlotResInternalDyn = 0; % 1: add plots of internal state (debug) / 0: do not add anything (better when plotting several results)
filename = 'test'; % name for printing
Nfig = length(ResList); % number of subfigures
Nplot = length(RoutesList); % number of plots per subfigure
Nresu = length(ResuList); % number of results

% Color and line options
sline = arrayextension(line_perso,Nresu,'column'); % line styles for results
LW = 2*ones(1,Nresu); % line widths for results
clight = linspace(0,(Nresu-1)/Nresu,Nresu); % color brightness (0 to 1) for results
% clight = zeros(1,Nresu); % color brightness (0 to 1) for results
cmap = arrayextension(cmap_perso2,Nplot,'row'); % colors for routes
% cmap0 = arrayextension(cmap_perso2,Nresu,'row'); % color for total values (sum over the routes)
cmap0 = ones(Nresu,1)*[0 0 0]; % color for total values (sum over the routes)

% Figure and subfigure options
Ncol = 3; % number of columns in the figure
Nrow = (floor(Nfig/Ncol) < Nfig/Ncol).*(floor(Nfig/Ncol) + 1) + (floor(Nfig/Ncol) == Nfig/Ncol).*floor(Nfig/Ncol);
[colindex, rowindex] = ind2sub([Ncol Nrow],1:Nfig);
marginleft = 0.1; % figure relative left margin
marginbottom = 0.2; % figure relative bottom margin
intermarginright = 0.15; % subplot relative right margin
intermarginbottom = 0.05; % subplot relative bottom margin
intermargintop = 0.2; % subplot relative top margin
figwidth = 6.5; % whole figure width [inches]
figheight = Nrow*2.5; % whole figure height [inches]

% Entry production
%--------------------------------------------------------------------------
figure
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    hp01 = zeros(1,Nresu);
    hp1 = zeros(Nresu,Nplot);
    strleg01 = cellstr(int2str(zeros(Nresu,1)));
    strleg1 = cellstr(int2str(zeros(Nplot,Nresu)));
    hold on
    for iresu = 1:Nresu
        Temp_prodintot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).TripLengthPerRoute * ...
            Results(ResuList(iresu)).Reservoir(ResList(ifig)).InflowPerRoute;
        Temp_V = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MeanSpeed;
        Temp_ntot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Acc;
        Temp_nj = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MaxAcc;
        Temp_MFDparam = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MFDfctParam;
        Temp_Entryparam = Results(ResuList(iresu)).Reservoir(ResList(ifig)).EntryfctParam;
        n_ = linspace(0,Temp_nj,100);
        plot(n_,Simulation.MFDfct(n_,Temp_MFDparam),'-k','linewidth',0.5)
        plot(n_,Simulation.Entryfct(n_,Temp_Entryparam),'-k','linewidth',0.5)
        if PlotResTotalVal == 1
            hp01(iresu) = plot(Temp_ntot,Temp_prodintot,...
                'linestyle',sline{iresu},'color',lightencolor(cmap0(iresu,:),clight(iresu)),'linewidth',LW(iresu));
            strleg01{iresu} = [Results(ResuList(iresu)).Name2 ' Total'];
            if PlotResInternalDyn == 1
                plot(Temp_ntot,Temp_ntot.*Temp_V,...
                    'linestyle','--','color',lightencolor(cmap0(iresu,:),0.4),'linewidth',3);
            end
        end
        for iplot = 1:Nplot
            i_r = Results(ResuList(iresu)).Route(RoutesList(iplot)).ResRouteIndex(ResList(ifig));
            if i_r > 0
                Temp_qin = Results(ResuList(iresu)).Reservoir(ResList(ifig)).InflowPerRoute(i_r,:);
                Temp_n = Results(ResuList(iresu)).Reservoir(ResList(ifig)).AccPerRoute(i_r,:);
                Temp_Ltrip = Results(ResuList(iresu)).Reservoir(ResList(ifig)).TripLengthPerRoute(i_r);
            else
                Temp_qin = zeros(1,length(Temp_ntot));
                Temp_n = zeros(1,length(Temp_ntot));
                Temp_Ltrip = zeros(1,length(Temp_ntot));
            end
            hp1(iresu,iplot) = plot(Temp_n,Temp_Ltrip.*Temp_qin,...
                'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
            strleg1{iplot,iresu} = [Results(ResuList(iresu)).Name2 ' Route ' int2str(RoutesList(iplot))];
            if PlotResInternalDyn == 1
                plot(Temp_n,Temp_n.*Temp_V,...
                    'linestyle','--','color',lightencolor(cmap(iplot,:),0.4),'linewidth',3);
            end
        end
    end
    hold off
    
    axis([AccRange ProdRange])
    title(['(' figindex(ifig) ') - \itR_{\rm' int2str(ResList(ifig)) '}'],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    if ifig + Ncol > Nfig
        if isempty(RoutesList)
            xlabel('accumulation \itn^r(t) \rm[veh]','FontName',fontname,'FontSize',FS)
        else
            xlabel('accumulation \itn_p^r(t) \rm[veh]','FontName',fontname,'FontSize',FS)
        end
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        if isempty(RoutesList)
            ylabel('production \itL^r(t) q_{\rmin}^r(t) \rm[veh.m/s]','FontName',fontname,'FontSize',FS)
        else
            ylabel('production \itL_p^r q_{\rmin,\itp}^r(t) \rm[veh.m/s]','FontName',fontname,'FontSize',FS)
        end
    end
    xj = marginleft + (colindex(ifig) - 1)*(1 - marginleft)/Ncol; % left pos
    yi = marginbottom + (Nrow - rowindex(ifig) + intermarginbottom)*(1 - marginbottom)/Nrow; % bottom pos
    w = (1 - intermarginright)*(1 - marginleft)/Ncol; % width
    h = (1 - intermargintop - intermarginbottom)*(1 - marginbottom)/Nrow; % height
    set(gca,'Position',[xj yi w h],'FontName',fontname,'FontSize',FS)
    if ifig == 1
        hp = [];
        strleg = [];
        for iresu = 1:Nresu
            if PlotResTotalVal == 1
                if ~isempty(hp1)
                    hp = [hp hp01(iresu) hp1(iresu,:)];
                    strleg = [strleg; strleg01(iresu); strleg1(:,iresu)];
                else
                    hp = [hp hp01(iresu)];
                    strleg = [strleg; strleg01(iresu)];
                end
            else
                if ~isempty(hp1)
                    hp = [hp hp1(iresu,:)];
                    strleg = [strleg; strleg1(:,iresu)];
                end
            end
        end
        hleg = gridLegend(hp,Nresu,strleg,'Location','North','FontName',fontname,'FontSize',FS);
        legpos = get(hleg,'Position');
        set(hleg,'Position',[xj 0.01 legpos(3) legpos(4)])
    end
end
set(gcf,'Position',[10 10 1000 700])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_inprod_' filename '.pdf'])

% Exit production
%--------------------------------------------------------------------------
figure
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    hp01 = zeros(1,Nresu);
    hp1 = zeros(Nresu,Nplot);
    strleg01 = cellstr(int2str(zeros(Nresu,1)));
    strleg1 = cellstr(int2str(zeros(Nplot,Nresu)));
    hold on
    for iresu = 1:Nresu
        Temp_prodouttot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).TripLengthPerRoute * ...
            Results(ResuList(iresu)).Reservoir(ResList(ifig)).OutflowPerRoute;
        Temp_V = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MeanSpeed;
        Temp_ntot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Acc;
        Temp_nj = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MaxAcc;
        Temp_MFDparam = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MFDfctParam;
        n_ = linspace(0,Temp_nj,100);
        plot(n_,Simulation.MFDfct(n_,Temp_MFDparam),'-k','linewidth',0.5)
        if PlotResTotalVal == 1
            hp01(iresu) = plot(Temp_ntot,Temp_prodouttot,...
                'linestyle',sline{iresu},'color',lightencolor(cmap0(iresu,:),clight(iresu)),'linewidth',LW(iresu));
            strleg01{iresu} = [Results(ResuList(iresu)).Name2 ' Total'];
            if PlotResInternalDyn == 1
                plot(Temp_ntot,Temp_ntot.*Temp_V,...
                    'linestyle','--','color',lightencolor(cmap0(iresu,:),0.4),'linewidth',3);
            end
        end
        for iplot = 1:Nplot
            i_r = Results(ResuList(iresu)).Route(RoutesList(iplot)).ResRouteIndex(ResList(ifig));
            if i_r > 0
                Temp_qout = Results(ResuList(iresu)).Reservoir(ResList(ifig)).OutflowPerRoute(i_r,:);
                Temp_n = Results(ResuList(iresu)).Reservoir(ResList(ifig)).AccPerRoute(i_r,:);
                Temp_Ltrip = Results(ResuList(iresu)).Reservoir(ResList(ifig)).TripLengthPerRoute(i_r);
            else
                Temp_qout = zeros(1,length(Temp_ntot));
                Temp_n = zeros(1,length(Temp_ntot));
                Temp_Ltrip = zeros(1,length(Temp_ntot));
            end
            hp1(iresu,iplot) = plot(Temp_n,Temp_Ltrip.*Temp_qout,...
                'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
            strleg1{iplot,iresu} = [Results(ResuList(iresu)).Name2 ' Route ' int2str(RoutesList(iplot))];
            if PlotResInternalDyn == 1
                plot(Temp_n,Temp_n.*Temp_V,...
                    'linestyle','--','color',lightencolor(cmap(iplot,:),0.4),'linewidth',3);
            end
        end
    end
    hold off
    
    axis([AccRange ProdRange])
    title(['(' figindex(ifig) ') - \itR_{\rm' int2str(ResList(ifig)) '}'],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    if ifig + Ncol > Nfig
        if isempty(RoutesList)
            xlabel('accumulation \itn^r(t) \rm[veh]','FontName',fontname,'FontSize',FS)
        else
            xlabel('accumulation \itn_p^r(t) \rm[veh]','FontName',fontname,'FontSize',FS)
        end
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        if isempty(RoutesList)
            ylabel('production \itL^r(t) q_{\rmout}^r(t) \rm[veh.m/s]','FontName',fontname,'FontSize',FS)
        else
            ylabel('production \itL_p^r q_{\rmout,\itp}^r(t) \rm[veh.m/s]','FontName',fontname,'FontSize',FS)
        end
    end
    xj = marginleft + (colindex(ifig) - 1)*(1 - marginleft)/Ncol; % left pos
    yi = marginbottom + (Nrow - rowindex(ifig) + intermarginbottom)*(1 - marginbottom)/Nrow; % bottom pos
    w = (1 - intermarginright)*(1 - marginleft)/Ncol; % width
    h = (1 - intermargintop - intermarginbottom)*(1 - marginbottom)/Nrow; % height
    set(gca,'Position',[xj yi w h],'FontName',fontname,'FontSize',FS)
    if ifig == 1
        hp = [];
        strleg = [];
        for iresu = 1:Nresu
            if PlotResTotalVal == 1
                if ~isempty(hp1)
                    hp = [hp hp01(iresu) hp1(iresu,:)];
                    strleg = [strleg; strleg01(iresu); strleg1(:,iresu)];
                else
                    hp = [hp hp01(iresu)];
                    strleg = [strleg; strleg01(iresu)];
                end
            else
                if ~isempty(hp1)
                    hp = [hp hp1(iresu,:)];
                    strleg = [strleg; strleg1(:,iresu)];
                end
            end
        end
        hleg = gridLegend(hp,Nresu,strleg,'Location','North','FontName',fontname,'FontSize',FS);
        legpos = get(hleg,'Position');
        set(hleg,'Position',[xj 0.01 legpos(3) legpos(4)])
    end
end
set(gcf,'Position',[10 10 1000 700])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_outprod_' filename '.pdf'])



%% Travel time per route
%--------------------------------------------------------------------------

FS = 16; % font size

% Plot options
RoutesList = [1 2 3 4]; % list of routes
ResuList = [1 2]; % list of results
TimeRange = [0 Simulation.Duration]; % [s]
TTRange = [0 1000]; % [s]
PlotFreeflowTT = 0; % 1: add free-flow TT / 0: do not add (better when plotting several results)
filename = 'test'; % name for printing
Nplot = length(RoutesList); % number of plots per figure
Nresu = length(ResuList); % number of results

% Color and line options
sline = arrayextension(line_perso,Nresu,'column'); % line styles for results
LW = 2*ones(1,Nresu); % line widths for results
clight = linspace(0,(Nresu-1)/Nresu,Nresu); % color brightness (0 to 1) for results
cmap = arrayextension(cmap_perso2,Nplot,'row'); % colors for routes

figure
strleg = cellstr(int2str(zeros(Nplot*Nresu,1)));
hp1 = zeros(1,Nplot*Nresu);
hold on
i = 1;
for iresu = 1:Nresu
    for iplot = 1:Nplot
        Temp_t = Results(ResuList(iresu)).SimulTime;
        Temp_TT = Results(ResuList(iresu)).Route(RoutesList(iplot)).TravelTime;
        Temp_TTfree = Results(ResuList(iresu)).Route(RoutesList(iplot)).FreeFlowTravelTime;
        if PlotFreeflowTT == 1
            plot([0 Simulation.Duration],Temp_TTfree*[1 1],...
                'linestyle','--','color',lightencolor(cmap(iplot,:),0.5),'linewidth',4);
        end
        hp1(i) = plot(Temp_t,Temp_TT,...
            'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
        strleg{i} = [Results(ResuList(iresu)).Name2 ' Route ' int2str(RoutesList(iplot)) ': [' int2str(Results(ResuList(iresu)).Route(RoutesList(iplot)).ResPath) ']'];
        i = i + 1;
    end
end
hold off
axis([TimeRange TTRange])
xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
ylabel('travel time \itT(t) \rm[s]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','best','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])

% set(gcf,'PaperSize',[30 15],'PaperPosition',[0 0 30 15])
% print('-dpng',['UserNetworks/' Results(1).Network '/img/TT_comp_' filename '.png'])


%% MFD of reservoirs, with data
%--------------------------------------------------------------------------

LW = 1; % line width
FS = 9; % font size
FS1 = 12; % title font size
MS = 2; % marker size
figindex = 'abcdefghijklmnopqrstuvwxyz';

load(['UserNetworks/' Results(iresu).Network '/networkdata/Grid_9res_traffic_SCref.mat'])
Ndata = length(MFDData{1}.n);

% Plot options
ResList = 1:NumRes; % list of reservoirs
iresu = 1; % result
AccRange = [0 3000]; % [veh]
ProdRange = [0 3500]; % [veh.m/s]
SpeedRange = [0 15]; % [m/s]
filename = 'SCref'; % name for printing
Nfig = length(ResList); % number of plots per figure

% Color and line options
color1 = [0 0 0];
color2 = 0.7*[1 1 1];
sline = line_perso;

% Figure and subfigure options
Ncol = 3; % number of columns in the figure
Nrow = (floor(Nfig/Ncol) < Nfig/Ncol).*(floor(Nfig/Ncol) + 1) + (floor(Nfig/Ncol) == Nfig/Ncol).*floor(Nfig/Ncol);
[colindex, rowindex] = ind2sub([Ncol Nrow],1:Nfig);
marginleft = 0.1; % figure relative left margin
marginbottom = 0.1; % figure relative bottom margin
intermarginright = 0.15; % subplot relative right margin
intermarginbottom = 0.05; % subplot relative bottom margin
intermargintop = 0.2; % subplot relative top margin
figwidth = 6.5; % whole figure width [inches]
figheight = Nrow*2; % whole figure height [inches]

% Prod-MFD
%--------------------------------------------------------------------------
figure
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    hold on
    Temp_nj = Results(iresu).Reservoir(ResList(ifig)).MaxAcc;
    Temp_nc = Results(iresu).Reservoir(ResList(ifig)).CritAcc;
    Temp_param = Results(iresu).Reservoir(ResList(ifig)).MFDfctParam;
    n_ = linspace(0,Temp_nj,200);
    hp2 = plot(MFDData{ifig}.n(1:Ndata),MFDData{ifig}.P(1:Ndata),...
        'linestyle','none','marker','o','markersize',MS,'MarkerFaceColor',color2,'color',color2,'linewidth',0.5*LW);
    hp1 = plot(n_,MFDfct(n_,Temp_param),...
        'linestyle',sline{1},'marker','none','markersize',MS,'color',color1,'linewidth',LW);
    hold off
    %axis([AccRange ProdRange])
    %xlim([0 0.9*Temp_nc])
    title(['(' figindex(ifig) ') - \itR_{\rm' int2str(ResList(ifig)) '}'],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    if ifig + Ncol > Nfig
        xlabel('accumulation \itn^r \rm[veh]','FontName',fontname,'FontSize',FS)
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        ylabel('prod. \itP^r \rm[veh.m/s]','FontName',fontname,'FontSize',FS)
    end
    if ifig == Nfig
        hleg = legend([hp1 hp2],'theoretical fit','empirical data');
        legpos = get(hleg,'Position');
        %set(hleg,'Position',[0.6 marginbottom legpos(3) legpos(4)],'FontName',fontname,'FontSize',FS)
        %set(hleg,'Position',[0.6 0.01 legpos(3) legpos(4)],'FontName',fontname,'FontSize',FS)
        set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
    end
    xj = marginleft + (colindex(ifig) - 1)*(1 - marginleft)/Ncol; % left pos
    yi = marginbottom + (Nrow - rowindex(ifig) + intermarginbottom)*(1 - marginbottom)/Nrow; % bottom pos
    w = (1 - intermarginright)*(1 - marginleft)/Ncol; % width
    h = (1 - intermargintop - intermarginbottom)*(1 - marginbottom)/Nrow; % height
    set(gca,'Position',[xj yi w h],'FontName',fontname,'FontSize',FS)
end
set(gcf,'Position',[10 10 1000 700])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(1).Network '_pMFD_' filename '.pdf'])


% Speed-MFD
%--------------------------------------------------------------------------
figure
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    hold on
    Temp_nj = Results(iresu).Reservoir(ResList(ifig)).MaxAcc;
    Temp_param = Results(iresu).Reservoir(ResList(ifig)).MFDfctParam;
    n_ = linspace(0,Temp_nj,200);
    hp2 = plot(MFDData{ifig}.n(1:Ndata),MFDData{ifig}.V(1:Ndata),...
        'linestyle','none','marker','o','markersize',MS,'MarkerFaceColor',color2,'color',color2,'linewidth',0.5*LW);
    hp1 = plot(n_,MFDfct(n_,Temp_param)./n_,...
        'linestyle',sline{1},'marker','none','markersize',MS,'color',color1,'linewidth',LW);
    hold off
    axis([AccRange SpeedRange])
    title(['(' figindex(ifig) ') - \itR_{\rm' int2str(ResList(ifig)) '}'],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    if ifig + Ncol > Nfig
        xlabel('accumulation \itn^r \rm[veh]','FontName',fontname,'FontSize',FS)
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        ylabel('mean speed \itV^r \rm[m/s]','FontName',fontname,'FontSize',FS)
    end
    if ifig == Nfig
        hleg = legend([hp1 hp2],'theoretical fit','empirical data');
        set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
    end
    xj = marginleft + (colindex(ifig) - 1)*(1 - marginleft)/Ncol; % left pos
    yi = marginbottom + (Nrow - rowindex(ifig) + intermarginbottom)*(1 - marginbottom)/Nrow; % bottom pos
    w = (1 - intermarginright)*(1 - marginleft)/Ncol; % width
    h = (1 - intermargintop - intermarginbottom)*(1 - marginbottom)/Nrow; % height
    set(gca,'Position',[xj yi w h],'FontName',fontname,'FontSize',FS)
end
set(gcf,'Position',[10 10 1000 700])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(1).Network '_vMFD_' filename '.pdf'])



%% MFD of reservoirs
%--------------------------------------------------------------------------

FS = 16; % font size

% Plot options
ResList = 1:NumRes; % list of reservoirs
iresu = 1; % result
AccRange = [0 3000]; % [veh]
ProdRange = [0 3500]; % [veh.m/s]
SpeedRange = [0 15]; % [m/s]
filename = 'test'; % name for printing
Nplot = length(ResList); % number of plots per figure
Nresu = length(ResuList); % number of results

% Color and line options
sline = arrayextension(line_perso,Nplot,'column'); % line styles for reservoirs
LW = linspace(2,5,Nplot); % line widths for reservoirs
cmap = arrayextension(cmap_perso2,Nplot,'row'); % colors for reservoirs


% Production-MFD
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nplot,1)));
hp1 = zeros(1,Nplot);
hold on
for iplot = 1:Nplot
    Temp_nj = Results(iresu).Reservoir(ResList(iplot)).MaxAcc;
    Temp_param = Results(iresu).Reservoir(ResList(iplot)).MFDfctParam;
    Temp_n = linspace(0,Temp_nj,100);
    hp1(iplot) = plot(Temp_n,MFDfct(Temp_n,Temp_param),...
        'linestyle',sline{iplot},'color',cmap(iplot,:),'linewidth',LW(iplot));
    strleg{iplot} = ['\itR_{\rm' int2str(ResList(iplot)) '}'];
end
hold off
axis([AccRange ProdRange])
xlabel('accumulation \itn^r \rm[veh]','FontName',fontname,'FontSize',FS)
ylabel('production \itP^r(n^r) \rm[veh.m/s]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])


% Speed-MFD
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nplot,1)));
hp1 = zeros(1,Nplot);
hold on
for iplot = 1:Nplot
    Temp_nj = Results(iresu).Reservoir(ResList(iplot)).MaxAcc;
    Temp_param = Results(iresu).Reservoir(ResList(iplot)).MFDfctParam;
    Temp_n = linspace(0,Temp_nj,100);
    hp1(iplot) = plot(Temp_n,MFDfct(Temp_n,Temp_param)./Temp_n,...
        'linestyle',sline{iplot},'color',cmap(iplot,:),'linewidth',LW(iplot));
    strleg{iplot} = ['\itR_{\rm' int2str(ResList(iplot)) '}'];
end
hold off
axis([AccRange SpeedRange])
xlabel('accumulation \itn^r \rm[veh]','FontName',fontname,'FontSize',FS)
ylabel('mean speed \itV^r(n^r) \rm[m/s]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])


%% Video of reservoir states
%--------------------------------------------------------------------------

% Result to plot
iresu = 2;
ResList = 1:NumRes; % list of reservoirs
RoutesList = [1 2]; % list of routes
filename = 'test';

% Time frame
tfinal = Simulation.Duration; % final time [s]
Dt_frame = 200; % time step between two frames [s]
Nframes = floor(tfinal/Dt_frame) + 1; % number of frames
t_frame = 0:Dt_frame:tfinal;

% Plot options
opts.fontname = 'Arial';
opts.fontsize = 18;
opts.linewidth = 0.8;
opts.colormap = cmap_perso2;
tref = datenum([2017 1 1 5 00 00]);

% Make a succession of images
% mkdir(['UserNetworks/' Results(1).Network '/img'],filename)

% Make a video: method 1 (old)
% strfile = ['UserNetworks/' Results(1).Network '/img/' filename '.avi'];
% aviobj = avifile(strfile, 'fps',5);

% Make a video: method 2 (old)
% strfile = ['UserNetworks/' Results(1).Network '/img/' filename '.avi'];
% mov(1:Nframes) = struct('cdata',[], 'colormap',[]);

% Make a video: method 3 (recommended)
% strfile = ['UserNetworks/' Results(1).Network '/img/' filename '.avi'];
% video = VideoWriter(strfile);
% video.FrameRate = 5; % video fps
% open(video)

figure
set(gcf,'Position',[50 50 800 500])
set(gcf,'PaperSize',[30 18],'PaperPosition',[0 0 30 18])
% set(gcf,'PaperSize',[15 11],'PaperPosition',[0 0 15 11])
set(gcf,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
for i = 1:Nframes
    clf
    t0 = t_frame(i);
    %opts.title = datestr(tref+t0./(24*3600),'HH:MM');
    opts.title = ['\bft = ' num2str(t0) ' s'];
    opts.showleg = 1;
    plotResBallAcc(t0,Results(iresu).Reservoir,ResList,Results(iresu).SimulTime,0.5,[1.5 1.5],'trafficolor',opts)
    %plotResBallAccPerRoute(t0,Results(iresu).Reservoir,ResList,Results(iresu).Route,RoutesList,Results(iresu).SimulTime,0.5,[1.5 1.5],opts)
    %plotResNetAcc(t0,[],Results(iresu).Reservoir,ResList,Results(iresu).SimulTime,'trafficolor',opts)
    %plotResNetAccPerRoute(t0,[],Results(iresu).Reservoir,ResList,Results(iresu).Route,RoutesList,Results(iresu).SimulTime,opts)
    
    % Succession of images
    %print('-dpdf','-painters',['UserNetworks/' Results(1).Network '/img/' filename '/img_' int2str(i) '.pdf'])
    
    % Video method 1 (old)
    %aviobj = addframe(aviobj, getframe(gcf));
    
    % Video method 2 (old)
    %mov(i) = getframe(gcf);
    
    % Video method 3 (recommended)
    %frame = getframe(gcf);
    %writeVideo(video,frame);
    
    pause(0.001)
    %pause
    
end

% Video method 1 (old)
% aviobj = close(aviobj);

% Video method 2 (old)
% movie2avi(mov, strfile, 'compression','Indeo5', 'fps',5);

% Video method 3 (recommended)
% close(video)



%% Video of reservoir states, 2 sub figures
%--------------------------------------------------------------------------

% Result to plot
ResList = 1:NumRes; % list of reservoirs
RoutesList = [1 2]; % list of routes
ResuList = [1 2]; % list of results
filename = 'test';

% Time frame
tfinal = Simulation.Duration; % final time [s]
Dt_frame = 200; % time step between two frames [s]
Nframes = floor(tfinal/Dt_frame) + 1; % number of frames
t_frame = 0:Dt_frame:tfinal;

% Plot options
opts.fontname = 'Arial';
opts.fontsize = 18;
opts.linewidth = 0.8;
opts.colormap = cmap_perso2;
tref = datenum([2017 1 1 5 00 00]);

% Make a succession of images
% mkdir(['UserNetworks/' Results(1).Network '/img'],filename)

% Make a video: method 1 (old)
% strfile = ['UserNetworks/' Results(1).Network '/img/' filename '.avi'];
% aviobj = avifile(strfile, 'fps',5);

% Make a video: method 2 (old)
% strfile = ['UserNetworks/' Results(1).Network '/img/' filename '.avi'];
% mov(1:Nframes) = struct('cdata',[], 'colormap',[]);

% Make a video: method 3 (recommended)
% strfile = ['UserNetworks/' Results(1).Network '/img/' filename '.avi'];
% video = VideoWriter(strfile);
% video.FrameRate = 5; % video fps
% open(video)

Nfig = length(ResuList); % number of subfigures
Ncol = 2; % number of columns in the figure
Nrow = (floor(Nfig/Ncol) < Nfig/Ncol).*(floor(Nfig/Ncol) + 1) + (floor(Nfig/Ncol) == Nfig/Ncol).*floor(Nfig/Ncol);

figure
set(gcf,'Position',[50 50 800 500])
set(gcf,'PaperSize',[30 18],'PaperPosition',[0 0 30 18])
% set(gcf,'PaperSize',[15 11],'PaperPosition',[0 0 15 11])
set(gcf,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
for i = 1:Nframes
    clf
    t0 = t_frame(i);
    for ifig = 1:Nfig
        subplot(Nrow,Ncol,ifig)
        %opts.title = datestr(tref+t0./(24*3600),'HH:MM');
        opts.title = ['\bft = ' num2str(t0) ' s'];
        opts.showleg = 1;
        %plotResBallAcc(t0,Results(ResuList(ifig)).Reservoir,ResList,Results(ResuList(ifig)).SimulTime,0.5,[1.5 1.5],'trafficolor',opts)
        %plotResBallAccPerRoute(t0,Results(ResuList(ifig)).Reservoir,ResList,Results(ResuList(ifig)).Route,RoutesList,Results(ResuList(ifig)).SimulTime,0.5,[1.5 1.5],opts)
        %plotResNetAcc(t0,[],Results(ResuList(ifig)).Reservoir,ResList,Results(ResuList(ifig)).SimulTime,'trafficolor',opts)
        plotResNetAccPerRoute(t0,[],Results(ResuList(ifig)).Reservoir,ResList,Results(ResuList(ifig)).Route,RoutesList,Results(ResuList(ifig)).SimulTime,opts)
        if ifig == 1
            set(gca,'Position',[0 0.01 0.47 0.95])
        elseif ifig == 2
            set(gca,'Position',[0.55 0.01 0.47 0.95])
        end
    end
    % Succession of images
    %print('-dpdf','-painters',['UserNetworks/' Results(1).Network '/img/' filename '/img_' int2str(i) '.pdf'])
    
    % Video method 1 (old)
    %aviobj = addframe(aviobj, getframe(gcf));
    
    % Video method 2 (old)
    %mov(i) = getframe(gcf);
    
    % Video method 3 (recommended)
    %frame = getframe(gcf);
    %writeVideo(video,frame);
    
    pause(0.001)
    
end

% Video method 1 (old)
% aviobj = close(aviobj);

% Video method 2 (old)
% movie2avi(mov, strfile, 'compression','Indeo5', 'fps',5);

% Video method 3 (recommended)
% close(video)


