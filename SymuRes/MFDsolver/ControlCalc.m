%% CONTROLCALC: Control action calculation
%--------------------------------------------------------------------------

itime = floor(Assignment.CurrentTime/TimeStep) + 1;
iperiod = Assignment.CurrentPeriodID;
Temp_Nperiods = length(Assignment.Periods) - 1;

NumControls = length(Simulation.Control);

for ictrl = 1:NumControls
    
    Temp_reslist = Simulation.Control(ictrl).ResID;
    Temp_nodelist = Simulation.Control(ictrl).MacroNodesID;
    Temp_minbound = Simulation.Control(ictrl).MinInflow;
    Temp_maxbound = Simulation.Control(ictrl).MaxInflow;
    
    if iperiod == 1
        
        % Initialization
        for inode = Temp_nodelist
            MacroNode(inode).Capacity.Time = Assignment.Periods(1:(end-1));
            MacroNode(inode).Capacity.Data = Temp_maxbound*ones(1,Temp_Nperiods);
        end
        
    else
        
        % PID controler
        %--------------------------------------------------------------------------
        if strcmp(Simulation.Control(ictrl).Type,'PID') && Simulation.Control(ictrl).Active == 1
            
            % Measured accumulation
            Temp_n = 0;
            if Simulation.Solver == 1
                for r = Temp_reslist
                    Temp_n = Temp_n + Reservoir(r).Acc(itime);
                end
            elseif Simulation.Solver == 2
                for r = Temp_reslist
                    Temp_n = Temp_n + 1/Simulation.TripbasedSimuFactor*Reservoir(r).CurrentAcc;
                end
            end
            
            % Accumulation errors
            %Temp_nerror = Temp_n - Simulation.Control(ictrl).SetPointAcc;
            Temp_nerror = Simulation.Control(ictrl).SetPointAcc - Temp_n;
            Temp_nintegral = Simulation.Control(ictrl).IntError + Temp_nerror*TimeStep;
            Temp_nderivative = (Temp_nerror - Simulation.Control(ictrl).PrevError)/TimeStep;
            Simulation.Control(ictrl).PrevError = Temp_nerror;
            Simulation.Control(ictrl).IntError = Temp_nintegral;
            
            % Gains
            Temp_Kp = Simulation.Control(ictrl).Param(1);
            Temp_Ki = Simulation.Control(ictrl).Param(2);
            Temp_Kd = Simulation.Control(ictrl).Param(3);
            
            % Control action
            Temp_ucontrol = Temp_Kp*Temp_nerror + Temp_Ki*Temp_nintegral + Temp_Kd*Temp_nderivative;
            
            % Macro node allowed flow update
            for inode = Temp_nodelist
                % Measured inflow through the node
                r = MacroNode(inode).ResID(2); % res receiving the macro node inflow
                i_n = 1;
                while inode ~= Reservoir(r).MacroNodesID(i_n)
                    i_n = i_n + 1;
                end
                %i_n = max([1 i_n-1]);
                Temp_indexes = Reservoir(r).NodeRoutesIndex{i_n}; % index of routes in r crossing inode
                if Simulation.Solver == 1
                    Temp_inflow = sum(Reservoir(r).InflowPerRoute(Temp_indexes,itime-1));
                elseif Simulation.Solver == 2
                    % TODO
                end
                
                % Allowed macro node flow
                %Temp_controlsupply = max([Temp_minbound (1 - Temp_ucontrol)*Temp_inflow]);
                Temp_controlsupply = Temp_ucontrol/length(Temp_nodelist);
                Temp_controlsupply = max([Temp_minbound Temp_controlsupply]);
                Temp_controlsupply = min([Temp_maxbound Temp_controlsupply]);
                MacroNode(inode).Capacity.Data(iperiod) = Temp_controlsupply;
            end
            
            % Print control info
            fprintf('%s %1.3f %s %1.3f \n','Error:',Temp_nerror,' control action:',Temp_ucontrol)
            fprintf('%s \n','Flow limit per macro node:')
            for inode = Temp_nodelist
                fprintf('%1.3f \t',MacroNode(inode).Capacity.Data(iperiod))
            end
            fprintf('%s\n',' ')
        end
        
    end
    
end


