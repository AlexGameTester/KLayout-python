classdef SonnetGoalSetBlock < handle
    %SonnetOptimizationSetBlock Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Enabled
        OptimizationSweeps
        OptimizationGoals
    end
    
    methods
        function obj = SonnetGoalSetBlock(theFid)
            aTempString = fgetl(theFid);
            if contains(aTempString, 'Y')       % if 'Y' not in the rest of the line, assume 'N'
                obj.Enabled = 'Y';
            else
                obj.Enabled = 'N';
            end
            
            NumberOfOptimizationSweeps = 0;     % keep track of sweeps
            NumberOfOptimizationGoals = 0;      % keep track of goals
            
            while 1==1
                aTempString = fscanf(theFid,' %s',1);
                switch aTempString
                    case 'FREQ'                                                                             % check if it is a sweep sweep
                        NumberOfOptimizationSweeps=NumberOfOptimizationSweeps+1;                            % Increase the counter for the number of sweeps we have in the block
                        obj.OptimizationSweeps{NumberOfOptimizationSweeps}=SonnetFrequencyParse(fgetl(theFid));     % construct the new sweep and store in the cell array

                    case 'GOAL'                                                                          % If it is the optimization line
                        NumberOfOptimizationGoals = NumberOfOptimizationGoals + 1;
                        obj.OptimizationGoals{NumberOfOptimizationGoals}=SonnetGoalParse(fgetl(theFid));     % construct the new goal and store in the cell array
                    case 'END'
                        break
                end
            end
        end
        
        function index = addSweep(obj, sweep)
            index = length(obj.OptimizationSweeps)+1;
            sweep.AdaptiveEnabled = 'X';
            obj.OptimizationSweeps{index} = sweep;
        end
        
        function aString = stringSignature(obj, ~)
            aString = ['ENABLED ' obj.Enabled '\n'];
            for iCounter =  1:length(obj.OptimizationSweeps)
                aString = [aString obj.OptimizationSweeps{iCounter}.stringSignature()];
            end
            for iCounter =  1:length(obj.OptimizationGoals)
                aString = [aString char(obj.OptimizationGoals{iCounter}.stringSignature())];
            end
            aString = sprintf([aString 'END\n']);
        end
        
        function writeObjectContents(obj, theFid, ~)
            fprintf(theFid, obj.stringSignature());
        end
    end
end

