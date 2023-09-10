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
            if nargin > 0
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
            else
                initialize(obj)
            end
        end
        
        function initialize(obj)
            obj.Enabled = 'Y';
            obj.OptimizationSweeps = {};
            obj.OptimizationGoals = {};
        end
        
        function enable(obj)
            obj.Enabled = 'Y';
        end
        
        function disable(obj)
            obj.Enabled = 'N';
        end
        
        function toggle(obj)
            if obj.Enabled == 'Y'
                obj.Enabled = 'N';
            else
                obj.Enabled = 'Y';
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetGoalSetBlock();
            SonnetClone(obj,aNewObject);
        end
        
        function index = addSweep(obj, sweep)   % add frequency sweep <sweep> to this object
            index = length(obj.OptimizationSweeps)+1;
            sweep.AdaptiveEnabled = 'X';
            obj.OptimizationSweeps{index} = sweep;
        end
        
        function index = addGoal(obj, goal)     % add goal <goal> to this object
            index = length(obj.OptimizationGoals)+1;
            obj.OptimizationGoals{index} = goal;
        end
        
        function index = rmSweep(obj, index)   % remove sweep at <index> from this object
            obj.OptimizationSweeps = [obj.OptimizationSweeps(1:index-1), obj.OptimizationSweeps(index+1:end)];
        end
        
        function index = rmGoal(obj, index)     % remove goal at <index> from this object
            obj.OptimizationGoals = [obj.OptimizationGoals(1:index-1), obj.OptimizationGoals(index+1:end)];
        end
        
        function index = setSweep(obj, index, sweep)    
            % set the frequency sweep at the specified index to sweep, overwriting
            % what is already there (if anything)
            if ~contains(class(sweep), 'SonnetFrequency')
            	error('sweep must be a SonnetFrequency* object')
            end
            obj.OptimizationSweeps{index} = sweep;
        end
        
        function index = setGoal(obj, index, goal)      
            % set the frequency sweep at the specified index to sweep, overwriting
            % what is already there (if anything)
            if ~contains(class(sweep), 'SonnetGoal')
            	error('sweep must be a SonnetGoal* object')
            end
            obj.OptimizationGoals{index} = goal;
        end
        
        function sweep = getSweep(obj, index)    
            % get the frequency sweep at the specified index
            sweep = obj.OptimizationSweeps{index};
        end
        
        function goal = getGoal(obj, index)      
             % get the goal at the specified index
            goal = obj.OptimizationGoals{index};
        end
        
        function aString = stringSignature(obj, ~)
            %stringSignature Convert this class to a string as it would
            %appear in a file
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

