classdef SonnetSweepSetBlock  < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class defines sweep set subsection of the VARSWP portion of a SONNET project file.
    % This class is a container for the sweep set information that is obtained
    % from the SONNET project file.
    %
    % SonnetLab, all included documentation, all included examples
    % and all other files (unless otherwise specified) are copyrighted by Sonnet Software
    % in 2011 with all rights reserved.
    %
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS". ANY AND
    % ALL EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED. UNDER NO CIRCUMSTANCES AND UNDER
    % NO LEGAL THEORY, TORT, CONTRACT, OR OTHERWISE, SHALL THE COPYWRITE HOLDERS,  CONTRIBUTORS,
    % MATLAB, OR SONNET SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
    % CONSEQUENTIAL DAMAGES OF ANY CHARACTER INCLUDING, WITHOUT LIMITATION, DAMAGES FOR LOSS OF
    % GOODWILL, WORK STOPPAGE, COMPUTER FAILURE OR MALFUNCTION, OR ANY AND ALL OTHER COMMERCIAL
    % DAMAGES OR LOSSES, OR FOR ANY DAMAGES EVEN IF THE COPYWRITE HOLDERS, CONTRIBUTORS, MATLAB,
    % OR SONNET SOFTWARE HAVE BEEN INFORMED OF THE POSSIBILITY OF SUCH DAMAGES, OR FOR ANY CLAIM
    % BY ANY OTHER PARTY.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties
	
		Enabled;
		ArrayOfSweeps;
		ArrayOfVariables;
	end
	
	methods
        function obj = SonnetSweepSetBlock(theFid)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % The constructor for SweepSetBlock.
            %     this will be passed the file ID from the
            %     SONNET project constructor.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            NumberOfSweepsInSet=0;   % Keeps track of the number of sweeps we have in the set
            NumberOfVariablesInSet=0;   % Keeps track of the number of variables we have in the set
            
            if nargin == 1
                
                aTempString=fscanf(theFid, '%s', 1);     % Try to read the the first type of the sweep from the file.
                                         %because SonnetVariableSweepBlock
                                         %already read ENABLED, pointer is
                                         %now before Y or N
                
                if ~isempty(aTempString)       % aTempString should now be Y or N, set ENABLED to that value
                    obj.Enabled = strtrim(aTempString);
                else                            % File pointer should be just after ENABLED
                                                % if it is empty, throw error
                    error("Attempted reading if sweep set is enabled, nothing found")
                end
                fgetl(theFid);      %Finish the line
                aTempString=fgetl(theFid);     % Try to read the the first type of the sweep from the file.
                while (1==1)                           % Keep looping until we get to the end of the block. This is when we have read in all the sweeps and read in 'END'
                    [aTempSplit, aTempRemain] = strtok(aTempString);       %split the string so it is easier to parse, strtrim for stray whitespaces in undefined sweeps
                    switch aTempSplit
                        case 'FREQ'
                            NumberOfSweepsInSet=NumberOfSweepsInSet+1;         % increment the sweep counter by one
                            obj.ArrayOfSweeps{NumberOfSweepsInSet}=SonnetFrequencyParse(aTempString);    	% construct the new sweep and store in the cell array
                        case 'VAR'
                            NumberOfVariablesInSet=NumberOfVariablesInSet+1;         % increment the sweep counter by one
                            obj.ArrayOfVariables{NumberOfVariablesInSet}=SonnetVariableDeclaration(aTempRemain);    	% construct the new sweep and store in the cell array
                        case 'END'			% check if we reached the end of the block
                            break;
                    end
                    aTempString=fgetl(theFid);    % read the next sweep name from the file, if it is END then we are done
                end
            else
                initialize(obj);
            end
        end
        
        function initialize(obj)
            obj.Enabled = 'Y';
            obj.ArrayOfSweeps = [];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetSweepSetBlock();
            SonnetClone(obj,aNewObject);
        end
        
        function index = addFrequencySweep(obj, sweep, enabled, adaptive, arg1, arg2, arg3)
            %addFrequencySweep adds a sweep to the sweep set. Can be passed the
            %arguments to create a sweep with the type as the sweep argument,
            %or a sweep itself (created for example by CreateSweep).
            %Will add a variable sweep if sweep is a SonnetVariableSweep
            %   Types and arguments are as follows:
            %
            %       SWEEP       StartFrequency,EndFrequency,StepFrequency
            %       ABS_ENTRY   StartFrequency,EndFrequency, FrequencyCount*
            %       DC          Mode**,Frequency**
            %       ESWEEP      StartFrequency,EndFrequency,AnalysisFrequencies
            %       LIST        List
            %       STEP        StepFrequency
            %
            %       *   For ABS_ENTRY:  FrequencyCount defaults to 300
            %       **  For a DC sweep: mode is either 'AUTO' for automatic or 'MAN' for manual.
            %       **  For a DC sweep: when mode is 'AUTO' the frequency does not need to
            %                           be supplied. The frequency is required when the DC
            %                           mode is manual.
            switch nargin
                case 2
                    if isa(sweep, 'SonnetVariableSweep')
                        obj.addVariableSweep(sweep);
                    elseif ~isa(sweep, 'SonnetFrequency')
                        error('sweep must be a SonnetFrequency* object')
                    end
                    index = length(obj.ArrayOfSweeps) + 1;
                    if sweep.AdaptiveEnabled == 'X'
                        sweep.AdaptiveEnabled = 'N';
                    end
                    obj.ArrayOfSweeps{index} = sweep;
                case 5
                    index = length(obj.ArrayOfSweeps) + 1;
                    obj.ArrayOfSweeps{index} = SonnetFrequencyParse(sweep, enabled, adaptive, arg1);
                case 6
                    index = length(obj.ArrayOfSweeps) + 1;
                    obj.ArrayOfSweeps{index} = SonnetFrequencyParse(sweep, enabled, adaptive, arg1, arg2);
                case 7
                    index = length(obj.ArrayOfSweeps) + 1;
                    obj.ArrayOfSweeps{index} = SonnetFrequencyParse(sweep, enabled, adaptive, arg1, arg2, arg3);
            end
        end
        
        function index = addSweep(obj, sweep, enabled, adaptive, arg1, arg2, arg3)
            %addSweep adds a sweep to the sweep set. Can be passed the
            %arguments to create a sweep with the type as the sweep argument,
            %or a sweep itself (created for example by CreateSweep).
            %Will add a variable sweep if sweep is a SonnetVariableSweep
            %   Types and arguments are as follows:
            %
            %       SWEEP       StartFrequency,EndFrequency,StepFrequency
            %       ABS_ENTRY   StartFrequency,EndFrequency, FrequencyCount*
            %       DC          Mode**,Frequency**
            %       ESWEEP      StartFrequency,EndFrequency,AnalysisFrequencies
            %       LIST        List
            %       STEP        StepFrequency
            %
            %       *   For ABS_ENTRY:  FrequencyCount defaults to 300
            %       **  For a DC sweep: mode is either 'AUTO' for automatic or 'MAN' for manual.
            %       **  For a DC sweep: when mode is 'AUTO' the frequency does not need to
            %                           be supplied. The frequency is required when the DC
            %                           mode is manual.
            switch nargin
                case 2
                    if isa(sweep, 'SonnetVariableSweep')
                        obj.addVariableSweep(sweep);
                    elseif ~isa(sweep, 'SonnetFrequency')
                        error('sweep must be a SonnetFrequency* object')
                    end
                    index = length(obj.ArrayOfSweeps) + 1;
                    if sweep.AdaptiveEnabled == 'X'
                        sweep.AdaptiveEnabled = 'N';
                    end
                    obj.ArrayOfSweeps{index} = sweep;
                case 5
                    index = length(obj.ArrayOfSweeps) + 1;
                    obj.ArrayOfSweeps{index} = SonnetFrequencyParse(sweep, enabled, adaptive, arg1);
                case 6
                    index = length(obj.ArrayOfSweeps) + 1;
                    obj.ArrayOfSweeps{index} = SonnetFrequencyParse(sweep, enabled, adaptive, arg1, arg2);
                case 7
                    index = length(obj.ArrayOfSweeps) + 1;
                    obj.ArrayOfSweeps{index} = SonnetFrequencyParse(sweep, enabled, adaptive, arg1, arg2, arg3);
            end
        end
        
        function sweep = getSweep(obj, index)
            sweep = obj.ArrayOfSweeps{index};
        end
        
        function sweep = getFrequencySweep(obj, index)
            sweep = obj.ArrayOfSweeps{index};
        end
        
        function index = setSweep(obj, index, sweep)  
            % set the frequency sweep at the specified index to sweep, overwriting
            % what is already there (if anything)
            if ~is(sweep, 'SonnetFrequency')
            	error('sweep must be a SonnetFrequency* object')
            end
            obj.ArrayOfSweeps{index} = sweep;
        end
        
        function sweep = rmSweep(obj, start_index, end_index)  
            % remove the frequency sweep at the specified index
            % Arguments:
            %   1) start_index - the index of the first sweep to remove
            %   2) end_index - the index of the last sweep to
            %       remove, romoving all between start and end index, inclusive
            % WARNING: this will change the index of all higher-indexed
            %   frequency sweeps
            sweep = obj.ArrayOfSweeps(start_index:end_index);
            obj.ArrayOfSweeps = [obj.ArrayOfSweeps(1:start_index-1) obj.ArrayOfSweeps(end_index+1:end)];
        end
        
        function sweep = getVariableSweep(obj, index)
            sweep = obj.ArrayOfVariables{index};
        end
        
        function index = addVariableSweep(obj, variable_sweep)
            %addVariableSweep adds a variable sweep to the sweep set.
            %Passed the result of a call to, for example, CreateVariableSweep
            %Any SonnetVariableSweep object will work
            if ~isa(variable_sweep, 'SonnetVariableSweep')
                error('variable_sweep must be a SonnetVariableSweep')
            end
            index = length(obj.ArrayOfVariables) + 1;
            obj.ArrayOfVariables{index} = variable_sweep;
        end
        
        function activateVariableSweep(obj, index)
            %activateVariableSweep activates a variable sweep in the sweep set.
            obj.ArrayOfVariables{index}.enable();
        end
        
        function deactivateVariableSweep(obj, index)
            %activateVariableSweep deactivates a variable sweep in the sweep set.
            obj.ArrayOfVariables{index}.disable();
        end
        
        function index = setVariableSweep(obj, index, sweep)  
            % set the variable sweep at the specified index to sweep, overwriting
            % what is already there (if anything)
            if ~isa(sweep, 'SonnetVariableSweep')
            	error('sweep must be a SonnetVariableSweep object')
            end
            obj.ArrayOfVariables{index} = sweep;
        end
        
        function sweep = rmVariableSweep(obj, start_index, end_index)
            % remove the variable sweeps between the specified start and end index
            % Arguments:
            %   1) start_index - the index of the first sweep to remove
            %   2) end_index - the index of the last sweep to
            %       remove, romoving all between start and end index, inclusive
            % WARNING: this will change the index of all higher-indexed
            %   variable sweeps
            if nargin == 2
                end_index = start_index;
            end
            sweep = obj.ArrayOfVariables(start_index:end_index);
            obj.ArrayOfVariables = [obj.ArrayOfVariables(1:start_index-1) obj.ArrayOfVariables(end_index+1:end)];
        end
        
        function enable(obj)
            obj.Enabled = 'Y';
        end
        
        function disable(obj)
            obj.Enabled = 'N';
        end
        
        function toggle(obj)
            if obj.Enabled == 'N'
                obj.Enabled = 'Y';
            else
                obj.Enabled = 'N';
            end
        end
        
        function aSignature = stringSignature(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Returns the data in the block as it would be printed to a
            % file, unformatted (contains raw identifiers like \n)
            % Does not print nicely to console, use stringSignaturef
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aSignature = ['ENABLED ' obj.Enabled '\n'];
            for iCounter = 1:length(obj.ArrayOfSweeps)
                aSignature = [aSignature obj.ArrayOfSweeps{iCounter}.stringSignature()]; %#ok<AGROW>
            end
            for iCounter = 1:length(obj.ArrayOfVariables)
                aSignature = [aSignature obj.ArrayOfVariables{iCounter}.stringSignature()]; %#ok<AGROW>
            end
            aSignature = sprintf([aSignature 'END\n']);
        end
        
        function aSignature = stringSignaturef(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Returns the data in the block as it would be printed to a
            % file, formatted (contains no raw identifiers like \n).
            % Prints nicely to console
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aSignature = stringSignature(obj);
            aSignature = sprintf(aSignature);
        end
        
        function aString = toString(obj)
            aString = stringSignature(obj);
        end
    end
end