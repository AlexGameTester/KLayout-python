classdef SonnetVariableSweepBlock  < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class defines the VARSWP portion of a SONNET project file.
    % This class is a container for the VARSWP information that is obtained
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
        
        SweepVariables;                 % This property determines whether parameter sweeps are on
        ArrayOfSweepSets;				% This property stores the sweep sets that we had in the file.
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = SonnetVariableSweepBlock(theFid)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % The constructor for VARSWP.
            %     the VARSWP will be passed the file ID from the
            %     SONNET project constructor.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            NumberOfSweepSetsInFile=0;   % Keeps track of the number of sweeps we have in the file
            
            initialize(obj);          % Initialize the values of the properties using the initializer function
            
            if nargin == 1
                aTempString=fscanf(theFid,'%s',1);     % Try to read the the first type of the sweep from the file.
                obj.SweepVariables = 0;
                while (1==1)                           % Keep looping till we get to the end of the block. This is when we have read in all the sweeps and read in 'END'
                    switch aTempString
                        
                        case 'ENABLED'                                                                                   % check if it is a simple sweep
                            NumberOfSweepSetsInFile=NumberOfSweepSetsInFile+1;                                                         % construct an sweep and pass it to the constructor for the VARSWPSWEEP
                            obj.ArrayOfSweepSets{NumberOfSweepSetsInFile}=SonnetSweepSetBlock(theFid);    	% construct the new sweep and store in the cell array, give the sweep the string so it knows what type of sweep it is
                        case 'SWPVARS'
                            if strcmp(fgetl(theFid), ' ON')
                                obj.SweepVariables = 1;
                            end
                        case 'END'   			% check if we reached the end of the block
                            fgetl(theFid);		% read the rest of the line we are on to get the theFid ready for the next block.
                            break;
                    end
                    
                    aTempString=fscanf(theFid,'%s',1);    % read the next sweep name from the file, if it is END then we are done
                    
                end
                
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % we come here when we didn't recieve a file ID as an argument
                % which means that we are going to create a default VARSWP block with
                % default values by calling the function's initialize method.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                initialize(obj);
                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initialize(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function initializes the VARSWP properties to some default
            %   values. This is called by the constructor and can
            %   be called by the user to reinitialize the object to
            %   default values.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            aBackup=warning();
            warning off all
            aProperties = properties(obj);
            for iCounter = 1:length(aProperties)
                obj.(aProperties{iCounter}) = [];
            end
            warning(aBackup);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetVariableSweepBlock();
            SonnetClone(obj,aNewObject);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function index = addSweepSet(obj, enabled)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function adds an empty sweep set to the block.
            % Enabled by default
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 1
                enabled = 'Y';
            end
            
            index = length(obj.ArrayOfSweepSets)+1;
            obj.ArrayOfSweepSets{index} = SonnetSweepSetBlock();
            if enabled == 'N' || enabled == 0
                obj.ArrayOfSweepSets{index}.disable()
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aSweepSet = getSweepSet(obj, index)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function gets a sweep set from the block. by index
            % in ArrayOfSweepSets
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(obj.ArrayOfSweepSets)
                obj.addSweepSet();
            end
            aSweepSet = obj.ArrayOfSweepSets{index};
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function writeObjectContents(obj, theFid, ~)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function writes the values from the object to a file.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Print the string signature to the file, more intensive but
            %less verbose
            fprintf(theFid, obj.stringSignature());
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aSignature=stringSignature(obj, ~)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function writes the values from the object to a string.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aSignature = 'VARSWP\n';
            
            if obj.SweepVariables
                aSignature = [aSignature 'SWPVARS ON\n'];
            end
            
            % Call the stringSignature function in each of the objects that we have in our cell array.
            for iCounter= 1:length(obj.ArrayOfSweepSets)
                aSignature = [aSignature obj.ArrayOfSweepSets{iCounter}.stringSignature()]; %#ok<AGROW>
            end
            
            aSignature = sprintf([aSignature 'END VARSWP\n']);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function activateVariableSweepParameter(obj,theVariableName,theSweepIndex)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   This function will set the parameter in
            %   use value for the specified parameter in the
            %   variable sweep.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin == 2
                theSweepIndex=1;
            end
            
            obj.ArrayOfSweeps{theSweepIndex}.activateVariableSweepParameter(theVariableName);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function deactivateVariableSweepParameter(obj,theVariableName,theSweepIndex)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   This function will set the parameter in
            %   use value for the specified parameter in the
            %   variable sweep.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin == 2
                theSweepIndex=1;
            end
            
            obj.ArrayOfSweeps{theSweepIndex}.deactivateVariableSweepParameter(theVariableName);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function changeVariableSweepParameterState(obj,theVariableName,theStatus,theSweepIndex)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   This function will modify the parameter in
            %   use value for the specified parameter in the
            %   variable sweep.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin == 3
                theSweepIndex=1;
            end
            
            obj.ArrayOfSweeps{theSweepIndex}.changeVariableSweepParameterState(theVariableName,theStatus);
            
        end
        
    end
end

