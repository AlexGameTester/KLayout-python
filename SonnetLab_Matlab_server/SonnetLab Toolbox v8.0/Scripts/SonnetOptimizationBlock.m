classdef SonnetOptimizationBlock  < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class defines the OPT portion of a SONNET project file.
    % This class is a container for the optimization information that
    % is obtained from the SONNET project file.
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
        MaxIterations
        GoalSets                    % This is an array of goal sets
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = SonnetOptimizationBlock(theFid)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % The constructor for OPT.
            %     the OPT will be passed the file ID from the
            %     SONNET project constructor.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin == 1                  % This checks if we got the file ID as an argument
                
                initialize(obj);            % Initialize the values of the properties using the initializer function
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Read a string from the file.
                % This String drives a switch
                % statement to determine what
                % values are going to be changed.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                aTempString=fscanf(theFid,' %s',1);
                
                NumberOfSets = 0;
                
                while (1==1)                % Loop forever till we get to the end of the paramter list for this sweep, there can be an undefined number of parameters
                    
                    switch aTempString
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % We want to read in a the MaxIterations value which will be
                        % after the keyword 'MAX'; we want to read in (and ignore) the
                        % term 'VARS'.  Then we want to read in an undisclosed number
                        % of optimization variables. Then we will read in a Frequency
                        % which will mean that we have completed reading in the
                        % variables.  Then we will read in an optimization line that
                        % begins with 'NET=GEO'.
                        %
                        % We can run this all in a case statement by reading in a
                        % String, then checking if it is MAX or VARS. If it is then
                        % Do the appropriate action.  Otherwise check if it is the
                        % same as a sweep name.  If it isn't a sweep then it is an
                        % optimization variable.  If it is a sweep that means we
                        % read in all of the optmization variables.  We can construct
                        % the appropriate sweep object based on the keyword we read.
                        % We then can read in the optimization parameter.
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        case 'MAX'                                                    % If we are reading in the MAX number of iterations
                            obj.MaxIterations=fscanf(theFid,' %d', 1);
                            
                        case 'VARS'                                                   % If we are reading in the VARS line then just ignore it, we dont need to save it.
                        
                        case 'ENABLED'
                            NumberOfSets = NumberOfSets + 1;
                            obj.GoalSets{NumberOfSets} = SonnetGoalSetBlock(theFid);
                            
                        case 'END'   % Check if it is END indicating we are done reading OPT
                            % Some projects have an END for vars and then a END OPT line
                            % We need to be sure we are done with the OPT block when we
                            % finish this constructor
                            aTempString=fgetl(theFid);
                            if ~isempty(aTempString)
                                break;
                            end

                        case 'VAR'
                            aTempString=fscanf(theFid,' %s',1);
                            NumberOfParamters=NumberOfParamters+1;                                            % Increment the parameter counter by one
                            obj.VarsArray{NumberOfParamters}=SonnetOptimizationVariable(theFid,aTempString);  % construct the new parameter and store in the cell array, give the parameter the string so it knows its name
                            
                        otherwise                                                                             % Otherwise it is a parameter and we should make a parameter object
                            NumberOfParamters=NumberOfParamters+1;                                            % Increment the parameter counter by one
                            obj.VarsArray{NumberOfParamters}=SonnetOptimizationVariable(theFid,aTempString);  % construct the new parameter and store in the cell array, give the parameter the string so it knows its name
                            
                    end
                    
                    aTempString=fscanf(theFid,' %s',1);   % read the next sweep name from the file, if it is END then we are done
                    
                end
                
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % we come here when we didn't recieve a file ID as an argument
                % which means that we are going to create a default OPT block with
                % default values by calling the function's initialize method.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                initialize(obj);
                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initialize(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function initializes the OPT properties to some default
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
            
            obj.MaxIterations=100;

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function index = addGoalSet(obj, enabled)
            index = length(obj.GoalSets)+1;
            enabled = enableFix(enabled);
            obj.GoalSets{index} = SonnetGoalSetBlock();
            obj.GoalSets{index}.Enabled = enabled;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetOptimizationBlock();
            SonnetClone(obj,aNewObject);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function writeObjectContents(obj, theFid, theVersion)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function writes the values from the object to a file.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf(theFid,'OPT\n');
            
            fprintf(theFid,'MAX %d\n',obj.MaxIterations);
            
            % We want to loop and print out the values for all the sweeps and optimizations
            for iCounter= 1:length(obj.GoalSets)
                
                % Call the writeguts method for the frequency type
                obj.GoalSets{iCounter}.writeObjectContents(theFid,theVersion);
                
            end
            
            fprintf(theFid,'END OPT\n');
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aSignature=stringSignature(obj, ~)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function writes the values from the object to a string.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            aSignature = sprintf('OPT\nMAX %d\n', obj.MaxIterations);
                        
            % We want to loop and print out the values for all the sweeps and optimizations
            for iCounter= 1:length(obj.GoalSets)
                
                % Call the writeguts method for the frequency type
                aSignature = [aSignature obj.GoalSets{iCounter}.stringSignature()];
                
            end
            
            aSignature = [aSignature 'END OPT\n'];
        end
        
    end
end
