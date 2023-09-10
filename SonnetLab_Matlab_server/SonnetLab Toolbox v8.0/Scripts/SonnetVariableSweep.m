classdef SonnetVariableSweep < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class defines the sweep for a variable in the VARSWP block
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
        Options
        VariableName
        Arguments
        
    end
        
    methods
        function obj = SonnetVariableSweep(name, options, arguments)
            %SonnetVariableSweep(name, options, arguments) creates a sonnet variable
            %for a sweep set by parsing text from the file (VAR head not
            %included)
            if nargin > 0
                obj.VariableName = name;
                obj.Options = options;
                obj.Arguments = trystr2double(arguments);
                if iscolumn(obj.Arguments)
                    obj.Arguments = obj.Arguments.';
                end
            else
                initialize(obj)
            end
        end
        
        function initialize(obj)
            obj.Options = 'Y';
            obj.VariableName = 'Default';
            obj.Arguments = [0 10 1];
        end
        
        function enable(obj)
            obj.Options(1) = 'Y';
        end
        
        function disable(obj)
            obj.Options(1) = 'N';
        end
        
        function setOptions(obj, options)
            if checkOptions(options)
                obj.Options = options;
            else
                error("Invalid option: %s", options)
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetVariableSweep();
            SonnetClone(obj,aNewObject);
        end
        
        function aString = argumentString(obj) % string representation of arguments
            aString = regexprep(trynum2str(obj.Arguments), '( +)', ' '); % regex removes extra spaces
        end
        
        function aString = stringSignature(obj, theVersion)
            if nargin == 1
                theVersion = 17;
            end
            if theVersion >= 17
                aString = ['VAR ' obj.VariableName ' ' obj.Options ' ' obj.argumentString() '\n'];
            else
                aString = '';
            end
        end
        
        function writeObjectContents(obj, theFid, theVersion)
            fprintf(theFid, obj.stringSignature(theVersion));
        end
    end
end

function output = trystr2double(try_str)
if ~isa(try_str, 'double')
    output = str2double(try_str);
    if isnan(output)
        output = try_str;
    end
else
    output = try_str;
end
end

function output = trynum2str(try_num)
if isa(try_num, 'double')
    output = num2str(try_num);
else
    output = try_num;
end
end

function ok = checkOptions(option)
    valid = ["N" "Y" "YN" "YS" "YE" "YL" "YC"];
    ok = 0;
    for iCounter = 1:length(valid)
        if strcmp(option, valid(iCounter))
            ok = 1;
        end
    end
end