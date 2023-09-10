classdef SonnetGeometryVariableOption < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class defines the constraints for a variable (VCNSTR) in a Sonnet
    %   Geometry Block which defines geometry variables.
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
        Option
        VariableName
        Argument1
        Argument2
        Argument3
        
    end
        
    methods
        function obj = SonnetGeometryVariableOption(theFid)
            aTempString = split(strtrim(fgetl(theFid)));
            obj.Option = aTempString{1};
            obj.VariableName = aTempString{2};
            obj.Argument1 = trystr2double(aTempString{3});
            obj.Argument2 = trystr2double(aTempString{4});
            obj.Argument3 = trystr2double(aTempString{5});
        end
        
        function aString = stringSignature(obj, theVersion)
            if nargin == 1
                theVersion = 17;
            end
            if theVersion >= 17
                aString = ['VCNSTR ' obj.Option ' ' obj.VariableName ' ' trynum2str(obj.Argument1) ' ' trynum2str(obj.Argument2) ' ' trynum2str(obj.Argument3) '\n'];
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
output = str2double(try_str);
if isnan(output)
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