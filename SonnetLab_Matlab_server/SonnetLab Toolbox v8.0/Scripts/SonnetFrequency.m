classdef SonnetFrequency < handle
    %SonnetFrequency A subclass for all v17 Sonnet frequency sweeps
    %   All v17 sweeps can be enabled/disabled, and can be adaptive or not
    %   Enabled can be passed as 'Y', 'N', 1, or 0
    %   Adaptive can be passed as 'Y', 'N', 1, or 0, or 'X' if in an
    %      optimization, and therefore undefined
    %   Type is to make it easier to see the type in editor, and isn't used
    %   for anything
    properties (SetAccess = immutable)
        Type;
    end
    properties
        Enabled;
        AdaptiveEnabled;
    end
    
    methods
        function obj = SonnetFrequency(enabled,adaptive, type)
            %SonnetFrequency(enabled,adaptive, type) Construct an instance of this class
            if nargin < 3
                type = "undefined";
            end
            obj.Type = type;
            obj.Enabled = enableFix(enabled);
            if strcmp(adaptive, 'default')
                adaptive = 'N';
            end
            if ischar(adaptive)
                adaptive = upper(adaptive);
                if adaptive == 'Y' || adaptive == 'N' || adaptive == 'X'
                    obj.AdaptiveEnabled = adaptive;
                end
            elseif adaptive == 1
                obj.AdaptiveEnabled = 'Y';
            elseif adaptive == 0
                obj.AdaptiveEnabled = 'N';
            else
                error("adaptive must be 'Y', 'N', 'X', 1, or 0");
            end
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
        
        function aString = stringHead(obj)
            %stringHead Beginning of string output for this frequency sweep
            %   Detailed explanation goes here
            if obj.AdaptiveEnabled ~= 'X'
                aString = ['FREQ ' obj.Enabled ' A' obj.AdaptiveEnabled ' '];
            else
                aString = ['FREQ ' obj.Enabled ' '];
            end
        end
        
        function aString = stringTail(obj)
            aString = '\n';
        end
        
        function aString = toString(obj)
            aString = [obj.stringHead() obj.stringTail()];
        end
        
        function aString = stringSignature(obj)
            aString = [obj.stringHead() obj.stringTail()];
        end
        
        function writeObjectContents(obj, theFid, ~)
            fprintf(theFid, obj.toString());
        end
    end
end

