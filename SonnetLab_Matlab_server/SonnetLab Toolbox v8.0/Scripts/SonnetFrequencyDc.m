classdef SonnetFrequencyDc < SonnetFrequency
    %SonnetFrequencyDc Defines the linear sweep frequency
    %   Requires one frequency, and auto identifier
    %   Printing methods already defined in parent class
    
    properties
        Auto;
        Frequency;
    end
    
    methods
        function obj = SonnetFrequencyDc(auto, frequency, enabled, adaptive)
            %SonnetFrequencyDc Construct an instance of this class
            %valid values of auto are 'AUTO' or 'MAN'
            if nargin < 3
                enabled = 'Y';
            end
            if nargin < 4
                adaptive = 'Y';
            end
            obj = obj@SonnetFrequency(enabled, adaptive, 'Dc Frequency');
            if nargin == 0
                initialize(obj)
                return
            end
            if strcmp(auto, 'AUTO')
                obj.Auto = 'AUTO';
            elseif strcmp(auto, 'MAN') || strcmp(auto, 'MANUAL')
                obj.Auto = 'MAN';
                obj.Frequency = frequency;
            else
                error("auto should be 'AUTO' or 'MAN'")
            end
        end
        
        function initialize(obj)
            %Set the object to its default values
            obj.Enabled = 'Y';
            obj.Auto = 'AUTO';  
        end
        
        function manual(obj)
            %Set to manual mode
            obj.Auto = 'MAN';
        end
        
        function auto(obj)
            %set to auto mode
            obj.Auto = 'AUTO';
        end
        
        function toggle_auto(obj)
            %Toggle auto mode
            if strcmp(obj.Auto, 'AUTO')
                obj.Auto = 'MAN';
            else
                obj.Auto = 'AUTO';
            end
        end
        
        function set.Frequency(obj, frequency)
            obj.Frequency = frequency;
        end
        
        function frequency = get.Frequency(obj)
            frequency = obj.Frequency;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetFrequencyDc();
            SonnetClone(obj,aNewObject);
        end
        
        function aString = stringTail(obj)
            %tailString End of string printout
            %   This comes after the head of the string, already defined in
            %   the parent class for this class
            if strcmp(obj.Auto, 'AUTO')
                aString = 'DC_FREQ AUTO\n';
            else
                aString = ['DC_FREQ MAN ' num2str(obj.Frequency) '\n'];
            end
        end
    end
end