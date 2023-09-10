classdef SonnetFrequencySingle < SonnetFrequency
    %SonnetFrequencySingle Defines the linear sweep frequency
    %   Requires one frequency
    %   Printing methods already defined in parent class
    
    properties
        Frequency;
    end
    
    methods
        function obj = SonnetFrequencySingle(frequency, enabled, adaptive)
            %SonnetFrequencySingle Construct an instance of this class
            
            if nargin < 2
                enabled = 'Y';
            end
            if nargin < 3
                adaptive = 'Y';
            end
            obj = obj@SonnetFrequency(enabled, adaptive, 'Single');
            if nargin == 0
                initialize(obj)
                return
            end
            obj.Frequency = frequency;
        end
        
        function initialize(obj)
            %Set the object to its default values
            obj.Enabled = 'Y';
            obj.Frequency = 1;
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
            aNewObject=SonnetFrequencySingle();
            SonnetClone(obj,aNewObject);
        end
        
        function aString = stringTail(obj)
            %tailString End of string printout
            %   This comes after the head of the string, already defined in
            %   the parent class for this class
            aString = ['STEP ' num2str(obj.Frequency) '\n'];
        end
        
    end
end