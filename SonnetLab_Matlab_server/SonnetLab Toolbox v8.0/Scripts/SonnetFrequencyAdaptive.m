classdef SonnetFrequencyAdaptive < SonnetFrequency
    %SonnetFrequencyAdaptive Defines the linear sweep frequency
    %   Requires start, end, and step frequencies
    %   Printing methods already defined in parent class
    
    properties
        StartFrequency;
        EndFrequency;
        ThirdValue;
        FrequencyCount;
    end
    
    methods
        function obj = SonnetFrequencyAdaptive(start_frequency, end_frequency, third, count, enabled, adaptive)
            %SonnetFrequencyAdaptive Construct an instance of this class
            %start_frequency, end_frequency, count should all be numbers
            %enabled and adaptive can be 'Y' 'N' 1 or 0
            if nargin < 5
                enabled = 'Y';
            end
            if nargin < 6 || strcmp(adaptive, 'default')
                adaptive = 'Y';
            end
            obj = obj@SonnetFrequency(enabled, adaptive, 'Adaptive');
            if nargin == 0
                initialize(obj)
                return
            end
            obj.StartFrequency = start_frequency;
            obj.EndFrequency = end_frequency;
            obj.ThirdValue = third;
            obj.FrequencyCount = count;
        end
        
        function initialize(obj)
            obj.Enabled = 'Y';
            obj.StartFrequency = 1;
            obj.EndFrequency = 10;
            obj.ThirdValue = -1;
            obj.FrequencyCount = 300;  
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetFrequencyAdaptive();
            SonnetClone(obj,aNewObject);
        end
        
        function set(obj, start_frequency, end_frequency)
            obj.StartFrequency = start_frequency;
            obj.EndFrequency = end_frequency;
        end
        
        function setStartFrequency(obj, frequency)
            obj.StartFrequency = frequency;
        end
        
        function setEndFrequency(obj, frequency)
            obj.EndFrequency = frequency;
        end
        
        function aString = stringTail(obj)
            %tailString End of string printout
            %   This comes after the head of the string, already defined in
            %   the parent class for this class
            aString = ['ABS_ENTRY ' num2str(obj.StartFrequency) ' ' num2str(obj.EndFrequency) ' ' num2str(obj.ThirdValue) ' ' num2str(obj.FrequencyCount) '\n'];
        end
    end
end

