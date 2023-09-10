classdef SonnetFrequencyExponential < SonnetFrequency
    %SonnetFrequencyExponential Defines the linear sweep frequency
    %   Requires start, end, and count of frequencies
    
    properties
        StartFrequency;
        EndFrequency;
        FrequencyCount;
    end
    
    methods
        function obj = SonnetFrequencyExponential(start_frequency, end_frequency, frequency_count, enabled, adaptive)
            %SonnetFrequencyExponential Construct an instance of this class
            %start_frequency, end_frequency, frequency_count should all be numbers
            %enabled and adaptive can be 'Y' 'N' 1 or 0
            if nargin < 4
                enabled = 'Y';
            end
            if nargin < 5
                adaptive = 'Y';
            end
            obj = obj@SonnetFrequency(enabled, adaptive, 'Exponential');
            if nargin == 0
                initialize(obj)
                return
            end
            obj.StartFrequency = start_frequency;
            obj.EndFrequency = end_frequency;
            obj.FrequencyCount = frequency_count;
        end
        
        function initialize(obj)
            %Set the object to its default values
            obj.Enabled = 'Y';
            obj.StartFrequency = 1;
            obj.EndFrequency = 10;
            obj.FrequencyCount = 300;  
        end
        
        function set(obj, start_frequency, end_frequency, number_of_points)
            %Set the parameters for the sweep
            obj.StartFrequency = start_frequency;
            obj.EndFrequency = end_frequency;
            obj.FrequencyCount = number_of_points;
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetFrequencyExponential();
            SonnetClone(obj,aNewObject);
        end
        
        
        function aString = stringTail(obj)
            %tailString End of string printout
            %   This comes after the head of the string, already defined in
            %   the parent class for this class
            aString = ['ESWEEP ' num2str(obj.StartFrequency) ' ' num2str(obj.EndFrequency) ' ' num2str(obj.FrequencyCount) '\n'];
        end
    end
end

