classdef SonnetFrequencyLinear < SonnetFrequency
    %SonnetFrequencyLinear Defines the linear sweep frequency
    %   Requires start, end, and step frequencies
    %   Defines a linear frequency sweep
    %       (Equal spacing between frequencies)
    
    properties
        StartFrequency
        EndFrequency
        Step
        IsFrequencyCount
    end
    
    methods
        function obj = SonnetFrequencyLinear(start_frequency, end_frequency, step, enabled, adaptive, is_frequency_count)
            %SonnetFrequencyLinear Construct an instance of this class
            %start_frequency, end_frequency, step should all be numbers
            %enabled and adaptive can be 'Y' 'N' 1 or 0
            %is_frequency_count determines whether step is the space, logical 1 or 0
            %   between frequencies or the number of frequencies
            if nargin < 4
                enabled = 'Y';
            end
            if nargin < 5
                adaptive = 'Y';
            end
            if nargin < 6
                is_frequency_count = 0;
            end
            obj = obj@SonnetFrequency(enabled, adaptive, 'Linear');
            if nargin == 0
                initialize(obj)
                return
            end
            obj.StartFrequency = start_frequency;
            obj.EndFrequency = end_frequency;
            obj.Step = step;
            obj.IsFrequencyCount = is_frequency_count;
        end
        
        function initialize(obj)
            %Set the object to its default values
            obj.Enabled = 'Y';
            obj.StartFrequency = 1;
            obj.EndFrequency = 10;
            obj.Step = 1;  
        end
        
        function set_with_step(obj, start_frequency, end_frequency, step)
            %Set the parameters for the sweep
            obj.StartFrequency = start_frequency;
            obj.EndFrequency = end_frequency;
            obj.Step = step;
            obj.IsFrequencyCount = 0;
        end
        function set_with_frequency_count(obj, start_frequency, end_frequency, frequency_count)
            %Set the parameters for the sweep
            obj.StartFrequency = start_frequency;
            obj.EndFrequency = end_frequency;
            obj.Step = frequency_count;
            obj.isFrequencyCount = 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetFrequencyLinear();
            SonnetClone(obj,aNewObject);
        end
        
        function aString = stringTail(obj)
            %tailString End of string printout
            %   This comes after the head of the string, already defined in
            %   the parent class for this class
            beginning = 'SWEEP ';
            if obj.IsFrequencyCount
                beginning = 'LSWEEP ';
            end
            aString = [beginning num2str(obj.StartFrequency) ' ' num2str(obj.EndFrequency) ' ' num2str(obj.Step) '\n'];
        end
        
    end
end

