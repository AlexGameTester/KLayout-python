classdef SonnetGoalValue < handle
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Value
    end
    
    methods
        function obj = SonnetGoalValue(value)
            %SonnetGoalValue(value) Construct an instance of this class
            %   value: the value to compare
            if nargin > 0
                obj.Value = value;
            else
                initialize(obj)
            end
        end
        
        function initialize(obj)
            obj.Value = -30;
        end
        
        function aString = stringSignature(obj)
            %stringSignature Convert this class to a string as it would
            %appear in a file
            aString = ['VALUE=' num2str(obj.Value)];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetGoalValue();
            SonnetClone(obj,aNewObject);
        end
    end
end

