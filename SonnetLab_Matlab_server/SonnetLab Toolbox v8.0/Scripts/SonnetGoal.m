classdef SonnetGoal < handle
    %SonnetGoal Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Enabled;
        Response1;
        Comparator;
        Response2;
        Weight;
    end
    
    methods
        function obj = SonnetGoal(enabled, response1, comparator, response2, weight)
            %SonnetGoal Construct an instance of this class
            %   It is recommended to use CreateGoal(...) to create an
            %   instance of this class
            % Arguments as follows:
            %   1) enabled      - Whether the goal will be used in simulation
            %   2) response1    - The response to optimize
            %   3) comparator   - The comparison string ('<', '=', '>')
            %   4) response2    - The response to compare
            %   5) weight       - The weight of this goal
            if nargin > 0
                if ischar(enabled) && (enabled == 'Y' || enabled == 'N')
                    obj.Enabled = enabled;
                elseif enabled == 1
                    obj.Enabled = 'Y';
                elseif enabled == 0
                    obj.Enabled = 'N';
                else
                    error("enabled must be 'Y', 'N', 1, or 0");
                end
                obj.Response1 = response1;
                obj.Comparator = strtrim(comparator);
                obj.Response2 = response2;
                obj.Weight = weight;
            else
                initialize(obj);
            end
        end
        
        function initialize(obj)
            obj.Enabled = 'Y';
            obj.Response1 = SonnetGoalNetwork();
            obj.Comparator = '=';
            obj.Response2 = SonnetGoalValue();
            obj.Weight = 1;            
        end
        
        function aString = stringSignature(obj)
            %SonnetGoal Summary of this method goes here
            %   Detailed explanation goes here
            aString = ['GOAL ' obj.Enabled ' ' obj.Response1.stringSignature() ' ' obj.Comparator ' ' obj.Response2.stringSignature() ' ' num2str(obj.Weight) '\n'];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetGoal();
            SonnetClone(obj,aNewObject);
        end
    end
end

