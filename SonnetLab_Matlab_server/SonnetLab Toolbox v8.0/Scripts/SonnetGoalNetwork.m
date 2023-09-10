classdef SonnetGoalNetwork < handle
    %SonnetGoalNetwork The network response goal
    %   Intended to be used as part of a SonnetGoal for optimizations
    %   Create with CreateResponse(...)
    
    properties
        Format;
        Parameter;
        Ports;
    end
    
    methods
        function obj = SonnetGoalNetwork(format, parameter, ports)
            %SonnetGoalNetwork Construct an instance of this class
            %   format: the format property
            %       (DB, ANG, MAG, RE, IM)
            %   parameter: the parameter to use
            %       (S, Z, Y, ...)
            %   ports: an array of input and output port
            %       ([1 1], [1 2], ...)
            if nargin > 0
                obj.Format = format;
                obj.Parameter = parameter;
                obj.Ports = ports;
            else
                initialize(obj)
            end
        end
        
        function initialize(obj)
            obj.Format = 'DB';
            obj.Parameter = 'S';
            obj.Ports = [1 1];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetGoalNetwork();
            SonnetClone(obj,aNewObject);
        end
        
        function aString = stringSignature(obj)
            %stringSignature Convert this class to a string as it would
            %appear in a file
            aString = ['NET=GEO ' obj.Format '[' obj.Parameter num2str(obj.Ports(1)) num2str(obj.Ports(2)) ']'];
        end
    end
end

