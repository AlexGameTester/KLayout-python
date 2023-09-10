classdef SonnetGoalFile < handle
    %SonnetGoalValue Summary of this class goes here
    %   Detailed explanation goes here
    
properties
        File;
        Format;
        Parameter;
        Ports;
    end
    
    methods
        function obj = SonnetGoalFile(file, format, parameter, ports)
            %SonnetGoalFile Construct an instance of this class
            %   file: the filename of the data to compare
            %   format: the format property
            %       (DB, ANG, MAG, RE, IM)
            %   parameter: the parameter to use
            %       (S, Z, Y, ...)
            %   ports: an array of input and output port
            %       ([1 1], [1 2], ...)
            if nargin > 0
                obj.File = file;
                obj.Format = format;
                obj.Parameter = parameter;
                obj.Ports = ports;
            else
                initialize(obj)
            end
        end
        
        function initialize(obj)
            obj.File = 'default';
            obj.Format = 'DB';
            obj.Parameter = 'S';
            obj.Ports = [1 1];
        end
        
        function aString = stringSignature(obj)
            %stringSignature Convert this class to a string as it would
            %appear in a file
            aString = ['FILE ' obj.File ' ' obj.Format '[' obj.Parameter num2str(obj.Ports(1)) num2str(obj.Ports(2)) ']'];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function aNewObject=clone(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function builds a deep copy of this object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aNewObject=SonnetGoalFile();
            SonnetClone(obj,aNewObject);
        end
    end
end