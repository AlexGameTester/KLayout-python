classdef SonnetFrequencyList < SonnetFrequency
    %SonnetFrequencyList Defines the linear sweep frequency
    %   Requires start, end, and step frequencies
    %   Printing methods already defined in parent class
    
    properties
        List;
    end
    
    methods
        function obj = SonnetFrequencyList(list, enabled, adaptive)
            %SonnetFrequencyList Construct an instance of this class
            %list is a 1-d array of numbers
            %   (can single be column or row)
            if nargin < 2
                enabled = 'Y';
            end
            if nargin < 3
                adaptive = 'Y';
            end
            obj = obj@SonnetFrequency(enabled, adaptive, 'List');
            if nargin == 0
                default(obj)
                return
            end
            shape = size(list);
            if shape(1) == 1
                obj.List = list;
            else
                obj.List = transpose(list);
            end
            shape = size(obj.List);
            if shape(1) ~= 1
                error("list must be a vector")
            end
        end
        
        function initialize(obj)
            %Set the object to its default values
            obj.List = [];
        end
        
        function aList = getList(obj)
            %Set the list for the sweep
            aList = obj.List;
        end
        
        function setList(obj, list)
            %Set the list for the sweep
            obj.List = list;
        end
        
        function aString = stringTail(obj)
            %tailString End of string printout
            %   This comes after the head of the string, already defined in
            %   the parent class for this class
            aString = ['LIST ' num2str(obj.List) '\n'];
        end
        
        function returnObj = defualt(obj)
            obj.List = [];
            returnObj = obj;
        end
    end
end