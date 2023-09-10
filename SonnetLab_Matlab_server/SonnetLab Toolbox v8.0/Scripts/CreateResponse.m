function response = CreateResponse(type, arg1, arg2, arg3, arg4)
%CreateResponse Create a response to be used in part of an optimization goal
%   Create a response of type <type> to be used in a constructor like CreateGoal. 
%   type can be ('network, 'value', 'file') case insensitive
%   Types and arguments are as follows:
%
%       Network     format, parameter, ports
%       Value       value
%       File        format, parameter, ports, file
%
%
%   Example usage:
%       CreateResponse('Network', 'DB', 'S', [1 1])
%       CreateResponse('VALUE', 13)
%       CreateResponse('File', 'ANG', 'Y', [2 1], 'Relative_path/to_the/data_file')
%       CreateResponse('Value', 13)
%       CreateResponse('value', 13)

switch lower(type)
    case 'network'
        response = SonnetGoalNetwork(arg1, arg2, arg3); % pass to SonnetGoalNetwork constructor
    case 'value'
        response = SonnetGoalValue(arg1); % pass to SonnetGoalValue constructor
    case 'file'
        response = SonnetGoalFile(arg4, arg1, arg2, arg3); % pass to SonnetGoalFile constructor (arg4 is file, goes first in constructor)
end
end

