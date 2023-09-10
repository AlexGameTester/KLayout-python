function variable_sweep = CreateVariableSweep(variable_name, type, enabled, arguments)
%CreateVariableSweep(variable_name, type, enabled, arguments) creates a
%variable sweep to be added to a sweep set
%   variable_name: Names for variables contain only letters, digits and underscores and must start with a letter or underscore. Also predefined names of constants and functions may not be used.
%   type: Type of sweep ('single', 'linear', 'number', 'exponential', 'list', 'corner', 'sensitivity')
%       'number' is linear for number of points
%   enabled: 'Y', 'N', 1, or 0
%   arguments: arguments for the sweep type, typically in the order 
%       (min, max, number)
%       linear: [min, max, step]
%       exponential:    [min, max, #points]
%       list:   1d array of arguments
%       single: value
%       corner: [min, max, nominal]
%       sensitivity: [min, max, nominal]
enabled = enableFix(enabled);
switch lower(type)
    case 'linear'
        type = '';
    case 'lin'
        type = '';
    case 'fixed'
        type = '';
    case 'f'
        type = '';
    case 'none'
        type = '';
    case 'single'
        type = '';
    case 'exp'
        type = 'E';
    case 'exponential'
        type = 'E';
    case 'e'
        type = 'E';
    case 'list'
        type = 'L';
    case 'ls'
        type = 'L';
    case 'corner'
        type = 'C';
    case 'sensitivity'
        type = 'S';
    case 'step'
        type = 'N';
    case 'number'
        type = 'N';
    case 'n'
        type = 'N';

end
variable_sweep = SonnetVariableSweep(variable_name, [enabled type], arguments);
end

