function aGoal = SonnetGoalParse(text)
%SonnetGoalParse Summary of this function goes here
%   Parses the goal based on the order elements appear, separated by spaces

strrep(text, 'GOAL', '');        %remove GOAL if there
splitText = split(strtrim(text));

enabled = splitText{1};

[format, parameter, ports] = extractFormatParameterPorts(splitText{3});
goal1 = SonnetGoalNetwork(format, parameter, ports);
comparator = splitText{4};
splitText = splitText(5:end);
switch strtok(splitText{1}, '=')
    case 'NET'
        [format, parameter, ports] = extractFormatParameterPorts(splitText{2});
        goal2 = SonnetGoalNetwork(format, parameter, ports);
    case 'VALUE'
        value = sscanf(splitText{1}, 'VALUE=%d');
        goal2 = SonnetGoalValue(value);
    case 'FILE'
        file = regexp(text, '(?:FILE )(.+)(?: [A-Z]{2,}\[)', 'tokens');
        file = file{1}{1};
        [format, parameter, ports] = extractFormatParameterPorts(splitText{end-1});
        goal2 = SonnetGoalFile(file, format, parameter, ports);
end
weight = str2double(splitText(end));
aGoal = SonnetGoal(enabled, goal1, comparator, goal2, weight);
end

function [format, remains, ports] = extractFormatParameterPorts(text)
[format, remains] = strtok(text, '[');
ports = [str2double(remains(end-2)) str2double(remains(end-1))];
remains = remains(2:end-3);
end
