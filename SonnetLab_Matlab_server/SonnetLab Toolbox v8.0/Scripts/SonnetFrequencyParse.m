function sweep = SonnetFrequencyParse(text, enabled, adaptive, arg1, arg2, arg3)
%sonnetFrequncyParse Return a sonnet frequency by parsing text
%   Will ignore FREQ and a space character at beginning
%   can also take a list of arguments
%   Types and arguments are as follows:
            %
            %       SWEEP       StartFrequency,EndFrequency,StepFrequency
            %       ABS_ENTRY    StartFrequency,EndFrequency, FrequencyCount*
            %       DC          Mode**,Frequency**
            %       ESWEEP      StartFrequency,EndFrequency,AnalysisFrequencies
            %       LSWEEP      StartFrequency,EndFrequency,AnalysisFrequencies
            %       SIMPLE      StartFrequency,EndFrequency,StepFrequency
            %       STEP        StepFrequency
            %
            %       *   For ABS_ENTRY:  FrequencyCount defaults to 300
            %       **  For a DC sweep: mode is either 'AUTO' for automatic or 'MAN' for manual.
            %       **  For a DC sweep: when mode is 'AUTO' the frequency does not need to
            %                           be supplied. The frequency is required when the DC
            %                           mode is manual.
            
if nargin == 1      %do this if parsing text
    splitText = split(strtrim(text));
    if strcmp(splitText(1), 'FREQ')
        splitText = splitText(2:end);
    end
    enabled = splitText{1};
    adaptive = splitText{2}(2);
    if ~strcmp(adaptive, 'Y') && ~strcmp(adaptive, 'N')
        adaptive = 'X';
    end
    if length(splitText) < 4
        sweep = SonnetFrequency(enabled, adaptive);
        return
    end
    args = [str2double(splitText(4:end)); zeros(4, 1)]; %zeros buffer for nonexistant values
    aTypeIndex = 3;
    while aTypeIndex >= 0
        switch splitText{aTypeIndex}
            case 'SWEEP'
                sweep = SonnetFrequencyLinear(args(1), args(2), args(3), enabled, adaptive);
                break
            case 'LSWEEP'
                sweep = SonnetFrequencyLinear(args(1), args(2), args(3), enabled, adaptive, 1);
                break
            case 'ESWEEP'
                sweep = SonnetFrequencyExponential(args(1), args(2), args(3), enabled, adaptive);
                break
            case 'ABS_ENTRY'
                sweep = SonnetFrequencyAdaptive(args(1), args(2), args(3), args(4), enabled, adaptive);
                break
            case 'LIST'
                sweep = SonnetFrequencyList(args(1:end-4), enabled, adaptive);
                break
            case 'STEP'
                sweep = SonnetFrequencySingle(args(1), enabled, adaptive);
                break
            case 'DC_FREQ'
                if length(args)>1
                    sweep = SonnetFrequencyDc('MAN', args(2), enabled, adaptive);
                    break
                else
                    sweep = SonnetFrequencyDc('AUTO', enabled, adaptive);
                    break
                end
            otherwise
                aTypeIndex = aTypeIndex - 1; %if we didn't find the type, we could be too far
        end
    end
else        %do this if arguments recieved
    switch nargin
        case 4
            sweep = frequencySweep(text, enabled, adaptive, arg1);
        case 5
            sweep = frequencySweep(text, enabled, adaptive, arg1, arg2);
        case 6
            sweep = frequencySweep(text, enabled, adaptive, arg1, arg2, arg3);
    end
end
end

function sweep = frequencySweep(type, enabled, adaptive, arg1, arg2, arg3)
%create and return sweeps 
if nargin == 4
    arg2 = 0;
    arg3 = 0;
elseif nargin == 5
    arg3 = 0;
end
switch upper(type) % Correct common aliases
    case 'ABS'
        type = 'ABS_ENTRY';
    case 'ADAPTIVE'
        type = 'ABS_ENTRY';
    case 'EXPONENTIAL'
        type = 'ABS_ENTRY';
    case 'DC'
        type = 'DC_FREQ';
    case 'LINEAR'
        type = 'LSWEEP';
    case 'SINGLE'
        type = 'LSWEEP';
end
switch type     % construct frequency sweeps based on value of type
    case 'SWEEP'
        sweep = SonnetFrequencyLinear(arg1, arg2, arg3, enabled, adaptive);
    case 'LSWEEP'
        sweep = SonnetFrequencyLinear(arg1, arg2, arg3, enabled, adaptive, 1);
    case 'ESWEEP'
        sweep = SonnetFrequencyExponential(arg1, arg2, arg3, enabled, adaptive);
    case 'ABS_ENTRY'
        if arg3 == 0
            arg3 = 300;
        end
        sweep = SonnetFrequencyAdaptive(arg1, arg2, -1, arg3, enabled, adaptive);
    case 'LIST'
        sweep = SonnetFrequencyList(arg1, enabled, adaptive);
    case 'STEP'
        sweep = SonnetFrequencySingle(arg1, enabled, adaptive);
    case 'DC_FREQ'
        if strcmp(arg1, 'MAN')
            sweep = SonnetFrequencyDc('MAN', arg2, enabled, adaptive);
        else
            sweep = SonnetFrequencyDc('AUTO', enabled, adaptive);
        end
    otherwise
        error("type invalid for frequency sweep")
end
end
