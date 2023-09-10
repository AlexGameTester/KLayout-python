function sweep = CreateSweep(type, enabled, adaptive, arg1, arg2, arg3)
%CreateSweep Creates a sonnet frequency sweep that can be added to
%sweep sets or optimizations
%   CreateSweep(type, enabled, adaptive, arg1, arg2, arg3) creates a frequency sweep
%   CreateSweep requires a string specifying the type
%   of frequency sweep to be added to the project and
%   all of the arguments necessary in order to construct
%   the sweep.
%
%   ONLY VALID v17+
%
%   Types and arguments are as follows:
%       
%       enabled and adaptive should each be either 'Y' or 'N' (default enabled = 'Y' adaptive = 'N')
%       type should be one of the following values (in single quotes) ex: 'SWEEP'
%       arg1, arg2, arg3 are what follow in the comma separated list
%
%       SWEEP       StartFrequency,EndFrequency,StepFrequency
%       ABS_ENTRY   StartFrequency,EndFrequency, FrequencyCount*
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

switch upper(type) % fixing common unambiguous aliases
    case 'EXPONENTIAL'
        type = 'ESWEEP';
    case 'SINGLE'
        type = 'STEP';
    case 'ABS'
        type = 'ABS_ENTRY';
    case 'ADAPTIVE'
        type = 'ABS_ENTRY';
    case 'DC_FREQ'
        type = 'DC';
end

switch upper(type) % create sweep for given type, intended to be used for sweep set
    case 'SWEEP'
        sweep = SonnetFrequencyLinear(arg1, arg2, arg3, enabled, adaptive);
    case 'ABS_ENTRY'
        if nargin < 6
            arg3 = 300;
        end
        sweep = SonnetFrequencyAdaptive(arg1, arg2, -1, arg3, enabled, adaptive);    
    case 'DC'
        sweep = SonnetFrequencyDc(arg1, arg2, enabled, adaptive);
    case 'ESWEEP'
        sweep = SonnetFrequencyExponential(arg1, arg2, arg3, enabled, adaptive);
    case 'LIST'
        sweep = SonnetFrequencyList(arg1, enabled, adaptive);
    case 'STEP'
        sweep = SonnetFrequencySingle(arg1, enabled, adaptive);
    otherwise
        error('type invalid. Must be (SWEEP, ABS_ENTRY, DC, ESWEEP, LIST, STEP')
end
end