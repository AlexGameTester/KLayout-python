function out_enabled = enableFix(enabled)
%takes an enabled argument and fixes it to 'Y' or 'N'
enabled = upper(enabled);
if ischar(enabled) && (enabled == 'Y' || enabled == 'N')
    out_enabled = enabled;
elseif enabled == 1
    out_enabled = 'Y';
elseif enabled == 0
    out_enabled = 'N';
else
    error("enabled must be 'Y', 'N', 1, or 0");
end
end

