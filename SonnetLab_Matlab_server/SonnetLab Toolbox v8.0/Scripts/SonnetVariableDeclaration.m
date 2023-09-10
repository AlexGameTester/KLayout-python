function variable_sweep = SonnetVariableDeclaration(theText)
%SonnetVariableDeclaration(theText) creates a sonnet variable sweep
%for a sweep set by parsing text from the file (VAR head not included)

aTempString = split(strtrim(theText));
variable_sweep = SonnetVariableSweep(aTempString{1}, aTempString{2}, aTempString(3:end));
end