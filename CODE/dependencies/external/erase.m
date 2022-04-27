function s = erase(str, match)
%ERASE Remove content from text
%   NEWSTR = ERASE(STR,PAT) removes occurrences of PAT from STR. If
%   STR contains multiple occurrences of PAT, then ERASE removes all
%   occurrences that do not overlap.
%
%   STR must be text, which means it can be a string array, a character
%   vector, or a cell array of character vectors. PAT can be text or a
%   pattern array. PAT and STR need not have the same size. 
%   If PAT is a string array, a cell array, or a pattern array, then ERASE removes each
%   occurrence of all the elements of PAT from STR.
%
%   Examples:
%       STR = "The quick brown fox";
%       PAT = " quick";
%       erase(STR,PAT)            
%
%       % returns "The brown fox"
%
%       STR = 'Hello World';
%       PAT = 'Hello ';
%       erase(STR,PAT)            
%
%       % returns 'World'
%
%       STR = "The answer is 42!";
%       PAT = whitespacePattern + digitsPattern;
%       erase(STR,PAT)            
%
%       % returns "The answer is!"
%
%   See also STRREP, REGEXPREP, REPLACE, ERASEBETWEEN, PATTERN
 
%   Copyright 2016-2020 The MathWorks, Inc.
 
    narginchk(2, 2);
    
    %if ~isTextStrict(str)
    %    firstInput = getString(message('MATLAB:string:FirstInput'));
    %    error(message('MATLAB:string:MustBeCharCellArrayOrString', firstInput));
    %end
 
    %try
        s = string(str);
        s = s.erase(match);
        s = convertStringToOriginalTextType(s, str);
        
    %catch E
    %    throw(E)
    %end
end