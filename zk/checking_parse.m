function checking_parse(script, varargin)
p = inputParser;   % Create an instance of the inputParser class.

p.addRequired('script', @ischar);
p.addOptional('format', 'csv', @(x)ischar(x) && ~any(strcmpi(x,{'outputDir','maxHeight','maxWidth'})));
p.addParamValue('outputDir', pwd, @ischar);
p.addParamValue('maxHeight', [], @(x)x>0 && mod(x,1)==0);
p.addParamValue('maxWidth', [], @(x)x>0 && mod(x,1)==0);

p.parse(script, varargin{:})
disp(p.Results)


end