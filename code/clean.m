function [varargout] = clean(varargin)
%clean Clean stray data wiht high deviation from mean
%   Data with high deviation from the mean value will get "cleaned"
input = varargin{1};
k = 2;
if ~isempty(varargin{2})
    k = varargin{2};
end
value = input.Data;
std_div = sqrt(var(value));
mean_v = mean(value);

place = abs(input.Data-mean_v)<= k*std_div;
input_clean = timeseries(value(place),input.Time(place));

if nargout==1
    varargout = input_clean;
else
    varargout{1} = input_clean;
    varargout{2} = place;
end


end

