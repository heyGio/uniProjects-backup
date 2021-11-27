function [a, e] = levinson_durbin(s, p)
%LEVINSON Computes the whitening filter coefficients via the
%Levinson-Durbin recursion.
%   The function parameters are
%       s: the short-time speech segment
%       p: the order of the filter
%   The function returns
%       a: the whitening filter coefficients 
%       e: the prediction error


% compute the whitening filter coefficients 
[r, lags] = xcorr(s);
r = r(lags>=0);
r = r(1:p);
a = levinson(r, p);
a = a';


% compute the prediction error
e = filter(a, 1, s);


end

% EOF