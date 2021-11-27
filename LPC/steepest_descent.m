function [a, e] = steepest_descent(s, p, n_steps)
%STEEPEST_DESCENT Computes the whitening filter coefficients via the
%steepest gradient descent method
%   The function parameters are
%       s: the short-time speech segment
%       p: the order of the filter
%       n_steps: the number of iterations
%   The function returns
%       a: the whitening filter coefficients 
%       e: the prediction error


% compute the autocorrelation function
[r, lags] = xcorr(s);


% compute the autocorrelation matrix R
R = r(lags >= 0);
R = R(1:p);
R = toeplitz(R);


% compute the autocorrelation vector r
r_v = r(lags >= 1);
r_v = r_v(1:p);

% determine the step-size parameter based on the stability constrain (take 
% the 0.9 of the maximum allowed value)   
lambda_max = max(eig(R));
mu = 0.9*(2/lambda_max);


% initialize prediction filter coeff to zero vector
a = zeros(p, 1);


for n = 1:n_steps
    % filter adaptation
    a = a + mu*(r_v - R*a); 
end

% define the whitening filter coefficients 
a = [1; -a];

% compute the prediction error
e = filter(a, 1, s);

end

% EOF