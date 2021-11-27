function [a, e] = lms_algorithm(s, p, mu)
%LMS Computes the whitening filter coefficients via the LMS algorithm.
%   The function parameters are
%       s: the short-time speech segment
%       p: the order of the filter
%       mu: the step-size parameter
%   The function returns
%       a: the whitening filter coefficients 
%       e: the prediction error

M = length(s);
e = zeros(M,1);
y = zeros(M,1);
a = zeros(p,1);
   
for m = 1:M
    
    % input vector
    u = s( m-1 : -1 : max(1,m-p));
    
    % zero pad if needed
    if length(u) < p
        u = [u; zeros(p-length(u),1)];
    end
        
    % filter output
    y(m) = conj(a') * u;

    
    % desired response
    d = s(m);

    
    % error signal
    e(m) = d - y(m);
    
    % stochastic gradient
    grad = mu * u * conj(e(m));
    
    % filter adaptation
    a = a + grad;
    
end

% define the whitening filter coefficients
a = [1; -a];


end

% EOF