function p = hnormpdf(x,sigma)
%This is a pdf for half normal distribution 
%   x is the argument, it can be a vector, must be non-negative
%   sigma is standard deviation, it can be a vector, must be positive
% WARNING: I do not check here if x>=0 and if sigma>0 becasue in
% my code I do this before executing this function. 
% Author:   Kamil Makiela
%           Department of econometrics and operations research 
%           Cracow University of Economics
p = sqrt(2/pi).*exp(-(x.^2)./(2.*sigma.^2))./sigma;

end

