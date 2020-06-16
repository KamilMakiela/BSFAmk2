function L = lgL_re(s_v, s_u, b, X, y, n ,t)
% likelihood for normal-half-normal with persitent efficiency
% based on Pitt and Lee (1981, p. 60 and 61)
sigs2 = s_v^2 + t*s_u^2;
s2v = s_v^2;
s2u = s_u^2;
e = y - X*b;
ee = reshape(e,t,n);
e_sr = mean(ee,1)';

tt = ones(t,t);         % ll' from p. 60 
%A = eye(t) - tt.*(s2u/sigs2);
A = (eye(t) - tt./t + (tt./t).*(s2v/sigs2));

%arg1 = n*log(2) - 0.5*n*t*log(2*pi) - 0.5*n*(t-1)*log(s2v) - 0.5*n*log(sigs2);
arg1  =          - 0.5*n*t*log(2*pi) - 0.5*n*(t-1)*log(s2v) - 0.5*n*log(sigs2) ;
%arg1 = -0.5*n*t*log(2*pi) - 0.5*n*log(sigs2*s2v^(t-1));

%arg2 = - e'*kron(eye(n),A)*e./(2*s2v);
arg2 =  - e'*kron(eye(n),A)*e./(2*s2v);
%arg3 = sum(log(1-normcdf(e_sr.*(t*s_u/(s_v*sqrt(sigs2))))));

L = arg1 + arg2;

end

