function L = lgL_nexp(s_v, s_u, b, X, y, n ,t)
% likelihood for normal-half-normal with persitent efficiency
% based on Pitt and Lee (1981, p. 60 and 61)
sigs2 = s_v^2 + t*s_u^2;
s2v = s_v^2;
s2u = s_u^2;
e = y - X*b;
%e_sr = y_sr - X_sr*b;
%tt = ones(t,t);         % ll' from p. 60 
ee = reshape(e,t,n);

e_sr = mean(ee,1)';
%disp(size(e_sr));
sse = diag(ee'*ee); % suma kwadratow reszt dla danego obiektu
%disp(size(sse));
%A = eye(t) - tt.*(s2u/sigs2);

%arg1 = n*log(2) - 0.5*n*t*log(2*pi) - 0.5*n*(t-1)*log(s2v) - 0.5*n*log(sigs2);
arg1 = - n*log(s_u) - 0.5*n*(t-1)*log(2*pi*s2v) - 0.5*n*log(t);

arg2 = - (t/(2*s2v))*sum(sse./t - (e_sr + s2v/(t*s_u)).^2);
%arg2 = - n*s_v^2/(2*s_u^2) - sum(- e - s_v^2/s_u)/s_u;
%arg3 = sum(log(normcdf(-e_sr.*(t*s_u/(s_v*sqrt(sigs2))))));
arg3 = sum(log(normcdf(-e_sr.*(t^0.5/s_v) -  s_v*(t^-0.5)/s_u)));

%{
spradzenie dla arg2
ec = reshape(e,t,n);
ECC = zeros(n,1);
for i = 1:n
    ECC(i) = ec(:,i)'*A*ec(:,i);
end
arg22 = sum(ECC)/(2*s2v);
%}
L = arg1 + arg2 + arg3;

%disp(arg1);
%disp(arg2);
%disp(arg3);
end

