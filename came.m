% Function starts here
%{ 
when file is loaded, run:

od = BSFA.settings.brnin_cycles+1;
do = BSFA.settings.tot_cycles;
mcmc = [gibbs.b(:,od:do)' gibbs.s_v(1,od:do)' gibbs.s_u(1,od:do)'];
ln_L = gibbs.lnL(1,od:do)';
pr_set = BSFA.settings;
X = BSFA.settings.x;
y = BSFA.settings.y;
t = 1;
E = came(mcmc, pr_set, X, y, t, ln_L);
%}

function E = came(mcmc, pr_set, x, y, t, Lgl)


[nt, dim_b] = size(x);
n = nt/t;
[iteracji, dim_all] = size(mcmc);

switch pr_set.name
    case {'nex t', 'nex p'} 
        ln_w = @lgL_nexp;
        p_su = @(p)log(invgampdf(exp(p), 1, -pr_set.ln_r0)*exp(p));
        %disp('next');
    case {'nhn t', 'nhn p'}
        ln_w = @lgL_nhnp;
        % KOS:
        p_su = @(p)log(invgampdf(exp(2*p), 5, 10*pr_set.ln_r0^2)*2*exp(2*p));
        % Tsionas: 
        %p_su = @(p)log(invgampdf(exp(2*p), 0.5*pr_set.Na, 0.5*pr_set.Qa)*2*exp(2*p));
    case 're'
        ln_w = @lgL_re;
        p_su = @(p)log(invgampdf(exp(2*p), 0.5*pr_set.Na, 0.5*pr_set.Qa)*2*exp(2*p));
    otherwise
        disp('dupa');
end

p_sv = @(p)log(invgampdf(exp(2*p), 0.5*pr_set.a0, 0.5*pr_set.n0)*2*exp(2*p));

%transformation for half-normal and exponential (s_v, s_u)
mcmc(:,dim_all-1) = log(mcmc(:,dim_all-1));
mcmc(:,dim_all) = log(mcmc(:,dim_all));

%disp('Wyznaczam min i max dla parametrow');
b = zeros(dim_all,4);
for i = 1:dim_all
    b(i,1) = min(mcmc(:,i));
    b(i,2) = max(mcmc(:,i));
    b(i,3) = mean(mcmc(:,i));
    b(i,4) = std(mcmc(:,i));
end

zakres_MCMC(:,1) = b(:,1); %minimum zbioru theta
zakres_MCMC(:,2) = b(:,2); %maximum zbioru theta
mju = b(:,3)';
Mcov = cov(mcmc);
% here I can increase dispersion
for i = 1:dim_all
    Mcov(i,i) = Mcov(i,i).*1;
end
%Mcor = corr(mcmc);

%disp('Wyznaczam min i max dla loglikelhood');
if nargin < 6
    Lgl = zeros(iteracji,1);
    for i=1:iteracji
        Lgl(i) = ln_w(exp(mcmc(i,dim_all-1)), exp(mcmc(i,dim_all)), mcmc(i,1:dim_b)', x, y, n, t);
    end
end
L_min = min(Lgl);
L_max = max(Lgl);
%disp(L_min);
%disp(L_max);

E=0;
licz_zera = 0;
iter = 20000;
%robie symulacje

disp('Calculating ln[p(y)] using CAME...');
liczba_pr = 1+1+1; %prior densities of: beta, s_v, s_u

pe = zeros(liczba_pr,1);
zero_flag = zeros(1,2);

for i = 1:iter
    draw = mvnrnd(mju,Mcov)';
    %draw = (mju + mvtrnd(Mcor,150))';
    L = ln_w(exp(draw(end-1)), exp(draw(end)), draw(1:dim_b), x, y, n, t);
    if isempty(draw(draw(:) <= zakres_MCMC(:,1) | zakres_MCMC(:,2) <= draw(:))) && L >= L_min && L <= L_max
        
        ze = log(mvnpdf(draw',mju, Mcov));
        %ze = log(mvtpdf(draw' - mju,Mcor,150));
        pe(1) = log(mvnpdf(draw(1:dim_b)',pr_set.b_prior,pr_set.Cinv_prior));
        pe(2) = p_sv(draw(end-1));
        
        pe(3) = p_su((draw(end)));
         
        E = E + exp(L+sum(pe)-ze);
        % zapisywanie
        %i_p(:,i) = pe;
        %i_L(:,i) = L;
    else
        licz_zera = licz_zera + 1;
        if ~isempty(draw(draw(:) <= zakres_MCMC(:,1) | zakres_MCMC(:,2) <= draw(:)))
            zero_flag(1) = zero_flag(1) +1;
        elseif L < L_min || L > L_max
            zero_flag(2) = zero_flag(2) +1;
        end
    end
end

E = log(E/iter);

disp([licz_zera zero_flag]);
disp(E);
end
