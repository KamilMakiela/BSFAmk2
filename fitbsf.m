function BSFA = fitbsf(x,y,name,n)
% Function for Bayesian estiamtion of several most popular 
% Stochstic Frontier models. 
% 
% (c) 2020 Kamil Makiela
%
% Syntax:
%   BSFA = fitbsf(x,y,name);
%   BSFA = fitbsf(x,y,name,n);
% Description:
%   BSFA = fitbsf(x,y,name) returns Bayesian fit of a Stochastic Frotnier 
%   production-oriented model, specified in 'name' using data in 
%   x (explanatory variables) and y (the dependent). In its simplest form 
%   the (logged) model is: y = intercept + x*beta + v - u; where u is the
%   nonnegative component of inefficiency, so that r=exp(-u) is the
%   standard efficiency measure known in SFA literature 
%   BSFA = fitbsf(x,y,name,n) returns Bayesian fit of an SF model for 
%   panel data (balanced panels only). Observations are to be stacked as
%   object:time (not period:time!); see example below. 
% inputs: 
%   x - explanatory variables in comlums (excluding column with "1"'s!!)
%   y - the dependent; column vector
%   NOTE: both y and x are assumed to have observations in logs of their 
%   original values so that there is a (log)linear dependency assumed
%   bewteen x's and y
%   name - name of SFA sampling model; you can choose from:
%        nhn t - normal-half-normal with transient efficiency
%        nhn p - normal-half-normal with persistent efficiency (needs "n")
%        nex t - normal-exponential with transient efficiency
%        nex p - normal-exponential with persistent efficiency
%        nhn true - "true" random effects model (with transient [half
%        normal] efficiency and a symmetric object sepcific effect (needs
%        "n")
%        nex true - "true" random effects model (with transient
%        [exponential] efficiency and a symmetric object sepcific effect  
%        (needs "n")
%        re - simple Bayesian random-effects model (non-SFA)
%   n - (optional) number of objects (for panel data); it is assumed
%       that observations are stacked in "object->time" order. 
% output:
% BSFA - matlab structure with results:
%   BSFA.coeffitients - posterior estimates and standard deviations of the
%   function coefficients; the first coefficient is the intercept; the rest
%   are in the order according to varaibles in matrics x
%   BSFA.sigma_v - std of the symmetric disturbance
%   BSFA.sigma_u - std of the inefficiency component 
%   BSFA.u - inefficiency estimates (or random effects in case of RE model)
%   BSFA.ef - efficiency estimates
%   BSFA.Loglikelihood - loglikelihood value of the Bayesian solution 
%   BSFA.model - name of the SFA model chosen
%   BSFA.settings - simulation and model settings
%   BSFA.sigma_a - (optional) std deviation of the additional symmetric
%                  object-specific disturbance (used in 'true' SFA models)
%   BSFA.a - (optional) estimated random effects in 'true' SFA models for
%            panel data 
% ADDITIONAL options: see "function settings" section below 
%
% Example
% estimating normal-half-normal SFA with persistent efficiency:
%   ExampleData;                    % load example data
%   bayes = fitbsf(x,y,'nhn p',n);  % run selected SFA model
% output structure 'bayes' contains results
% 
%% Function settings

beta_prior = 1; % if 1: informative prior on coefficients, uninformative otherwise
zapisz     = 0; % if 1: saves gibbs cycles, results and settings to current directory
make_plots = 0; % if 1: displays sequential and cusum plots for the intercept
calc_came  = 0; % if 1: the procedure also calculates the log of marginal data 
                %       density (aka marignal likelihood) based on:
                %       corrected arithemtic mean estimator (Pajor, 2017; Bayesian Analysis)
                %       and harmonic mean estimator (Geweke, 1992)
                %       use informative prior setting (beta_prior=1) when
                %       using this option!

tot_cycles =    25000; % total number of cycles: in theory, the more the better (smaller numeric error)
burnin_cycles =  5000; % burnin phase; usually 10k is enough for most applications

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Hiperparameters

k_par = size(x,2)+1;
a0 = 10^-4;
n0 = 10^-4;
b = zeros(k_par,1);
C = eye(k_par).*100^-2;
Cb = C*b;
% for u (if exists)
r0 = 0.75;
ln_r0 = log(r0);
% for a (if exists)
Qa = 10^-4;
Na = 1;
settings.ln_r0 = ln_r0;
settings.name = name;

if beta_prior == 1
    settings.a0 = a0;
    settings.n0 = n0;
    settings.b_prior = b';
    settings.Cinv_prior = inv(C);
    settings.Cb = Cb;
else
    settings.beta = 'uninformative';
end
settings.r0 = r0;
settings.brnin_cycles = burnin_cycles;
settings.tot_cycles = tot_cycles;
settings.ln_r0 = ln_r0;


%% preliminaries
nt = size(y,1);
t=1;
X = [ones(nt,1) x];
XtrX = X'*X;
XtrX_inv = inv(XtrX);
beta = XtrX_inv*(X'*y);
ols_res = y-X*beta;
beta(1) = beta(1) + max(ols_res);

skok = tot_cycles/100;
od = burnin_cycles +1;
%% For panel data
if nargin == 4
    t = nt/n;
    i_t = ones(t,1);
	i_n = ones(n,1);
    t_sqr = sqrt(t);
    
    %% Averages
    % the assumption is that we have observations in order of:
    % firm1 over 1 to t, then firm2 over 1 to t, and so on
    
    y_sr = mean(reshape(y,t,n))';
    X_sr = zeros(n,k_par);
    for i = 1:k_par
        X_sr(:,i) = mean(reshape(X(:,i),t,n))';
    end
    settings.y_sr = y_sr;
    settings.X_sr = X_sr;

end
settings.y = y;
settings.x = X;

%% main matrices
gibbs.lnL = zeros(1,tot_cycles);
gibbs.s_v = zeros(1,tot_cycles);
gibbs.s_u = zeros(1,tot_cycles);
gibbs.b = zeros(k_par,tot_cycles);
if strcmp(name,'nhn p') || strcmp(name,'nex p') || strcmp(name,'re') 
    gibbs.u = zeros(n,tot_cycles);

else
    gibbs.u = zeros(nt,tot_cycles);
end
settings.Qa = Qa;
settings.Na = Na;
%% Sampling
h = waitbar(0,'Running...');
switch name
    case 'nhn t'
        nhn_transient;
        BSFA.LogLikelihood = lgL_nhnp(BSFA.sigma_v(1), BSFA.sigma_u(1), BSFA.coefficients(:,1), X, y, nt, 1);
    case 'nhn p'
        nhn_persistent;
        BSFA.LogLikelihood = lgL_nhnp(BSFA.sigma_v(1), BSFA.sigma_u(1), BSFA.coefficients(:,1), X, y, n, t);
    case 'nex t'
        nex_transient;
        BSFA.LogLikelihood = lgL_nexp(BSFA.sigma_v(1), BSFA.sigma_u(1), BSFA.coefficients(:,1), X, y, nt, 1);
    case 'nex p'
        nex_persistent;
        BSFA.LogLikelihood = lgL_nexp(BSFA.sigma_v(1), BSFA.sigma_u(1), BSFA.coefficients(:,1), X, y, n, t);
    case 'nhn true'
        nhn_true;
    case 'nex true' %!!
        nex_true;
    case 're'
        re_panel;
        BSFA.LogLikelihood = lgL_re(BSFA.sigma_v(1), BSFA.sigma_u(1), BSFA.coefficients(:,1), X, y, n, t);

    otherwise
        disp('unknown model');
end

BSFA.model = name;
BSFA.settings = settings;
delete(h);

if calc_came == 1
    switch name
        case {'nhn t','nhn p','nex t','nex p', 're'}
            mcmc = [gibbs.b(:,od:end)' gibbs.s_v(1,od:end)' gibbs.s_u(1,od:end)'];
            BSFA.py_came = came(mcmc,settings,X,y,t,gibbs.lnL(1,od:end)');
            BSFA.py_hme = harmmean(gibbs.lnL(1,od:end));
            BSFA.BIC = (k_par+2)*log(nt) - 2*BSFA.LogLikelihood;
        otherwise
            disp('p(y) estimation is not avliable for this model; to be added in future release');
    end
end

if zapisz == 1
    plik = datestr(now);
    plik = strrep(plik,':','-');
    plik = strcat('bsfa-',plik);
    try
        save(plik,'BSFA','gibbs');
    catch
        disp('Unable to save file. It might be too big. Check Matlab settings and enable file comprsession');
    end
end