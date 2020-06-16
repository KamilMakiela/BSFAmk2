# BSFAmk2
A simple to use function for Bayesian Stochastic Frontier Analysis 

Description From fitbsf.m file: 

% Function for Bayesian estiamtion of several most popular Stochstic Frontier models. 
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
