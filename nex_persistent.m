%% Gibbs

%u = zeros(n,1);
u = max(ols_res) - ols_res;
u = reshape(u,t,n);
u = mean(u)';
ltn = -kron(u,i_t);

for i=1:tot_cycles
    % beta part
    calc_beta;
    
    % latent variables part
    s_u = 1/(randg(n+1)/(mean(u)*n-ln_r0));
    u = kmdraw2(X_sr*beta - y_sr - (s2v/s_u)/t, s_v/t_sqr, n);
    ltn = -kron(u,i_t); % here the only latent variable is inefficiency
    
    % saving part
    gibbs.s_v(1,i) = s_v;
    gibbs.s_u(1,i) = s_u;
    gibbs.b(:,i) = beta;
    gibbs.u(:,i) = u;
    gibbs.lnL(1,i) = lgL_nexp(s_v, s_u, beta, X, y, n, t);
    if rem(i,skok) == 0
        waitbar(i/tot_cycles,h);
    end
end

%% Statistics

waitbar(i/tot_cycles,h,'Calculating statistics');
stats;
