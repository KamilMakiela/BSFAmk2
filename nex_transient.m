%% gibbs

%u = zeros(nt,1);
u = max(ols_res) - ols_res;
ltn = -u;

for i=1:tot_cycles
    % beta part
    calc_beta;
    
    % latent variables part
    s_u = 1/(randg(nt+1)/(sum(u)-ln_r0));
    u = kmdraw2(X*beta - y - s2v/s_u, s_v, nt);
    ltn = -u; % here the only latent variable is inefficiency
    
    % saving part
    gibbs.s_v(1,i) = s_v;
    gibbs.s_u(1,i) = s_u;
    gibbs.b(:,i) = beta;
    gibbs.u(:,i) = u;
    gibbs.lnL(1,i) = lgL_nexp(s_v, s_u, beta, X, y, nt, 1);
    if rem(i,skok) == 0
        waitbar(i/tot_cycles,h);
    end
end

%% Statistics

waitbar(i/tot_cycles,h,'Calculating statistics');
stats;