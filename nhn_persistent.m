%% normal-half-normal SFA model - Gibbs sampling implementation

%u = zeros(n,1);
u = max(ols_res) - ols_res;
u = reshape(u,t,n);
u = mean(u)';
ltn = -kron(u,i_t);

for i=1:tot_cycles
    % beta part
    calc_beta;
    
    % latent variables part 
    % KOS:
    s2u = 1/(randg(0.5*n+5)/(0.5*sum(u.^2)+10*ln_r0.^2));
    s_u = sqrt(s2u);
    u = kmdraw2(s2u*(X_sr*beta-y_sr)/(s2u+s2v/t), sqrt((s2v*s2u/t)/(s2v/t+s2u)), n);
    ltn = -kron(u,i_t); % here the only latent variable is inefficiency
    
    % saving part
    gibbs.s_v(1,i) = s_v;
    gibbs.s_u(1,i) = s_u;
    gibbs.b(:,i) = beta;
    gibbs.u(:,i) = u;
    gibbs.lnL(1,i) = lgL_nhnp(s_v, s_u, beta, X, y, n, t);
    if rem(i,skok) == 0
        waitbar(i/tot_cycles,h);
    end
end

%% Statistics

waitbar(i/tot_cycles,h,'Calculating statistics');
stats;