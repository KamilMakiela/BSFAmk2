% normal-half-normal SFA model - Gibbs sampling implementation

u = zeros(n,1);
ltn = kron(u,i_t);

for i=1:tot_cycles
    % beta part
    calc_beta;
    
    % latent variables part
    % a
    s2u = (Qa+u'*u)/chi2rnd(n+Na,1);
    a_fala = mean(reshape(y-X*beta,t,n),1)';
    mi_a = (s2u/(s2v/t+s2u)).*a_fala;
    u = normrnd(mi_a,sqrt((s2u*s2v/t)/(s2v/t+s2u)));
    
    ltn = kron(u,i_t); % here the only latent variable is inefficiency
    
    % saving part
    gibbs.s_v(1,i) = s_v;
    gibbs.s_u(1,i) = sqrt(s2u);
    gibbs.b(:,i) = beta;
    gibbs.u(:,i) = u;
    gibbs.lnL(1,i) = lgL_re(s_v, sqrt(s2u), beta, X, y, n, t);
    if rem(i,skok) == 0
        waitbar(i/tot_cycles,h);
    end
end

%% Statistics

waitbar(i/tot_cycles,h,'Calculating statistics');
stats;