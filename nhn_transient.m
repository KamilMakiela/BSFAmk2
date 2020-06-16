%% normal-half-normal SFA model - Gibbs sampling implementation

%u = zeros(nt,1);
u = max(ols_res) - ols_res;
ltn = -u;
r010lnsqr = 10*ln_r0^2;
for i=1:tot_cycles
    % beta part
    calc_beta;
    
    % latent variables part
    % KOS:
    s2u = 1/(randg(0.5*nt+5)/(0.5*sum(u.^2)+r010lnsqr));
    % Tsionas: (G(0.5N0,0.5Q0), where N0=1 Q0=10^-4
    %s2u = 1/(randg(0.5*(nt+1))/(0.5*(sum(u.^2)+10^-4)));
    s_u = sqrt(s2u);
    u = kmdraw2(s2u*(X*beta-y)/(s2u+s2v),sqrt((s2v*s2u)/(s2v+s2u)),nt);

    ltn = -u; % here the only latent variable is inefficiency
    
    % saving part
    gibbs.s_v(1,i) = s_v;
    gibbs.s_u(1,i) = s_u;
    gibbs.b(:,i) = beta;
    gibbs.u(:,i) = u;
    gibbs.lnL(1,i) = lgL_nhnp(s_v, s_u, beta, X, y, nt, 1);
    if rem(i,skok) == 0
        waitbar(i/tot_cycles,h);
    end
end

%% Statistics

waitbar(i/tot_cycles,h,'Calculating statistics');
stats;