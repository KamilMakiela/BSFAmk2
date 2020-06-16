% normal-half-normal SFA model - Gibbs sampling implementation

gibbs.s_a = zeros(1,tot_cycles);
gibbs.a = zeros(n,tot_cycles);

%u = zeros(nt,1);
u = max(ols_res) - ols_res;
a = zeros(n,1);
ltn = kron(a,i_t)-u;

for i=1:tot_cycles
    % beta part
    calc_beta;
    
    % latent variables part
    % a 
    s2a = (Qa+a'*a)/chi2rnd(n+Na,1);
    a_fala = mean(reshape(y-X*beta+u,t,n),1)';
    mi_a = (s2a/(s2v/t+s2a)).*a_fala;
    a = normrnd(mi_a,sqrt((s2a*s2v/t)/(s2v/t+s2a)));
    
    % u
    s_u = 1/(randg(nt+1)/(sum(u)-ln_r0));
	u = kmdraw2(X*beta + kron(a,i_t)- y - s2v/s_u, s_v, nt);
        
    ltn = kron(a,i_t)-u; % here the only latent variable is inefficiency
    
    % saving part
	gibbs.s_v(1,i) = s_v;
	gibbs.s_u(1,i) = s_u;
	gibbs.s_a(1,i) = sqrt(s2a);
	gibbs.b(:,i) = beta;
	gibbs.u(:,i) = u;
	gibbs.a(:,i) = a;
    if rem(i,skok) == 0
        waitbar(i/tot_cycles,h);
    end
end

%% Statistics

waitbar(i/tot_cycles,h,'Calculating statistics');
stats;
