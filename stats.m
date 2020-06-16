% model statistics

BSFA.coefficients(:,1) = mean(gibbs.b(:,od:end),2);
BSFA.coefficients(:,2) = std(gibbs.b(:,od:end),0,2);

BSFA.sigma_v(:,1) = mean(gibbs.s_v(1,od:end),2);
BSFA.sigma_v(:,2) = std(gibbs.s_v(1,od:end),0,2);

BSFA.sigma_u(:,1) = mean(gibbs.s_u(1,od:end),2);
BSFA.sigma_u(:,2) = std(gibbs.s_u(1,od:end),0,2);

if strcmp(name,'nhn true') || strcmp(name,'nex true')
    BSFA.sigma_a(:,1) = mean(gibbs.s_a(1,od:end),2);
    BSFA.sigma_a(:,2) = std(gibbs.s_a(1,od:end),0,2);
    
    BSFA.a(:,1) = mean(gibbs.a(:,od:end),2);
    BSFA.a(:,2) = std(gibbs.a(:,od:end),0,2);
end

BSFA.u(:,1) = mean(gibbs.u(:,od:end),2);
BSFA.u(:,2) = std(gibbs.u(:,od:end),0,2);

BSFA.ef(:,1) = mean(exp(-gibbs.u(:,od:end)),2);
BSFA.ef(:,2) = std(exp(-gibbs.u(:,od:end)),0,2);



if make_plots == 1
    subplot(2,1,1);
    plot(gibbs.b(1,:));
    title('Intercept: sequential plot (all cycles)');
    subplot(2,1,2);
    rlanc = tot_cycles-burnin_cycles;
    pomoc2 = (1:rlanc);
    int = BSFA.coefficients(1,:);
    ss = gibbs.b(1,od:end);
    kusum = cumsum(ss)./pomoc2;
    CUSUM = (kusum-int(1))./int(2);
    
    bch_var = randn(1,rlanc);
    bch_var = (bch_var-mean(bch_var))./std(bch_var);
    bch_var = int(2).*bch_var+int(1);
    bch_kusum = cumsum(bch_var(1,1:rlanc))./pomoc2;
    bch_pth = (bch_kusum-int(1))./int(2);
    
    plot(pomoc2, CUSUM, pomoc2, bch_pth);
    title('Intercept: CUSUM plot for accepted draws')
    a2 = min(CUSUM(1,round(rlanc/8):rlanc));
    a1 = min(bch_pth(1,round(rlanc/8):rlanc));
    if a1<a2
        a2 = a1;
    end
    a2 = 2*a2;
    b = max(CUSUM(1,round(rlanc/8):rlanc));
    b1 = max(bch_pth(1,round(rlanc/8):rlanc));
    if b1>b
        b = b1;
    end
    b = 2*b;
    ylim([a2 b]);
    xlim([0, rlanc]);
end
