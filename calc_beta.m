%% this script calculates beta part of Gibbs

s2v = 1/(randg(0.5*(a0+ nt))/(0.5*(n0+(y-X*beta-ltn)'*(y-X*beta-ltn))));
s_v = sqrt(s2v);
if beta_prior == 1  % informative prior
    cov_b = inv(C+XtrX./s2v);
    mean_b = (cov_b*(Cb+X'*(y-ltn)./s2v))';
else % uninformative prior
    cov_b = s2v*XtrX_inv;
    mean_b = (XtrX_inv*X'*(y-ltn))';
end
beta = mvnrnd(mean_b,cov_b)';