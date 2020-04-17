function y = linear_approach(x)
[n,m] = size(x);
t = repmat((1:n)',1,m);
alpha = (mean(x.*t) - mean(t).*mean(x))./(mean(t.^2)-mean(t).^2);
beta = mean(x) - alpha.*mean(t);
alpha = repmat(alpha,n,1);
beta = repmat(beta,n,1);
y = alpha.*t+beta;

