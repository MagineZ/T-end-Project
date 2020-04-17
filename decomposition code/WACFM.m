function v = WACFM(X, qsi,m)
w_diff = inf;
v = mean(X,2);
[~,N] = size(X);
w_old = zeros(1,N);
while w_diff > qsi
    w_new = sqrt(sum((X-repmat(v,1,N)).^2)).^(2/(1-m))/sum(sqrt(sum((X-repmat(v,1,N)).^2)).^(2/(1-m)));
    v = (X*w_new')/sum(w_new);
    w_diff = norm(w_new-w_old);
    w_old = w_new;
end