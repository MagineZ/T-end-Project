function v = epiWACFM(X, epsi,qsi,m)
w_diff = inf;
v = mean(X,2);
[p,n] = size(X);
w_old = zeros(1,n);
options = optimoptions('linprog','Display','off');
while w_diff > qsi
    w_new = sum(max(abs((X-repmat(v,1,n)))-epsi,0)).^(1/(1-m))/sum(sum(max(abs((X-repmat(v,1,n)))-epsi,0)).^(1/(1-m)));
    for j = 1:p
        prob = optimproblem;
        lambda_j = optimvar('lambda_j',[2,n]);
        lambda_j.LowerBound = 0;
        lambda_j.UpperBound = repmat(w_new.^m,2,1);
        prob.Objective = lambda_j(1,:)*(X(j,:)'+epsi)+lambda_j(2,:)*(epsi-X(j,:)');
        cons1 = sum(lambda_j(1,:)) - sum(lambda_j(2,:)) == 0;
        prob.Constraints.cons1 = cons1;
        sol = solve(prob,'Options',options);
        lambda = sol.lambda_j;
        card_j = sum(((lambda(1,:)>0)&(lambda(1,:)<(w_new.^m)))|((lambda(2,:)>0) & (lambda(2,:)<(w_new.^m))));
        v(j) = (sum((X(j,:)+epsi).*((lambda(1,:)>0)&(lambda(1,:)<(w_new.^m)))) + sum((X(j,:)-epsi).*((lambda(2,:)>0)&(lambda(2,:)<(w_new.^m)))))/card_j;
    end
    w_diff = norm(w_new-w_old);
    w_old = w_new;
end