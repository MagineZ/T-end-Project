function  trend = optimal_detrend( x0, Om,Of,OptimalShrinkageOpt, num_nonlocal,wlen, ovlap_ratio, ifauto)
x = x0-Om-Of;
olen = round(wlen*ovlap_ratio);
gap = wlen - olen;
E0 = [];
E = [];
II = [];
n = floor(length(x0)/gap)-1;


for i = 1:n
    E0 = [E0;x0((i-1)*gap+(1:wlen))];
    E = [E;x((i-1)*gap+(1:wlen))];
    II = [II;(i-1)*gap+(1:wlen)];
end
[n_seg,~ ] = size(E);
n_cluster = round(n_seg/(10*num_nonlocal));
idx = kmeans(E,n_cluster);
E0 = E0';
E = E';
XN0_c = {};

for i = 1:n_cluster
    V = E0(:,idx == i);
    if sum(idx == i) ==1
        XN0_c{i} = V;
        continue;
    end
    
    [n_t, n_theta] = size(V);
    
    sigma = V - median(V,2)*(ones(1,size(V,2)));
    sigma = sqrt(var(sigma(:)));
    
    if n_theta>n_t
        beta0 = n_t/n_theta;
        [ae,be,ce] = svd(V./(sigma*sqrt(n_theta))) ;
        lambdaOSe = diag(be) ;
        if ifauto
            singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt);
        else
            singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt,1);
        end
        XN0 = sigma*sqrt(n_theta).*(ae*diag(singvals)*ce(:,1:n_t)');
    else
        beta0 = n_theta/n_t;
        [ae,be,ce] = svd(V'./(sigma*sqrt(n_t))) ;
        lambdaOSe = diag(be) ;
        if ifauto
            singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt);
        else
            singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt,1);
        end
        XN0 = sigma*sqrt(n_t).*(ae*diag(singvals)*ce(:,1:n_theta)')';
    end
    
    XN0_c{i} = XN0;
end


trend = zeros(1,length(x0));
for s = 1:n
    XN0 = XN0_c{idx(s)};
    nb_id = find(idx == idx(s));
    Z = sum((XN0 - XN0(:,nb_id == s)).^2);
    
    [~,ind] = sort(Z,'ascend');
    Nidx = nb_id(ind(1:min(num_nonlocal,length(ind))));
   
    XN = E0(:,Nidx);
    u1 = median(XN,2);
    %u1 = smooth(u1,round(wlen/5),'loess');
    XN_to_add = u1;
  
    %reconstruction
    
    if s == 1
        left_overlap = 2;
        right_overlap = length(intersect(II(s+1,:),II(s,:)));
    elseif s == n
        left_overlap = length(intersect(II(s-1,:),II(s,:)));
        right_overlap = 2;
    else
        left_overlap = length(intersect(II(s-1,:),II(s,:)));
        right_overlap = length(intersect(II(s+1,:),II(s,:)));
    end
    
    
    if left_overlap <= 1; left_overlap = 2; end
    if right_overlap <= 1; right_overlap = 2; end
    
    
    W = ones(wlen,1);
    W(1:left_overlap) = sin(linspace(0,pi/2,left_overlap)).^2;
    W(end:-1:end-right_overlap+1) = sin(linspace(0,pi/2,right_overlap)).^2;
    %Om_toadd = interp1(II(s,:),XN(:,loc),II(s,:));
    %Om_toadd = W'.*Om_toadd;
    trend_toadd = (W.*XN_to_add)';
       
    trend(II(s,:)) = trend(II(s,:)) + trend_toadd;
  
end
  

