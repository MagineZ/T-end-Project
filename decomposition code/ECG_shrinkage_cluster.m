function Om0_real = ECG_shrinkage_cluster(x0_real, current_beats, OptimalShrinkageOpt, num_nonlocal, ifauto)

% optimal shrinakge of ECG

RRI = diff(current_beats);
RRI = [RRI(1) RRI];

%MaximalQT = ceil(max(RRI)/2);
MaximalQTp = ceil(median(RRI)*4/8);
MaximalQTt = ceil(median(RRI)*4/8);
tmp = find(current_beats>MaximalQTp) ;
current_beats = current_beats(tmp) ;
RRI = RRI(tmp) ;

tmp = find(current_beats+MaximalQTt<=length(x0_real)) ;
current_beats = current_beats(tmp) ;
RRI = RRI(tmp) ;


V2 = [];
II = [];
count = 0;
for ii = 1:length(current_beats)
    count = count+1;
    % the QT interval is smaller than the RR interval
    % we should also use the information that the QT interval is
    % proportional to the RR interval
    idx = (current_beats(ii)-MaximalQTp): (current_beats(ii) + MaximalQTt);
    V2(:,count) = x0_real(idx);
    II = [II;idx];
    
end


[~,n_seg] = size(V2);
n_cluster = round(sqrt(n_seg/2));
idx = kmeans(RRI',n_cluster);

XN0_c = {};

for i = 1:n_cluster
    E = V2(:,idx == i);
    [n_t, n_theta] = size(E);
    
    sigma = (E - (median(E,2)*(ones(1,size(E,2))))).^2;
    sigma = sqrt(sum(sum(sigma)/n_t)/n_theta);
    if sum(idx == i) ==1
        XN0_c{i} = E;
        continue;
    end
    
    
    
    %sigma = E - median(E,2)*(ones(1,size(E,2)));
    %sigma = sqrt(var(sigma(:)));
    
    if n_theta>n_t
        beta0 = n_t/n_theta;
        [ae,be,ce] = svd(E./(sigma*sqrt(n_theta))) ;
        lambdaOSe = diag(be) ;
        if ifauto
            singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt);
        else
            singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt,1);
        end
        XN0 = sigma*sqrt(n_theta).*(ae*diag(singvals)*ce(:,1:n_t)');
    else
        beta0 = n_theta/n_t;
        [ae,be,ce] = svd(E'./(sigma*sqrt(n_t))) ;
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

Om0_real = zeros(1,length(x0_real));

for s = 1:length(current_beats)
    
    XN0 = XN0_c{idx(s)};
    
    nb_id = find(idx == idx(s));
   
    Z = sum((XN0 - XN0(:,nb_id == s)).^2);
    [~,ind] = sort(Z,'ascend');
    Nidx = nb_id(ind(1:min(num_nonlocal,length(ind))));
    XN2 = V2(:,Nidx);
    u2 = median(XN2,2);
    XN2_to_add = u2;
    
    %XN2_to_add = XN0(:,nb_id == s);%u2;
    if s == 1
        left_overlap = 2;
        right_overlap = length(intersect(II(s+1,:),II(s,:)));
    elseif s == length(current_beats)
        left_overlap = length(intersect(II(s-1,:),II(s,:)));
        right_overlap = 2;
    else
        left_overlap = length(intersect(II(s-1,:),II(s,:)));
        right_overlap = length(intersect(II(s+1,:),II(s,:)));
    end
    
    
    if left_overlap <= 1; left_overlap = 2; end
    if right_overlap <= 1; right_overlap = 2; end
    
    
    W = ones(MaximalQTp+MaximalQTt+1,1);
    W(1:left_overlap) = sin(linspace(0,pi/2,left_overlap)).^2;
    W(end:-1:end-right_overlap+1) = sin(linspace(0,pi/2,right_overlap)).^2;
    Om_real_toadd = (W.*XN2_to_add)';
    Om0_real(II(s,:)) = Om0_real(II(s,:)) + Om_real_toadd;
    
end






