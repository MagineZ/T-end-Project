function  [Om0,Om0_real] = ECG_shrinkage0( x0,x0_real, current_beats, OptimalShrinkageOpt, ifauto)

% optimal shrinakge of ECG
  
                
RRI = diff(current_beats);
RRI = [RRI(1) RRI];
w1 = 50;
w2 = 50;
MaximalQTp = ceil(quantile(RRI,0.95)*4/8);
MaximalQTt = ceil(quantile(RRI,0.95)*4/8);

tmp = find(current_beats>MaximalQTp) ;
current_beats = current_beats(tmp) ;
RRI = RRI(tmp) ;

tmp = find(current_beats+MaximalQTt<=length(x0)) ;
current_beats = current_beats(tmp) ;
RRI = RRI(tmp) ;



V=[];
V2 = [];
II = [];
II2 = [];
R = [];
S = [];
tRS = [];
count = 0;
res = x0_real - x0;
for ii = 1:length(current_beats)
    count = count+1;
    idx = (current_beats(ii)-w1): (current_beats(ii) + w2);
    idx2 = (current_beats(ii)-MaximalQTp): (current_beats(ii) + MaximalQTt);
    V(:,count) = x0(idx);
    V2(:,count) = res(idx2);
    II = [II;idx];
    II2 = [II2;idx2];
    
end

[n_t, n_theta] = size(V);
sigma = (V - (median(V,2)*(ones(1,size(V,2))))).^2;
sigma =  sqrt(sum(sum(sigma)/n_t)/n_theta);
%sigma = (V - (median(V,2)*(ones(1,size(V,2)))));
%sigma = sqrt(var(sigma(:)));
[n_t, n_theta] = size(V);

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

[n_t2, n_theta2] = size(V2);
sigma2 = (V2 - (median(V2,2)*(ones(1,size(V2,2))))).^2;
sigma2 =  sqrt(sum(sum(sigma2)/n_t2)/n_theta2);

if n_theta2>n_t2
    beta0 = n_t2/n_theta2;
    [ae,be,ce] = svd(V2./(sigma2*sqrt(n_theta2))) ;
    lambdaOSe = diag(be) ;
    if ifauto
        singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt);
    else
        singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt,1);
    end
    XN02 = sigma2*sqrt(n_theta2).*(ae*diag(singvals)*ce(:,1:n_t2)');
else
    beta0 = n_theta2/n_t2;
    [ae,be,ce] = svd(V2'./(sigma2*sqrt(n_t2))) ;
    lambdaOSe = diag(be) ;
    if ifauto
        singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt);
    else
        singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt,1);
    end
    XN02 = sigma2*sqrt(n_t2).*(ae*diag(singvals)*ce(:,1:n_theta2)')';
end






Om0 = zeros(1,length(x0));
Om0_real = zeros(1,length(x0_real));

for s = 1:length(current_beats)

    XN_to_add = XN0(:,s);
    XN2_to_add = XN02(:,s);
    %reconstruction
    
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
    %Om_toadd = interp1(II(s,:),XN(:,loc),II(s,:));
    %Om_toadd = W'.*Om_toadd;
    %Om_toadd = (W.*XN_to_add)';
    Om_toadd = (XN_to_add)';
    Om_real_toadd = (W.*XN2_to_add)';
    
    Om0(II(s,:)) = Om0(II(s,:)) + Om_toadd;
    Om0_real(II2(s,:)) = Om0_real(II2(s,:)) + Om_real_toadd;
    
end
Om0_real = Om0_real+Om0;
Om0_real(Om0_real == 0) = x0_real(Om0_real == 0);
Om0 = ECG_detrend(Om0_real,91,0,0);






