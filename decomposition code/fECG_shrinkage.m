function  [Om0,Om0_real] = fECG_shrinkage( x0,x0_real, current_beats, OptimalShrinkageOpt,num_nonlocal, ifauto)

% optimal shrinakge of ECG

RRI = diff(current_beats);
RRI = [RRI(1) RRI];
%MaximalQT = ceil(max(RRI)/2);
MaximalQTp = ceil(median(RRI)*4/8);
MaximalQTt = ceil(median(RRI)*4/8);

tmp = find(current_beats>MaximalQTp) ;
current_beats = current_beats(tmp) ;
RRI = RRI(tmp) ;

tmp = find(current_beats+MaximalQTt<=length(x0)) ;
current_beats = current_beats(tmp) ;
RRI = RRI(tmp) ;



V=[];
V2 = [];
II = [];
R = [];
S = [];
tRS = [];
count = 0;
for ii = 1:length(current_beats)
    count = count+1;
    % the QT interval is smaller than the RR interval
    % we should also use the information that the QT interval is
    % proportional to the RR interval
    idx = (current_beats(ii)-MaximalQTp): (current_beats(ii) + MaximalQTt);
    V(:,count) = x0(idx);
    V2(:,count) = x0_real(idx);
    II = [II;idx];
    
end

[n_t, n_theta] = size(V);

sigma = V2 - median(V2,2)*(ones(1,size(V2,2)));
sigma = 3*sqrt(var(sigma(:)));
%sigma = median(abs(x0(current_beats)))/10;
%sigma = 1;

[n_t, n_theta] = size(V);
%V = V - linear_approach(V);
%V2 = V2 - linear_approach(V2);
if n_theta>n_t
    beta0 = n_t/n_theta;
    [ae,be,ce] = svd(V2./(sigma*sqrt(n_theta))) ;
    lambdaOSe = diag(be) ;
    if ifauto
        singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt);
    else
        singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt,1);
    end
    XN0 = sigma*sqrt(n_theta).*(ae*diag(singvals)*ce(:,1:n_t)');
else
    beta0 = n_theta/n_t;
    [ae,be,ce] = svd(V2'./(sigma*sqrt(n_t))) ;
    lambdaOSe = diag(be) ;
    if ifauto
        singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt);
    else
        singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt,1);
    end
    XN0 = sigma*sqrt(n_t).*(ae*diag(singvals)*ce(:,1:n_theta)')';
end






Om0 = zeros(1,length(x0));
Om0_real = zeros(1,length(x0));
%XN0_trend = linear_approach(XN0);
%wlen = round(median(RRI)/5);
%wlen = round(abs(100*x0_RS/14));
%swl = round(size(V,1)/(20*beta));
for s = 1:length(current_beats)
    
   
    ZZZ = XN0 - XN0(:,s)*(ones(1,size(XN0,2)));
    ZZZ = sum(ZZZ.^2);
    %ZZZ = abs(RRI_m - RRI_m(s));
    ZZZ(s) = inf;
    %ZZZ(1:length(current_beats)) = inf;
    [~,ind] = sort(ZZZ,'ascend');
    
    
    Nidx =[s,ind(1:(num_nonlocal-1))];
    
    
    XN = V(:,Nidx);
    XN2 = V2(:,Nidx);
    sigma = XN2 - median(XN2,2)*(ones(1,size(XN2,2)));
    sigma = 2*sqrt(var(sigma(:)));
    [p, n] = size(XN2);
    if n>p
        beta0 = p/n;
        [ae,be,ce] = svd(XN2./(sigma*sqrt(n))) ;
        lambdaOSe = diag(be) ;
        if ifauto
            singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt);
        else
            singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt,1);
        end
        XN2 = sigma*sqrt(n).*(ae*diag(singvals)*ce(:,1:p)');
    else
        beta0 = n/p;
        [ae,be,ce] = svd(XN2'./(sigma*sqrt(p))) ;
        lambdaOSe = diag(be) ;
        if ifauto
            singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt);
        else
            singvals = optimal_shrinkage(lambdaOSe,beta0,OptimalShrinkageOpt,1);
        end
        XN2 = sigma*sqrt(p).*(ae*diag(singvals)*ce(:,1:n)')';
    end

    
    
    XN_to_add = XN2(:,1);
    XN2_to_add = XN2(:,1);
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
    W(end:-1:end-right_overlap+1) = cos(linspace(0,pi/2,right_overlap)).^2;
    %Om_toadd = interp1(II(s,:),XN(:,loc),II(s,:));
    %Om_toadd = W'.*Om_toadd;
    Om_toadd = (W.*XN_to_add)';
    Om_real_toadd = (W.*XN2_to_add)';
    
    Om0(II(s,:)) = Om0(II(s,:)) + Om_toadd;
    Om0_real(II(s,:)) = Om0_real(II(s,:)) + Om_real_toadd;
    
end






