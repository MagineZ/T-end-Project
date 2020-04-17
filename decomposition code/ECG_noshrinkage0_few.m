function  [Om0,Om0_real] = ECG_shrinkage0( x0,x0_real, current_beats, OptimalShrinkageOpt, sigma_coeff,num_singval)

% optimal shrinakge of ECG
  
                
RRI = diff(current_beats);
RRI = [RRI(1) RRI];
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
R = [];
S = [];
tRS = [];
count = 0;
for ii = 1:length(current_beats)
    count = count+1;
    idx = (current_beats(ii)-MaximalQTp): (current_beats(ii) + MaximalQTt);
    V(:,count) = x0(idx);
    V2(:,count) = x0_real(idx);
    II = [II;idx];
    
end

[n_t, n_theta] = size(V);

sigma = (V - (median(V,2)*(ones(1,size(V,2))))).^2;
sigma =  sigma_coeff*sqrt(sum(sum(sigma)/n_t)/n_theta);
%sigma = (V - (median(V,2)*(ones(1,size(V,2)))));
%sigma = sqrt(var(sigma(:)));
[n_t, n_theta] = size(V);

if n_theta>n_t
    beta0 = n_t/n_theta;
    [ae,be,ce] = svd(V./(sigma*sqrt(n_theta))) ;
    singvals = diag(be) ;
    singvals((num_singval+1):end) = 0;
    XN0 = sigma*sqrt(n_theta).*(ae*diag(singvals)*ce(:,1:n_t)');
else
    beta0 = n_theta/n_t;
    [ae,be,ce] = svd(V'./(sigma*sqrt(n_t))) ;
    singvals = diag(be) ;
    singvals((num_singval+1):end) = 0;
    XN0 = sigma*sqrt(n_t).*(ae*diag(singvals)*ce(:,1:n_theta)')';
end


sigma = (V2 - (median(V2,2)*(ones(1,size(V2,2))))).^2;
sigma = sigma_coeff*sqrt(sum(sum(sigma)/n_t)/n_theta);
if n_theta>n_t
    beta0 = n_t/n_theta;
    [ae,be,ce] = svd(V2./(sigma*sqrt(n_theta))) ;
    singvals = diag(be) ;
    singvals((num_singval+1):end) = 0;
    XN02 = sigma*sqrt(n_theta).*(ae*diag(singvals)*ce(:,1:n_t)');
else
    beta0 = n_theta/n_t;
    [ae,be,ce] = svd(V2'./(sigma*sqrt(n_t))) ;
    singvals = diag(be) ;
    singvals((num_singval+1):end) = 0;
    XN02 = sigma*sqrt(n_t).*(ae*diag(singvals)*ce(:,1:n_theta)')';
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
    Om_toadd = (W.*XN_to_add)';
    Om_real_toadd = (W.*XN2_to_add)';
    
    Om0(II(s,:)) = Om0(II(s,:)) + Om_toadd;
    Om0_real(II(s,:)) = Om0_real(II(s,:)) + Om_real_toadd;
    
end

Om0(Om0 == 0) = x0(Om0 == 0);
Om0_real(Om0_real == 0) = x0_real(Om0_real == 0);







