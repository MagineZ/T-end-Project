function  Om0 = ECG_shrinkage0_qual1( x0,current_beats,beats_todo,tsqi,q,num_nonlocal, OptimalShrinkageOpt, sigma_coeff,ifauto)

% optimal shrinakge of ECG
  
                
RRI = diff(current_beats);
RRI = [RRI(1) RRI];
MaximalQTp = ceil(quantile(RRI,0.95)*4/8);
MaximalQTt = ceil(quantile(RRI,0.95)*4/8);

tmp = find(current_beats>MaximalQTp) ;
current_beats = current_beats(tmp) ;
RRI = RRI(tmp) ;
tsqi = tsqi(tmp);
tmp = find(current_beats+MaximalQTt<=length(x0)) ;
current_beats = current_beats(tmp) ;
RRI = RRI(tmp) ;
tsqi = tsqi(tmp);


V=[];

II = [];
count = 0;
for ii = 1:length(current_beats)
    count = count+1;
    idx = (current_beats(ii)-MaximalQTp): (current_beats(ii) + MaximalQTt);
    V(:,count) = x0(idx);
    II = [II;idx];
    
end

%{
[t,~] = size(V);
t = 1:t;
for s = 1:40
    plot(t/250,V(:,s),'Color',[0,0,0]+s/40); hold on;
end
plot(t/250,median(V(:,1:40),2),'r', 'LineWidth',3)
hold off
pause

%}

Om0 = zeros(1,length(x0));
bsqi_Tend_pb = zeros(1,length(current_beats));

for s = 1:length(current_beats)
    
    if ismember(current_beats(s),beats_todo) && tsqi(s)<q
        bsqi_Tend_pb(s) = 1;
        Z = abs(RRI - RRI(s));
        Z(s) = inf;
        [~,ind] = sort(Z,'ascend');
        Nidx = [s, ind(1:min(num_nonlocal,length(Z)-1))];
        V2 = V(:,Nidx);
        [n_t, n_theta] = size(V2);
   
        sigma = (V2 - (median(V2,2)*(ones(1,size(V2,2))))).^2;
        sigma =  sigma_coeff*sqrt(sum(sum(sigma)/n_t)/n_theta);

   

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


        XN_to_add = XN0(:,1);
    
    else 
        XN_to_add = V(:,s);
    end
    
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
   
    Om0(II(s,:)) = Om0(II(s,:)) + Om_toadd;
   
end

Om0(Om0 == 0) = x0(Om0 == 0);
%Om0_real(Om0_real == 0) = x0_real(Om0_real == 0);







