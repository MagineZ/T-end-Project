function  Om0 = fECG_linear_detrend(x0,current_beats)

% optimal shrinakge of ECG

RRI = diff(current_beats);
RRI = [RRI(1) RRI];

%MaximalQT = ceil(max(RRI)/2);
MaximalQTp = ceil((median(RRI)+2*iqr(RRI))*4/8);
MaximalQTt = ceil((median(RRI)+2*iqr(RRI))*4/8);
tmp = find(current_beats>MaximalQTp) ;
current_beats = current_beats(tmp) ;
RRI = RRI(tmp) ;

tmp = find(current_beats+MaximalQTt<=length(x0)) ;
current_beats = current_beats(tmp) ;
RRI = RRI(tmp) ;

V = [];
II = [];
count = 0;
for ii = 1:length(current_beats)
    count = count+1;
    % the QT interval is smaller than the RR interval
    % we should also use the information that the QT interval is
    % proportional to the RR interval
    idx = (current_beats(ii)-MaximalQTp): (current_beats(ii) + MaximalQTt);
    V(:,count) = x0(idx);
    II = [II;idx];
    
end

V = V - linear_approach(V);

Om0 = zeros(1,length(x0));

for s = 1:length(current_beats)
    XN1_to_add = V(:,s);
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
    Om_toadd = (W.*XN1_to_add)';
    Om0(II(s,:)) = Om0(II(s,:)) + Om_toadd;
end






