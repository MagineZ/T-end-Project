function [O,O2] = compute_morph_nonlocal_median2(I,I2, beats, N, fs)
%
% I: column vector
%

O = zeros(size(I));
O2 = zeros(size(I2));
% IdeVnlm = I;
% Ev = fs*0.08;
	% step1: find R peaks (in variable beats)

	% step2: group beats
RRI = beats(2:end) - beats(1:end-1) ;	
RRI = [RRI(1) RRI] ;
%MaximalQT = ceil(median(RRI)/2) ;
MaximalQT = ceil(quantile(RRI,.95)/2) ;

tmp = find(beats>MaximalQT) ;
beats = beats(tmp) ;
RRI = RRI(tmp) ;

tmp = find(beats+MaximalQT<=length(I)) ;
beats = beats(tmp) ;
RRI = RRI(tmp) ;

V = [] ;
V2 = [];
II = [];

for ii = 1:length(RRI)
		% the QT interval is smaller than the RR interval
		% we should also use the information that the QT interval is
		% proportional to the RR interval
	idx = beats(ii)-MaximalQT: beats(ii) + MaximalQT ; %min(MaximalQT, RRI(ii)-Ev) ;
    V(:,ii) = [I(idx)] ;
    V2(:,ii) = [I2(idx)] ;
    II = [II; idx] ;
end

Nbeat = size(V,2) ;
%Z = mandist(V',V) ;
Z = squareform(pdist(V')) ; %, 'cityblock')) ;
%Z = squareform(pdist(RRI'));
for i = 1:Nbeat; Z(i,i)=Inf; end

for i = 1:Nbeat
	tmpI = max(1,i-N*6) : min(Nbeat,i+N*6) ;
    z = sort(Z(i,tmpI),'ascend');
    th = z(N-1);
    Nidx = [i tmpI(1)-1+find(Z(i,tmpI)<=th)];
    Vx = V(:,Nidx); 
    Vx2 = V2(:,Nidx); 
    m = median(Vx, 2) ;
    m2 = median(Vx2, 2) ;
    if i == 1
        left_overlap = 2;
        right_overlap = length(intersect(II(i+1,:),II(i,:)));
    elseif i == Nbeat
        left_overlap = length(intersect(II(i-1,:),II(i,:)));
        right_overlap = 2;
    else
        left_overlap = length(intersect(II(i-1,:),II(i,:)));
        right_overlap = length(intersect(II(i+1,:),II(i,:)));
    end
    
    if left_overlap <= 1; left_overlap = 2; end
    if right_overlap <= 1; right_overlap = 2; end
    
    W = ones(2*MaximalQT+1,1);
    W(1:left_overlap) = sin(linspace(0,pi/2,left_overlap)).^2;
    W(end:-1:end-right_overlap+1) = sin(linspace(0,pi/2,right_overlap)).^2;
    O(II(i,:)) = O(II(i,:)) + (W.*m)' ;
    O2(II(i,:)) = O(II(i,:)) + (W.*m2)' ;
end

