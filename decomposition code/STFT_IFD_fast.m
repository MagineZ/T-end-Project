function [tfr, ifd, tfrtic, t] = STFT_IFD_fast(x, alpha, Hop, h, Dh)
% Synchrosqueezing by Li Su, 2015

if size(h,2) > size(h,1)
    h = h';
end

if size(Dh,2) > size(Dh,1)
    Dh = Dh';
end

	% for tfr
N = length(-0.5+alpha:alpha:0.5);
Win_length = max(size(h));
TH = 7*N/size(h,1); 
tfrtic = linspace(0, 0.5, round(N/2))' ;
	% for tfrsq
% Lidx = ceil( (N/2)*(lowFreq/0.5) ) ; 
% Hidx = floor( (N/2)*(highFreq/0.5) ) ; 
% fLen = Hidx - Lidx + 1 ;

Overlap = Win_length-Hop;
Lh = floor((Win_length-1)/2);
t = Hop:Hop:floor(length(x)/Hop)*Hop;
x_Frame = zeros(N, length(t));
tf2 = zeros(N, length(t));
for ii = 1:length(t)
    ti = t(ii); 
    tau = -min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,length(x)-ti]);
    indices= rem(N+tau,N)+1;
    norm_h=norm(h(Lh+1+tau)); 

% 	tf0 = zeros(N, 1) ; tf1 = zeros(N, 1) ;
    x_Frame(indices,ii) = (x(ti+tau)-mean(x(ti+tau)))'.*conj( h(Lh+1+tau)) /norm_h;
    tf2(indices,ii) = (x(ti+tau)-mean(x(ti+tau)))'.*conj( Dh(Lh+1+tau)) /norm_h;
end
% x_Frame = buffer(x, Win_length, Overlap);
% x_Frame = x_Frame(:,ceil(Overlap/Hop/2)+1:end);
% x_Frame = x_Frame - repmat(mean(x_Frame),[size(x_Frame,1) 1]);
Stime = round(Hop-Win_length/2+1)+ceil(Overlap/Hop/2)*Hop;

% tfr = x_Frame.*repmat(h, [1 size(x_Frame,2)]); 	% for h
tfr = fft(x_Frame, N, 1);
tfr = tfr(1:round(N/2),:);

% tf2 = x_Frame.*repmat(Dh, [1 size(x_Frame,2)]); 	% for h
tf2 = fft(tf2, N, 1);
tf2 = tf2(1:round(N/2),:);

		% get the first order omega
omega = zeros(size(tf2)) ;
avoid_warn = find(tfr~=0);
omega(avoid_warn) = imag(N*tf2(avoid_warn)./tfr(avoid_warn)/(2.0*pi));
ifd = omega;
% omega = round(omega);
% % omega(abs(omega)>TH/2)=0;
% 
% OrigIndex = repmat((1:round(N/2))', [1 size(tfr,2)]);
% omega(OrigIndex - omega < 1 | OrigIndex - omega > round(N/2))=0;
% % ReasIndex = OrigIndex - omega;
% 
% Ex = mean(abs(x).^2);
% Threshold = 1.0e-8*Ex;	% originally it was 1e-6*Ex
% tfr(abs(tfr) < Threshold) = 0;
% 
% totLength = size(tfr,1)*size(tfr,2);
% rtfr = accumarray((1:totLength)'-omega(:),tfr(:));
% rtfr = [rtfr; zeros(totLength-length(rtfr),1)];
% 
% rtfr = reshape(rtfr,size(tfr,1),size(tfr,2));

end