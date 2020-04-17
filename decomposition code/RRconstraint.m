function y = RRconstraint(beats, xf, fs, RS)

RRI = diff(beats);
FP = [];

for ii = 1:length(RRI)
    if RRI(ii) <= RS*fs
			% we should take the "larger pick" as the R peak
		if abs(xf(beats(ii))) <= abs(xf(beats(ii+1)))
        	FP = [FP beats(ii)] ;
		else
			FP = [FP beats(ii+1)];
% 			RRI(ii+1) = RR(ii+1) + RR(ii) ;
		end
    end
end
y = setdiff(beats, FP);
