function [ceps0, tceps] = cepstrum_convert(tfr, tfrtic, g, fs, Tc, num_s, HighFreq, LowFreq)

if g~=0
    ceps = real(ifft(abs(tfr).^g,2*size(tfr,1),1));
else
    ceps = real(ifft(log(abs(tfr)),2*size(tfr,1),1));
end
for mi=1:size(ceps,2)
if exist('tra', 'var')
    if  tra>=1
        ceps(1:max([round(1/HighFreq) tra]),mi)=0;
    else
        ceps(1:round(1/HighFreq),mi)=0;
    end
else
    ceps(1:round(1/HighFreq),mi)=0;
end
end

ceps(isnan(ceps)|isinf(ceps))=0;

ceps = ceps(1:round(1/LowFreq),:);
ceps0=ceps;
ceps =  interp1(1:size(ceps,1), ceps, 1:0.1:size(ceps,1));
tceps=zeros(length(tfrtic),size(ceps,2));
freq_scale=10.*fs./(1:size(ceps,1)-1);
empt = [];
for ii=2:length(tfrtic)-1
    p_index = find(freq_scale > (tfrtic(ii-1)+tfrtic(ii))*fs/2 & freq_scale < (tfrtic(ii+1)+tfrtic(ii))*fs/2);
    if isempty(p_index)
        empt = [empt; ii];
    else
        % previous code: ii.^2. Changed to ii.^0 for SampTA2017 manuscript
        tceps(ii,:)=sum(ceps(p_index,:),1).*(ii.^0);
    end
end
tceps(tceps<Tc)=0;
tceps_new = tceps;
if num_s>1
    for fi = round(LowFreq*fs*num_s)+1:size(tceps,1)
        for kk = 2:num_s
            tceps_new(fi,:)=tceps(fi,:).*tceps(max([1 round(fi/kk)]),:);
        end
    end
end
tceps = tceps_new;
