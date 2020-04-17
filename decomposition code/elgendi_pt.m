function [pwaves,pwave_blocks,twaves,twave_blocks, blocks] = elgendi_pt(Om,mbeats,w1,w2)
y = Om;
%[b_lp,a_lp] = butter(3,[0.5,10]/(1000/2));
%Om = filtfilt(b_lp,a_lp,Om);
RRI = diff(mbeats)';
RRI = [RRI(1),RRI];
cal = 1000/median(RRI);
%cal = 1;
%w1 = round(55/cal);
%w2 = round(110/cal);
w1 = round(w1/cal);
w2 = round(w2/cal);
for i = 1 : length(mbeats)
    a = max(mbeats(i)-round(83/cal),1);
    b = min(mbeats(i)+round(166/cal),length(Om));
    if a==1
        aa = Om(a);
    else
        aa = Om(a-1);
    end
    
    if b == length(Om)
        bb = Om(end);
    else
        bb = Om(b+1);
    end
    
    Om(a:b) = min(aa,bb);
end
MApeak = medfilt2(Om', [1, w1]);
%MApeak = smooth(MApeak,10,'loess')';
MApwave = medfilt2(Om', [1, w2]);
%MApwave = smooth(MApwave,10,'loess')';
block_ind = zeros(1,length(MApeak));

block_ind(MApeak > MApwave) = 1;
blocks = [];

for i = 1:length(block_ind)
    if block_ind(i) == 0
        
        if i ~=1 & block_ind(i-1) == 1
            if i~=2
                blocks = [blocks,i-1];
            else
                continue;
            end
        else
            continue;
        end
    elseif block_ind(i) == 1
        
        if i == 1
            if block_ind(i+1) == 0
                blocks = [blocks,i,i];
            else
                blocks = [blocks,i];
            end
            
        elseif i == length(block_ind)
            
            if block_ind(i-1) == 1
                blocks = [blocks,i];
            else
                blocks = [blocks,i,i];
            end
            
        elseif block_ind(i-1) == 0
            blocks = [blocks,i];
            
        else
            continue;
        end
    end
end

blocks =  reshape(blocks,[2,round(length(blocks)/2)])';
%block_dur = diff(blocks');
%blocks(block_dur<w1/2,:) = [];
pwaves = [];
twaves = [];
pwave_blocks = [];
twave_blocks = [];

for i = 1:length(mbeats)
    Ri = mbeats(i);
    if i == 1
        Ri_1 = Ri-RRI(i);
    else
        Ri_1 = mbeats(i-1);
    end
    RR = Ri-Ri_1;
    %if RR > 1000
    %    RR = median(RRI);
    %    Ri_1 = mbeats(i)-RR;
    %end
    blocks_i = blocks(blocks(:,1)>Ri_1 & blocks(:,1)<Ri,:);
    block_i_w = diff(blocks_i');
    
    pwaves_i = [];
    twaves_i = [];
    pwave_blocks_i = [];
    twave_blocks_i = [];
    PRmin = 55*RR/1000;
    PRmax = 470*RR/1000;
    RTmin = 110*(RR/1000);
    RTmax = 860*(RR/1000);
    for j = 1:length(block_i_w)
        block_w = block_i_w(j);
        %if block_w <w1/2
        %    continue;
        %end
        if   Ri- blocks_i(j,2)>=PRmin &&  Ri-blocks_i(j,1) <=PRmax
            pwave_blocks_i = [pwave_blocks_i;blocks_i(j,:)];
            [~,ind] = max(Om(blocks_i(j,1):blocks_i(j,2)));
            pwaves_i = [pwaves_i, blocks_i(j,1)+ind-1];
            
        elseif  blocks_i(j,1)-Ri_1 >= RTmin && blocks_i(j,2) - Ri_1<=RTmax
            twave_blocks_i = [twave_blocks_i;blocks_i(j,:)];
            [~,ind] = max(Om(blocks_i(j,1):blocks_i(j,2)));
            twaves_i = [twaves_i, blocks_i(j,1)+ind-1];
        end
    end
    [~,ind] = max(y(pwaves_i));
    pwaves = [pwaves,pwaves_i(ind)];
    pwave_blocks = [pwave_blocks;pwave_blocks_i(ind,:)];
    
    [~,ind] = max(y(twaves_i));
    twaves = [twaves,twaves_i(ind)];
    twave_blocks = [twave_blocks;twave_blocks_i(ind,:)];
end
end