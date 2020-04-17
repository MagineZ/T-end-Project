function err1 = CompareWith2Cardios(Tend0, Tendindq1c)
% Compare automatically located T-wave ends with Tendindq1c, Tendindq2c


nTend = length(Tendindq1c);
Tcmp1 = zeros(nTend,1);
for ks = 1:nTend
    [dum, kc] = min(abs(Tend0-Tendindq1c(ks)));
    Tq1d0 = Tend0(kc);
    
    Tcmp1(ks) = Tq1d0;
    
    
end

err1 = Tcmp1-Tendindq1c;


