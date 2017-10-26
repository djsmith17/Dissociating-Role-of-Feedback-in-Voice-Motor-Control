function f0 = dfCalcf0NFFT(x, fs)
NFFT = pow2(nextpow2(length(x)/4));

[Pmic, f] = pwelch(x, [], [], NFFT, fs, 'onesided');
    
[~, ind] = max(Pmic);
f0 = f(ind);
end