function [Freq, P1] = getFFT(s, fs)
% Sampling_period = 0.01; % in ms
% fs = 1/(Sampling_period/1000); % in Hz

    Size = size(s);
    Len = max(Size);
    padding = 2^nextpow2(Len);
    if Size(1) == 1
        s = s'; % reverse dimensions
    end
    Y = fft(s,padding,1); % pad column for better fft
    P2 = abs(Y/Len); % double-sided spectrum
    P1 = P2(1:padding/2+1,:); 
    P1(2:end-1,:) = 2*P1(2:end-1,:); % one-sided spectrum

    Freq = 0:fs/padding:(fs/2-fs/padding);