function zc = makeZappCurrent(Iapp, fmin, fmax, tmax, fs)
% source https://nl.mathworks.com/matlabcentral/answers/217746-implementing-a-sine-wave-with-linearly-changing-frequency
    t = 0:1/fs:tmax;
    f_in = linspace(fmin, fmax, length(t));
    phase_in = cumsum(f_in/fs);
    zc = Iapp*sin(2*pi*phase_in); % approximation
    