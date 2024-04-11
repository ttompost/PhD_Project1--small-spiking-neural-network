function ConvTrace = convolve_whisking_signal_with_kernels(WhiskerStruct, KernelStruct, savename, plotyn)

if nargin == 3
    plotyn = 0;
end

[Ndimw, Ntrace] = size(WhiskerStruct.Recording);
[Nkernel, Ndimk] = size(KernelStruct.Kernels);
binsize_kernels = KernelStruct.kerneltime(2)-KernelStruct.kerneltime(1);
%% Check
if Nkernel < Ndimk
    ip = input([num2str(Nkernel) ' kernels with ' num2str(Ndimk) ' dimensions each. Is this ok, or transpose? (ok/tp)'], 's');
    if strcmp(ip, 'tp')
        KernelStruct.Kernels = KernelStruct.Kernels';
        [Nkernel, Ndimk] = size(KernelStruct.Kernels);
    end
end
        

if ~(Ndimk == Ndimw)
    error('Kernels should have the same dimensions as recordings')
end

for nt = 1:Ntrace
    for nd = 1:Ndimw
        if ~(length(WhiskerStruct.Recording{1,nt}) == length(WhiskerStruct.Recording{nd,nt}))
            error(['Error in trace ' num2str(nt) ', dimension ' num2str(nd) ': Make sure the length of the recording is the same for all dimensions'])
        end
    end
end

for nd = 1:Ndimk
    for nk = 1:Nkernel
        if ~(length(KernelStruct.Kernels{nk,nd}) == length(KernelStruct.kerneltime))
            error(['Error in kernel ' num2str(nk) ': Make sure all kernels have the same length as the kerneltime array'])
        end
    end
end


if WhiskerStruct.binsize == binsize_kernels
    disp('Sampling rate kernels the same as whisker recordings')
else
    if WhiskerStruct.binsize > binsize_kernels
        disp('Sampling rate kernels higher than whisker recordings: upsample whisker recordings')
    elseif WhiskerStruct.binsize < binsize_kernels
        disp('Sampling rate kernels lower than whisker recordings: downsample whisker recordings')
    end
    for nd = 1:Ndimw
        for nt = 1:Ntrace
            f1 = 1/binsize_kernels;
            f2 = 1/WhiskerStruct.binsize; 
            while ~((rem(f1,1)==0) && (rem(f2,1)==0))
                f1 = f1*10;
                f2 = f2*10;
            end
            WhiskerStruct.Recording{nd,nt} = resample(WhiskerStruct.Recording{nd,nt},f1 ,f2);
        end
    end
    WhiskerStruct.binsize = binsize_kernels;
end

%% Convolve recordings with relevant kernel
la = sum(KernelStruct.kerneltime<0);
lc = sum(KernelStruct.kerneltime>0);
ConvTrace = cell(Nkernel, Ntrace);
for nk = 1:Nkernel
    for nt = 1:Ntrace
        ConvTrace{nk,nt} = nan*ones(Ndimw, length(WhiskerStruct.Recording{nd,nt}));
        for nd = 1:Ndimw  
            ConvTrace{nk,nt}(nd,:) = convolve_kernel_acausal( WhiskerStruct.Recording{nd,nt}, KernelStruct.Kernels{nk, nd}, la, lc);
        end
    end
end

%% Save
if ~isempty(savename)
    save(savename, 'ConvTrace');
end

%% plot to check
maxct = 0;
minct = 0;
colorvec = {'b','r'};
for nk = 1:Nkernel
    for nt = 1:Ntrace
        maxt = max(max(ConvTrace{nk,nt}));
        mint = min(min(ConvTrace{nk,nt}));
        if maxt>maxct
            maxct = maxt;
        end
        if mint<minct
            minct = mint;
        end
    end
end
if maxct<1
    maxct = (10^floor(log10(maxct)))*round(maxct/(10^floor(log10(maxct))));
    minct = -maxct;
else
    maxct = ceil(maxct);
    minct = floor(minct);
end
kk = minct:(maxct-minct)/100:maxct;
nkk = length(kk);

if plotyn
    % Kernels
    figure
    plot_thalamic_kernels(KernelStruct) 
    
end

%% Helper functions
function s = convolve_kernel_acausal( signal, filter, la, lc)
% Convolve signal with kernel with causal and acausal part, give back on same time
% scale.

% Function input:
%   signal(vector) = spike train or input signal
%   filter (vector) = filter
%   la = length acausal part filter (tbegin(negative):0-dt)
%   lc = length causal part filter (dt:tmax/dt)
%   NB lc+la+1 = length(filter) for filters with both causal and acausal
%   parts

% Function output: convolved signal s

    cf = 0;
    ca = 0;
    c0 = 0;

    [~, Ny] = size(signal);
    if Ny == 1
        signal = signal';
    end

    if ~(la+lc+1==length(filter))
        if la == 0 && length(filter)==lc
            disp('purely causal filter')
            cf = 1;
        elseif lc == 0 && length(filter)==la
            disp('purely acausal filter')
            ca = 1;
        elseif (la+lc == length(filter))
            disp('filter does not contain t=0 value, making small delay')
            c0 = 1;
        else
            error('give correct lengths of causal and acausal part filter')
        end
    end

    if cf
        s = conv(signal, [0 filter]);
        s = s(1:length(signal));
    elseif ca
        s = conv([signal zeros(1,la)], [filter 0], 'valid');
    elseif c0
        s = conv([zeros(1,lc) signal zeros(1,la-1)], filter, 'valid');
    else
        s = conv([zeros(1,lc) signal zeros(1,la)], filter, 'valid');
    end
end

function plot_thalamic_kernels(KernelStruct)
    [Nbx, Nby] = size(KernelStruct);
    Nbarrel = Nbx*Nby;
    nb = 0;
    for nbx = 1:Nbx
        for nby=1:Nby
            nb = nb+1;
            [Nkernel, Ndimk]  = size(KernelStruct.Kernels);
            for nd = 1:Ndimk
                subplot(Nbarrel,Ndimk, (nb-1)*Ndimk+nd)
                title(['Dimension ', num2str(nd) ', barrel ' num2str(nb)])
                hold all
                for nk = 1:Nkernel            
                    plot(-1*fliplr(KernelStruct.kerneltime), fliplr(KernelStruct.Kernels{nk, nd}))
                    xlabel('time before spike (ms)')
                    ylabel('stimulus amplitude')
                    grid on
                    box on
                end
            end
        end
    end
end

end