function [Nc, Nb, Ncu, CF] = calccofac(eventtimes1, eventtimes2, rate2, precision)
    N1=length(eventtimes1);  
    N2=length(eventtimes2);  

    coinc=[];
    Ncoinc=[];
    Nbar=[];
    Ncurl=[];
    CoFac=[];    
    for n=1:N1
        clear a
        a=find(eventtimes2>=eventtimes1(n)-precision&eventtimes2<=eventtimes1(n)+precision);
        if isempty(a)
        elseif length(a)==1
                coinc=[coinc eventtimes1(n)];
        elseif length(a)>1
                disp('two coincident spikes, make bin smaller')
                %continue
                Nc=[];
                Nb=[];
                Ncu=[];
                CF=[];
                return
        end
    end
    Ncoinc=length(coinc);
    Nbar=2*rate2*precision*N1;
    Ncurl=1-2*rate2*precision;
    CoFac=(Ncoinc-Nbar)/((1/2)*(N1+N2))*(1/Ncurl);
    Nc=Ncoinc;
    Nb=Nbar;
    Ncu=Ncurl;
    CF=CoFac;
    
return