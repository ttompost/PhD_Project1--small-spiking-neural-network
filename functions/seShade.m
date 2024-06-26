function seShade(amatrix,alpha,acolor,F,smth, lst)
% usage: stdshading(amatrix,alpha,acolor,F,smth)
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23
if exist('lst','var')==0 || isempty(lst)
    lst='-'; 
end
if exist('acolor','var')==0 || isempty(acolor)
    acolor='r'; 
end
if exist('F','var')==0 || isempty(F)
    F=1:size(amatrix,2);
end
if exist('smth','var'); if isempty(smth); smth=1; end
else 
    smth=1; %no smoothing by default
end  
if ~size(F,1)
    F=F';
end
amean = nanmean(amatrix,1); %get man over first dimension
if smth > 1
    amean = boxFilter(nanmean(amatrix,1),smth); %use boxfilter to smooth data
end
ase = nanstd(amatrix,[],1) ./ sqrt(size(amatrix,1)); % to get SE shading

ase(isnan(ase)) = 0;
amean(isnan(amean)) = 0;

% astd = nanstd(amatrix,[],1)/sqrt(size(amatrix,1)); % to get sem shading
if exist('alpha','var')==0 || isempty(alpha) 
    fill([F fliplr(F)],[amean+ase fliplr(amean-ase)],acolor,'linestyle','none','marker','none');
    acolor='k';
else 
    fill([F fliplr(F)],[amean+ase fliplr(amean-ase)],acolor, 'FaceAlpha', alpha,'linestyle','none', 'marker','none');    
end
if ishold==0
    check=true; 
else 
    check=false;
end
hold on;plot(F,amean,'color',acolor,'linewidth',1.2,'linestyle',lst,'marker','none'); %% change color or linewidth to adjust mean line
if check
    hold off;
end
end
function dataOut = boxFilter(dataIn, fWidth)
% apply 1-D boxcar filter for smoothing
dataStart = cumsum(dataIn(1:fWidth-2),2);
dataStart = dataStart(1:2:end) ./ (1:2:(fWidth-2));
dataEnd = cumsum(dataIn(length(dataIn):-1:length(dataIn)-fWidth+3),2);
dataEnd = dataEnd(end:-2:1) ./ (fWidth-2:-2:1);
dataOut = conv(dataIn,ones(fWidth,1)/fWidth,'full');
dataOut = [dataStart,dataOut(fWidth:end-fWidth+1),dataEnd];
end
