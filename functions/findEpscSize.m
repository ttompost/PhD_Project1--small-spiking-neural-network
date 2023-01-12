function epscs = findEpscSize(s)
    
    minDist = 3/0.01;
    minHeight = 0.2;
    minProm = 0.5;
    
    [epscs(:,1), epscs(:,2)] = findpeaks(s,'sortstr','desc','minpeakdistance',minDist,'minpeakheight',minHeight,'minpeakprominence',minProm);
end