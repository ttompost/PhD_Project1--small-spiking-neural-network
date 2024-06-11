function a = howManyAsterisks(p)
    if p > 0.05
        a=0;
    elseif p <= 0.05 && p > 0.01
        a=1;
    elseif p <= 0.01 && p > 0.001
        a=2;
    elseif p <= 0.001 && p > 0.0001
        a=3;
    elseif p <= 0.0001
        a=4;
    else
        a=NaN;
    end
end
