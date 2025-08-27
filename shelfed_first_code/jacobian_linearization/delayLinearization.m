function [realBase, imagBase] = delayLinearization(q, Period_code)

    realBase = [];
    imagBase = [];
    Band = [-1/(2*Period_code) : 1/(Period_code*5001) : 1/(2*Period_code)];
    PROP = 1./(pi^2*Band.^2*Period_code);
    SIN = sin(pi*Band.*Period_code);
    EXPCOEFF = 2*pi*Band;
    for l = 0 : q
        for step = -q : q
            COSD = cos(2*pi*Band.*(l - step)*(Period_code/q));
            Function = EXPCOEFF.*PROP.*SIN.^2.*COSD;
            imagBase(step + q + 1, l + 1) = trapz(Band, Function);   
        end
    end
    for l = 0 : q
        for step = -q : q
            SIND = sin(2*pi*Band.*(l - step)*(Period_code/q));
            Function = -EXPCOEFF.*PROP.*SIN.^2.*SIND;
            realBase(step + q + 1, l + 1) = trapz(Band, Function);   
        end
    end
end