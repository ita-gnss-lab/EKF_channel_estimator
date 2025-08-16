function [realBase, imagBase] = ambiguityVector(q, Period_code)

    realBase = [];
    imagBase = [];
    Band = [-1/(2*Period_code) : 1/(Period_code*5001) : 1/(2*Period_code)];
    PROP = 1./(pi^2*Band.^2*Period_code);
    SIN = sin(pi*Band.*Period_code);
    for l = 0 : q
        for step = -q : q
            COSD = cos(2*pi*Band.*(l - step)*(Period_code/q));
            Function = PROP.*SIN.^2.*COSD;
            realBase(step + q + 1, l + 1) = trapz(Band, Function);   
        end
    end
    for l = 0 : q
        for step = -q : q
            SIND = sin(2*pi*Band.*(l - step)*(Period_code/q));
            Function = PROP.*SIN.^2.*SIND;
            imagBase(step + q + 1, l + 1) = trapz(Band, Function);   
        end
    end
end