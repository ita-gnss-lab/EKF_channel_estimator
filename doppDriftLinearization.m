function [realDoppDriftSensibilityBase, imagDoppDriftSensibilityBase] = doppDriftLinearization(q, Period_code)

    Base = [];
    Band = [-1/(2*Period_code) : 1/(Period_code*5001) : 1/(2*Period_code)];
    PROP = 1./(pi^3*Band.^3*Period_code);
    MULT = -pi*Band.*Period_code;
    SIN = sin(pi*Band.*Period_code);
    COS = cos(pi*Band.*Period_code);
    TOP = (MULT.*COS.*SIN - SIN.^2);
    for l = 0 : q
        for step = -q : q
            COSD = cos(2*pi*Band.*(l - step)*(Period_code/q));
            Function = TOP.*PROP.*COSD;
            realDoppDriftSensibilityBase(step + q + 1, l + 1) = trapz(Band, Function);   
        end
    end
    for l = 0 : q
        for step = -q : q
            SIND = sin(2*pi*Band.*(l - step)*(Period_code/q));
            Function = TOP.*PROP.*SIND;
            imagDoppDriftSensibilityBase(step + q + 1, l + 1) = trapz(Band, Function);   
        end
    end
end

