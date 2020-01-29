function [lineshape2D, lineshape1D] = getLineshape(t1s, t2, t3s, ...
    lineshapeForm, lineshapeParams)

T1 = repmat(t1s, [length(t3s), 1]);
T3 = repmat(t3s, [1, length(t1s)]);


lineshape = chooseLineshapeFunction(lineshapeForm, lineshapeParams);
g =@(t) lineshape.g(t, lineshapeParams);

gT1 = g(t1s);
expGT1 = exp(-1*gT1);
% lineshape1D = num2cell(expGT1);
lineshape1D = expGT1;
%%

expGTot = exp(-g(T1) + g(t2) - g(T3) - g(T1+t2) - g(t2 + T3) + g(T1+t2+T3));

% lineshape2D = num2cell(expGTot);
lineshape2D = expGTot;