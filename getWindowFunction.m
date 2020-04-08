function [windowFunction2D, windowFunction1D] = getWindowFunction(windowType, windowParams)

switch windowType
    case 'hamming'
        windowFunction1D = @(t1, t1_max) (0.5 + 0.5 .* cos(pi.*t1./t1_max));
        windowFunction2D = @(t1, t3, t1_max, t3_max) (0.5 + 0.5 .* cos(pi.*t1./t1_max)) .* ( 0.5 + 0.5 .* cos(pi.*t3./t3_max));
    
    case 'shifted-sine-bell'
        t1_0 = windowParams.t1_0;
        t3_0 = windowParams.t3_0;
        
        windowFunction1D = @(t1, t1_max) sin(pi.*(t1 + t1_0)./(t1_max+t1_0));
        windowFunction2D = @(t1, t3, t1_max, t3_max) sin(pi.*(t1 + t1_0)./(t1_max+t1_0))...
            .*sin(pi.*(t3 + t3_0)./(t3_max + t3_0));
    otherwise
        error('Window type is unknown.');
end