%% getLineshape: This function generates the one-dimensional and
%       2-dimensional linshapes for calculating linear and non-linear
%       spectra.
%
%   *Inputs*
%       *t1s*: A 1xN array of desired time-steps through t1.
%
%       *t2*: A scalar holding the desired duration of t2.
%
%       *t3s*: A 1xM array of desired time-steps through t3.
%
%       *lineshapeForm*: A string holding one of the available lineshape
%               forms in @chooseLineshapFunction.
%
%       *linshapeParams*: A structure holding the parameters needed to
%               define the exact lineshape form. The structure elements
%               vary depending on the lineshape form selected and can be
%               better determine by viewing @chooseLineshapeFunctions
%               documentation.
%
%   *Outputs*
%       *lineshape2D*: An MxN matrix holding a value of the two-dimensional
%               lineshape in each element. The value of the lineshape
%               coincides with the same indeces of t3 and t1 time points
%               (e.g. lineshape2D(ii,jj) gives the lineshape at t3s(ii) and
%               t1s(jj)).
%
%       *lineshape1D*:  A 1xN array holding a value of the one-dimensional
%               lineshape in each element. The index of each lineshape
%               value corresponds to the same index of the t1s array (e.g.
%               lineshape1D(ii) gibes the lineshape at t1s(ii)).


function [lineshape2D, lineshape1D] = getLineshape(t1s, t2, t3s, ...
    lineshapeForm, lineshapeParams)
%% Setting Up Mess Grids Of Time
T1 = repmat(t1s, [length(t3s), 1]);
T3 = repmat(t3s, [1, length(t1s)]);

%% Defining The Lineshape Function
lineshape = chooseLineshapeFunction(lineshapeForm, lineshapeParams);
g =@(t) lineshape.g(t, lineshapeParams);

%% Calculating The One-Dimensional Lineshape
gT1 = g(t1s); % Plugging t1s into the lineshape function
expGT1 = exp(-1*gT1); % Raising e to the power of the negative lineshape function
lineshape1D = expGT1; % Defining the 1D lineshape return array

%% Calculating the Two-Dimensional Lineshape

% Calculate the 2D rephasing lineshape
expG_Rr = exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3));

% Calculate the 2D non-rephasing lineshape
expG_Rnr = exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3));

% Store both lineshapes in a structure for return
lineshape2D.Rr = expG_Rr;
lineshape2D.Rnr = expG_Rnr;