function roVibHamiltonianMatrix = ...
    getRoVibHamiltonianMatrix(v, j, gasConst, varargin)

%% Setting Up Options

% Define initial options to false.
rotWaveShift = 0;

% If options are passed in what are they and if they are one of the
% accepted options attribute their value to a variable inside this function
while length(varargin)>=2
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
        case 'rotating_wave_shift'
            rotWaveShift = val;
        otherwise
            error(['getRoVibHamiltonianMatrix: unknown option ',arg])
    end
    varargin = varargin(3:end);
end

B_e = gasConst.B_e;
D = gasConst.D;
a_e = gasConst.a_e;
% v_e = sqrt(4*B_e^3/abs(D));
v_ex_e = gasConst.v_ex_e;
v_e = gasConst.v_e;

% finding where the vibrational and rotational levels are the same
v = diag(v);
j = diag(j);

% Creating hamiltonian matrix

roVibHamiltonianMatrix = diag(v_e .* (v + 1/2) ...
    - v_ex_e .* (v + 1/2).^2 ...
    + B_e .* j .* (j+1) ...
    - a_e .* (v+1/2).*j.*(j+1) ...
    - D .* j.^2 .* (j+1).^2 - v.*rotWaveShift);
