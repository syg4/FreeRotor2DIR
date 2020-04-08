%% getTimePropagationCell: This function generates  the time evolution of
%                   the Hamiltonian matrix in accordance with the time-dependent
%                   Schodinger equation.
%
%   *Inputs*:
%       *HamiltonianMatrix*: An NxN matrix with elements equal to the
%               energy of the corresponding states of the density matrix.
%
%       *ts*: A 1xM array containing the desired time-steps to be used to
%               generate the time-dependent matrices
%
%   *Options*:
%       *flag_parallel_pool*: This is boolean stating whether or not the
%                   time dependence matracies should be calculated using
%                   parallel computing methods. Check the availability of
%                   your computer to conduct parallel pool calculations
%                   before setting this option to true.
%
%       *flag_gpu*: This is boolean stating whether or not the time
%                   dependence matracies should be defined as arrays on the
%                   gpu instead of the CPU. Check your computers ability to
%                   conduct GPU computing before setting this option to true.
%
%               Options are to be passed into this function using a two
%               input system. This first variable should be a string of the
%               option name (i.e. 'flag_gpu'). The second variable
%               should be the value of the option (i.e. true or false). If
%               an option is not passed into the function the option will
%               default to false.
%
%   *Outputs*
%       *timePropagationCell*: A 1xlength(ts) cell array. Each cell
%               contains an NxN matrix representing the time-evolution of
%               the Hamiltonian at the cell's corresponding time point.

function timePropagationCell = getTimePropagationCell(HamiltonianMatrix, ts, varargin)

%% Setting Up Options

% Define initial options to false.
flag_parallel_pool = false;
flag_gpu = false;

% If options are passed in what are they and if they are one of the
% accepted options attribute their value to a variable inside this function
while length(varargin)>=2
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
        case 'flag_parallel_pool'
            flag_parallel_pool = val;
        case 'flag_gpu'
            flag_gpu = val;
        otherwise
            error(['my2dPlot: unknown option ',arg])
    end
    varargin = varargin(3:end);
end
%% Defining Constants

h = 6.626E-34*5.034E22*10^15; % Define Planck's constant in [cm^-1*fs]
hbar = h/(2*pi); % Define hbar
i = sqrt(-1); % Define imaginary number

%% Calculate Time-Dependence Matraces

% Define a cell array to store the matrices
timePropagationCell = cell(1, length(ts));

if flag_parallel_pool % If the matraces are to be calculated in parallel do this
    parfor ii = 1:length(ts) % Loop over all ts in parallel
        
        % Calculate the the time-dependent matrix at t
        timePropagationMatrix = expm(-1.*i.*HamiltonianMatrix.*ts(ii)./hbar);
        % Convert matrix to a sparse matrix
        timePropagationMatrix = sparse(timePropagationMatrix);
        
        if flag_gpu % If the matrix should be stored on GPU do this
            % Convert matrix into GPU array
            timePropagationMatrix = gpuArray(timePropagationMatrix);
        end
        % Add matrix to cell array
        timePropagationCell(ii) = {timePropagationMatrix};
    end
else
    % This loop works the same way as the one above but loops in a linear
    % fashion.
    for ii = 1:length(ts)
        timePropagationMatrix = expm(-1.*i.*HamiltonianMatrix.*ts(ii)./hbar);
        timePropagationMatrix = sparse(timePropagationMatrix);
        
        if flag_gpu
            timePropagationMatrix = gpuArray(timePropagationMatrix);
        end
        
        timePropagationCell(ii) = {timePropagationMatrix};
    end
end
