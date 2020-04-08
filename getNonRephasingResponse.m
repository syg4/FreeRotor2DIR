% getNonRephasingResponse: This function generates all non-rephasing first-
%       and third-order response (R_4, R_5, R_6) needed to calculate a 
%       one-dimensional and two-dimensional infrared spectrum.
%
%   *Inputs*:
%       *densityMatrix*: An NxN matrix determined to represent your system at
%                       t = -inf.
%
%       *transitionDipoleMatrix*: An NxN matrix determined to represent the
%                               the allowed transition of your system.
%
%       *t1TimeProp*: A 1xlength(t_1) cell array where each cell is an NxN
%                   matrix describing the time evolution of the system at
%                   its corresponding t_1 time.
%
%       *t2TimeProp*: An NxN matrix determined to represent additional
%                   dynamics/kinentics during t_2 that are not accounted
%                   for in the lineshape function.
%
%       *t3TimeProp*: A 1xlength(t_3) cell array where each cell is an NxN
%                   matrix describing the time evolution of the system at
%                   its corresponding t_3 time.
%
%       *lineshape1D*: A 1xlength(t_1) array with each element
%                   corresponding to the value of the lineshape function at
%                   the elements corresponding t_1 time.
%
%       *lineshape2D*: A length(t_3)xlength(t_1) matrix with each element
%                   corresponding to the value of the 2D lineshape for the
%                   non-rephasing response at the elements corresponing
%                   t_1 and t_3 times.
%
%       *block*: An NxN block diagonal matrix with elementss equal to 1 when
%               the corresponding row and column of the element have equal
%               vibrational quantum numbers. All other values are equal to
%               zero.
%
%   *Options*:
%       *flag_t1_parallel_pool*: This is boolean stating whether or not the
%                   t_1 time dependence of the reponse calculation should be
%                   conducted using parallel computing methods. Check the
%                   availability of your computer to conduct parallel pool
%                   calculations before setting this option to true.
%
%       *flag_t1_gpu*: This is boolean stating whether or not the t_1 time
%                   dependence of the response calculation should be
%                   conducted using the GPU instead of the CPU. Check your
%                   computers ability to conduct GPU computing before
%                   setting this option to true.
%
%       *flag_t3_parallel_pool*: This is boolean stating whether or not the
%                   t_3 time dependence of the reponse calculation should be
%                   conducted using parallel computing methods. Check the
%                   availability of your computer to conduct parallel pool
%                   calculations before setting this option to true.
%
%       *flag_t3_gpu*: This is boolean stating whether or not the t_3 time
%                   dependence of the response calculation should be
%                   conducted using the GPU instead of the CPU. Check your
%                   computers ability to conduct GPU computing before
%                   setting this option to true.
%
%               Options are to be passed into this function using a two
%               input system. This first variable should be a string of the
%               option name (i.e. 'flag_t3_gpu'). The second variable
%               should be the value of the option (i.e. true or false). If
%               an option is not passed into the function the option will
%               default to false.
%
%   *Outputs*:
%       *Rnr*: This is a length(t_3)xlength(t_1) matrix where each element
%               is equal to the non-rephasing third-order response for the
%               elements corresponding t_1 and t_3 values.
%
%       *Jnr*: This is a 1xlength(t_1) array where each element is equsl to
%               the non-rephasing first-order response for the elements
%               corresponding t_1 value.

function [Rnr, Jnr] = getNonRephasingResponse(densityMatrix, transitionDipoleMatrix,...
    t1TimeProp, t3TimeProp, t2TimeProp, lineshape1D, lineshape2D, block, varargin)
%% Setting Up Options

% Define initial options to false.
flag_t1_parallel_pool = false;
flag_t1_gpu = false;
flag_t3_parallel_pool = false;
flag_t3_gpu = false;

% If options are passed in what are they and if they are one of the
% accepted options attribute their value to a variable inside this function
while length(varargin)>=2
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
        case 'flag_t1_parallel_pool'
            flag_t1_parallel_pool = val;
        case 'flag_t1_gpu'
            flag_t1_gpu = val;
        case 'flag_t3_parallel_pool'
            flag_t3_parallel_pool = val;
        case 'flag_t3_gpu'
            flag_t3_gpu = val;
        otherwise
            error(['my2dPlot: unknown option ',arg])
    end
    varargin = varargin(3:end);
end

%% Calculate The Non-Rephasing First-Order Response

% Calculate the density matrix evolution after the first pulse interaction
J_456 = transitionDipoleMatrix*densityMatrix; 
clear densityMatrix % clear the density matrix. We don't need it anymore and it's taking up space.

t1Length = numel(t1TimeProp); % Determine the length of t1
J_46 = cell(1, t1Length); % Define a cell array to hold matracies for J4 and J6 which remain the same for the first order response.
J_5 = cell(1, t1Length); % Define a cell array to hold matracies for J5
Jnr = zeros(1, t1Length); % Define an array to hold non-rephasing first-order response 

tElapsedTot = tic; % Start timer to track total elapsed time of all calculations
tElapsed1 = tic; % Start timer to track total time to calculate first-order response
fprintf('\nCalculating t1 Interactions ...\n\n');
if flag_t1_parallel_pool % If t1 evolution is to be calculated in parallel do this
    parfor ii = 1:t1Length % Loop over all t1 time propagation
        tLoop1 = tic; % Start timer to track how long each loop takes
        % Finish first-order non-rephasing response calculation for J4 and 6
        J_46(ii) = {t1TimeProp{ii}*J_456*conj(t1TimeProp{ii})*transitionDipoleMatrix};
        % Finish first-order non-rephasing response calculation for J5
        J_5(ii) = {transitionDipoleMatrix*t1TimeProp{ii}*J_456*conj(t1TimeProp{ii})}; 
        
        if flag_t1_gpu % If this calculation was being done using the GPU do this
            
            % Calculate the trace of the density matrix after the first
            % two pulse interaction and multiply it by the lineshape then 
            % gather it from the gpu back into RAM
            Jnr(ii) = gather(lineshape1D(ii).*(trace(J_46{ii})+trace(J_5{ii})));
            
        else
            % Calculate the trace of the density matrix after the first
            % two pulse interaction and multiply it by the lineshape
            Jnr(ii) = lineshape1D(ii).*(trace(J_46{ii})+trace(J_5{ii}));
        end
         % Print the percent progress elapsed time to calculate first-order
         % response and the time to calculate this singular loop
        fprintf('Percent Complete: %0.2f %%\t Elapsed Time For Section: %0.2f s\t Time For Loop: %0.2f s\n',...
            ii/t1Length*100, toc(tElapsed1), toc(tLoop1));
    end
else % If t1 evolution is not to be calculated in parallel do this
    % This loop is the same as the above loop except it is using a for loop
    % instead of a parfor loop to conduct the calculation linearly. See
    % above loop for comments on the interworkings of this loop.
    for ii = 1:t1Length
        tLoop1 = tic;
        J_46(ii) = {t1TimeProp{ii}*J_456*conj(t1TimeProp{ii})*transitionDipoleMatrix};
        J_5(ii) = {transitionDipoleMatrix*t1TimeProp{ii}*J_456*conj(t1TimeProp{ii})};
        if flag_t1_gpu
            Jnr(ii) = gather(lineshape1D(ii).*(trace(J_46{ii})+trace(J_5{ii})));
        else
            Jnr(ii) = lineshape1D(ii).*(trace(J_46{ii})+trace(J_5{ii}));
        end
        fprintf('Percent Complete: %0.2f %%\t Elapsed Time For Section: %0.2f s\t Time For Loop: %0.2f s\n',...
            ii/t1Length*100, toc(tElapsed1), toc(tLoop1));
    end
end

% Print total elapsed time of simulation at the end of the first-order
% response calculation
fprintf('\n=============================================\n');
fprintf('\t\tTotal Elapsed Time: %0.2f s\n', toc(tElapsedTot)); 
fprintf('\n=============================================\n\n');

clear t1TimeProp lineshape1D % Clear these! We won't be needing them anymore...
%% Calculate t_2 Time-Dependence

% Setting up cell arrays to hold matrices for each time point and each
% response pathway
rho4 = cell(1, t1Length);
rho5 = cell(1, t1Length);
rho6 = cell(1, t1Length);

fprintf('\nCalculating t2 Interactions...\n\n');
tElapsed2 = tic; % Start time to track elapsed time to calculate t_2 interactions
if flag_t1_parallel_pool % If t1 was calculated in parallel calculate t2 in parallel also
    parfor ii = 1:t1Length % Start a parallel loop over first-order density matrices
        tLoop2 = tic; % Start a timer to track how long each loop takes
        
        % Multiply each density matrix by a block diagonal matrix to eliminate
        % the vibrational coherences during t2
        rho46 = block.*J_46{ii}; 
        rho5t2 = block.*J_5{ii};
        
        % Apply t2 interaction matrices
        rho46 = t2TimeProp*rho46*conj(t2TimeProp);
        rho5t2 = t2TimeProp*rho5t2*conj(t2TimeProp);

        % Again apply block diagonal matracies to eliminate any vibrational
        % coherences caused by t2 time propagation
        rho46 = block.*rho46;
        rho5t2 = block.*rho5t2;

        % Depending on the options for t1 and t3 gpu calculations transform
        % variables into gpuArrays from regular arrays or vice versa.
        if (flag_t1_gpu && flag_t3_gpu) || (~flag_t1_gpu && ~flag_t3_gpu)
            rho4(ii) = {rho46*transitionDipoleMatrix};
            rho5(ii) = {transitionDipoleMatrix*rho5t2};
            rho6(ii) = {transitionDipoleMatrix*rho46};
        elseif ~flag_t1_gpu && flag_t3_gpu
            rho4(ii) = {gpuArray(rho46*transitionDipoleMatrix)};
            rho5(ii) = {gpuArray(transitionDipoleMatrix*rho5t2)};
            rho6(ii) = {gpuArray(transitionDipoleMatrix*rho46)};
        elseif flag_t1_gpu && ~flag_t3_gpu
            rho4(ii) = {gather(rho46*transitionDipoleMatrix)};
            rho5(ii) = {gather(transitionDipoleMatrix*rho5t2)};
            rho6(ii) = {gather(transitionDipoleMatrix*rho46)};
        end
%         clear rho46 rho5t2 % Variables can't be cleared inside parallel
%         for loops

        % Print progress, loop time, and total time for this section
        fprintf('Percent Complete: %0.2f %%\t Elapsed Time: %0.2f s\t Time For Loop: %0.2f s\n',...
            ii/t1Length*100, toc(tElapsed2), toc(tLoop2));

    end
    
else
    % This loop operates identically to the one above but without
    % parallelization.
    for ii = 1:t1Length
        tLoop2 = tic;
        rho46 = block.*J_46{ii};
        rho5t2 = block.*J_5{ii};

        rho46 = t2TimeProp*rho46*conj(t2TimeProp);
        rho5t2 = t2TimeProp*rho5t2*conj(t2TimeProp);

        rho46 = block.*rho46;
        rho5t2 = block.*rho5t2;

        if (flag_t1_gpu && flag_t3_gpu) || (~flag_t1_gpu && ~flag_t3_gpu)
            rho4(ii) = {rho46*transitionDipoleMatrix};
            rho5(ii) = {transitionDipoleMatrix*rho5t2};
            rho6(ii) = {transitionDipoleMatrix*rho46};
        elseif ~flag_t1_gpu && flag_t3_gpu
            rho4(ii) = {gpuArray(rho46*transitionDipoleMatrix)};
            rho5(ii) = {gpuArray(transitionDipoleMatrix*rho5t2)};
            rho6(ii) = {gpuArray(transitionDipoleMatrix*rho46)};
        elseif flag_t1_gpu && ~flag_t3_gpu
            rho4(ii) = {gather(rho46*transitionDipoleMatrix)};
            rho5(ii) = {gather(transitionDipoleMatrix*rho5t2)};
            rho6(ii) = {gather(transitionDipoleMatrix*rho46)};
        end
        clear rho46 rho5t2
        fprintf('Percent Complete: %0.2f %%\t Elapsed Time: %0.2f s\t Time For Loop: %0.2f s\n',...
            ii/t1Length*100, toc(tElapsed2), toc(tLoop2));
    end
end

fprintf('\n=============================================\n');
fprintf('\t\tTotal Elapsed Time: %0.2f s\n', toc(tElapsedTot));
fprintf('\n=============================================\n\n');
%% Calculate Non-Rephasing Third-Order Response

t3Length = numel(t3TimeProp); % Determine length of t3
% Setup matrices to hold response values
R4 = zeros(t3Length, t1Length);
R5 = zeros(t3Length, t1Length);
R6 = zeros(t3Length, t1Length);

% Transforming the transition dipole matrix to allow for use of GPU or CPU
% depending on options
if ~flag_t3_gpu && flag_t1_gpu
    transitionDipoleMatrix = gather(transitionDipoleMatrix);
elseif flag_t3_gpu && ~flag_t1_gpu
    transitionDipoleMatrix = gpuArray(transitionDipoleMatrix);
end

fprintf('\nCalculating t3 Interactions...\n\n');
tElapsed = tic; % Start timer to track how long it takes to calculate t3 time dependence
for ii = 1:t1Length % Loop over each density matrix
    tLoop = tic; % Start time to track how long it takes to run each loop
    r4t2 = rho4{ii};
    r5t2 = rho5{ii};
    r6t2 = rho6{ii};
    if flag_t3_parallel_pool % If t3 evolution is to be calculated in parallel do this
        parfor jj = 1:t3Length % Start parallel loop over t3 time propgators
            
            % Calculate the t3 time dependence of each density matrix
            r4 = transitionDipoleMatrix*t3TimeProp{jj}*r4t2*conj(t3TimeProp{jj});
            r5 = transitionDipoleMatrix*t3TimeProp{jj}*r5t2*conj(t3TimeProp{jj});
            r6 = transitionDipoleMatrix*t3TimeProp{jj}*r6t2*conj(t3TimeProp{jj});
            
            if flag_t3_gpu % If calculations were completed on GPU gather them back into RAM
                % Take the trace of the density matrix and multiply it by the corresponding 2D linshape
                R4(jj,ii) = gather(lineshape2D(jj, ii).*trace(r4));
                R5(jj,ii) = gather(lineshape2D(jj, ii).*trace(r5));
                R6(jj,ii) = gather(lineshape2D(jj, ii).*trace(r6));
            else % Do the same but you don't have to gather output back into RAM
                R4(jj,ii) = lineshape2D(jj, ii).*trace(r4);
                R5(jj,ii) = lineshape2D(jj, ii).*trace(r5);
                R6(jj,ii) = lineshape2D(jj, ii).*trace(r6);
            end
        end        
    else % If not parallel:
        % This loop operates the same way as its parallel counterpart but in a linear fashion
        for jj = 1:t3Length
            r4 = transitionDipoleMatrix*t3TimeProp{jj}*r4t2*conj(t3TimeProp{jj});
            r5 = transitionDipoleMatrix*t3TimeProp{jj}*r5t2*conj(t3TimeProp{jj});
            r6 = transitionDipoleMatrix*t3TimeProp{jj}*r6t2*conj(t3TimeProp{jj});
            
            if flag_t3_gpu
                R4(jj,ii) = gather(lineshape2D(jj, ii).*trace(r4));
                R5(jj,ii) = gather(lineshape2D(jj, ii).*trace(r5));
                R6(jj,ii) = gather(lineshape2D(jj, ii).*trace(r6));
            else
                R4(jj,ii) = lineshape2D(jj, ii).*trace(r4);
                R5(jj,ii) = lineshape2D(jj, ii).*trace(r5);
                R6(jj,ii) = lineshape2D(jj, ii).*trace(r6);
            end
            clear r4 r5 r6
        end
    end
    clear r4t2 r5t2 r6t2
    % Print progress, loop time, and time to calculate t3 dependence
    fprintf('Percent Complete: %0.2f %%\t Elapsed Time: %0.2f s\t Time For Loop: %0.2f s\n',...
        ii/t1Length*100, toc(tElapsed), toc(tLoop)); 
end

% Print total time of simulation
fprintf('\n=============================================\n');
fprintf('\t\tTotal Elapsed Time: %0.2f s\n', toc(tElapsedTot));
fprintf('\n=============================================\n\n');

Rnr = R6-R4-R5; % Combine responses

% Devide the first row and column of response by 2
Rnr(:,1) = Rnr(:,1)./2;
Rnr(1,:) = Rnr(1,:)./2;
