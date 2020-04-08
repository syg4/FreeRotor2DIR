% getRephasingResponse: This function generates all rephasing first-
%       and third-order response (R_1, R_2, R_3) needed to calculate a 
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
%                   rephasing response at the elements corresponing
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
%       *Rr*: This is a length(t_3)xlength(t_1) matrix where each element
%               is equal to the rephasing third-order response for the
%               elements corresponding t_1 and t_3 values.
%
%       *Jr*: This is a 1xlength(t_1) array where each element is equsl to
%               the rephasing first-order response for the elements
%               corresponding t_1 value.

function [Rr, Jr] = getRephasingResponse(densityMatrix, transitionDipoleMatrix,...
    t1TimeProp, t3TimeProp, t2TimeProp, lineshape1D, lineshape2D, block, varargin)
%% Setting Up Options

% Define initial options to equal false.
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
%% Calculate The Rephasing First-Order Response

% Calculate the density matrix evolution after the first pulse interaction
J_123 = densityMatrix*transitionDipoleMatrix;
clear densityMatrix % clear the density matrix. We don't need it anymore and it's taking up space.

t1Length = numel(t1TimeProp); % Determine the length of t1
J_13 = cell(1, t1Length); % Define a cell array to hold matracies for J1 and J3 which remain the same for the first order response.
J_2 = cell(1, t1Length); % Define a cell array to hold matracies for J2
Jr = zeros(1, t1Length); % Define an array to hold rephasing first-order response 

tElapsedTot = tic; % Start timer to track total elapsed time of all calculations
tElapsed1 = tic; % Start timer to track total time to calculate first-order response
if flag_t1_parallel_pool % If t1 evolution is to be calculated in parallel do this
    parfor ii = 1:t1Length % Loop over all t1 time propagation
        tLoop1 = tic; % Start timer to track how long each loop takes
        % Finish first-order rephasing response calculation for J1 and 3
        J_13(ii) = {transitionDipoleMatrix*t1TimeProp{ii}*J_123*conj(t1TimeProp{ii})};
        % Finish first-order non-rephasing response calculation for J2
        J_2(ii) = {t1TimeProp{ii}*J_123*conj(t1TimeProp{ii})*transitionDipoleMatrix};
        
        if flag_t1_gpu % If this calculation was being done using the GPU do this
            
            % Calculate the trace of the density matrix after the first
            % two pulse interaction and multiply it by the lineshape then 
            % gather it from the gpu back into RAM
            Jr(ii) = gather(lineshape1D(ii).*(trace(J_13{ii})+trace(J_2{ii})));
        else
            % Calculate the trace of the density matrix after the first
            % two pulse interaction and multiply it by the lineshape
            Jr(ii) = lineshape1D(ii).*(trace(J_13{ii})+trace(J_2{ii}));
        end
        % Print the percent progress elapsed time to calculate first-order
        % response and the time to calculate this singular loop
        fprintf('Percent Complete: %0.2f %%\t Elapsed Time For Section: %0.2f s\t Time For Loop: %0.2f s\n', ii/t1Length*100, toc(tElapsed1), toc(tLoop1));
    end
else % If t1 evolution is not to be calculated in parallel do this
    % This loop is the same as the above loop except it is using a for loop
    % instead of a parfor loop to conduct the calculation linearly. See
    % above loop for comments on the interworkings of this loop.
    for ii = 1:t1Length
        tLoop1 = tic;
        J_13(ii) = {transitionDipoleMatrix*t1TimeProp{ii}*J_123*conj(t1TimeProp{ii})};
        J_2(ii) = {t1TimeProp{ii}*J_123*conj(t1TimeProp{ii})*transitionDipoleMatrix};
        if flag_t1_gpu
            Jr(ii) = gather(lineshape1D(ii).*(trace(J_13{ii})+trace(J_2{ii})));
        else
            Jr(ii) = lineshape1D(ii).*(trace(J_13{ii})+trace(J_2{ii}));
        end
        fprintf('Percent Complete: %0.2f %%\t Elapsed Time For Section: %0.2f s\t Time For Loop: %0.2f s\n', ii/t1Length*100, toc(tElapsed1), toc(tLoop1));
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
rho1 = cell(1, t1Length);
rho2 = cell(1, t1Length);
rho3 = cell(1, t1Length);


tElapsed2 = tic; % Start time to track elapsed time to calculate t_2 interactions
if flag_t1_parallel_pool % If t1 was calculated in parallel calculate t2 in parallel also
    parfor ii = 1:t1Length % Start a parallel loop over first-order density matrices
        tLoop2 = tic; % Start a timer to track how long each loop takes
        
        % Multiply each density matrix by a block diagonal matrix to eliminate
        % the vibrational coherences during t2
        rho13 = block.*J_13{ii};
        rho2t2 = block.*J_2{ii};

        % Apply t2 interaction matrices
        rho13 = t2TimeProp*rho13*conj(t2TimeProp);
        rho2t2 = t2TimeProp*rho2t2*conj(t2TimeProp);

        % Again apply block diagonal matracies to eliminate any vibrational
        % coherences caused by t2 time propagation
        rho13 = block.*rho13;
        rho2t2 = block.*rho2t2;

        % Depending on the options for t1 and t3 gpu calculations transform
        % variables into gpuArrays from regular arrays or vice versa.
        if (flag_t1_gpu && flag_t3_gpu) || (~flag_t1_gpu && ~flag_t3_gpu)
            rho1(ii) = {rho13*transitionDipoleMatrix};
            rho2(ii) = {transitionDipoleMatrix*rho2t2};
            rho3(ii) = {transitionDipoleMatrix*rho13};
        elseif ~flag_t1_gpu && flag_t3_gpu
            rho1(ii) = {gpuArray(rho13*transitionDipoleMatrix)};
            rho2(ii) = {gpuArray(transitionDipoleMatrix*rho2t2)};
            rho3(ii) = {gpuArray(transitionDipoleMatrix*rho13)};
        elseif flag_t1_gpu && ~flag_t3_gpu
            rho1(ii) = {gather(rho13*transitionDipoleMatrix)};
            rho2(ii) = {gather(transitionDipoleMatrix*rho2t2)};
            rho3(ii) = {gather(transitionDipoleMatrix*rho13)};
        end
%         clear rho46 rho5t2 % Variables can't be cleared inside parallel
%         for loops

        % Print progress, loop time, and total time for this section
        fprintf('Percent Complete: %0.2f %%\t Elapsed Time: %0.2f s\t Time For Loop: %0.2f s\n', ii/t1Length*100, toc(tElapsed2), toc(tLoop2));

    end
    
else
    % This loop operates identically to the one above but without
    % parallelization.
    for ii = 1:t1Length
        tLoop2 = tic;
        rho13 = block.*J_13{ii};
        rho2t2 = block.*J_2{ii};

        rho13 = t2TimeProp*rho13*conj(t2TimeProp);
        rho2t2 = t2TimeProp*rho2t2*conj(t2TimeProp);

        rho13 = block.*rho13;
        rho2t2 = block.*rho2t2;

        if (flag_t1_gpu && flag_t3_gpu) || (~flag_t1_gpu && ~flag_t3_gpu)
            rho1(ii) = {rho13*transitionDipoleMatrix};
            rho2(ii) = {transitionDipoleMatrix*rho2t2};
            rho3(ii) = {transitionDipoleMatrix*rho13};
        elseif ~flag_t1_gpu && flag_t3_gpu
            rho1(ii) = {gpuArray(rho13*transitionDipoleMatrix)};
            rho2(ii) = {gpuArray(transitionDipoleMatrix*rho2t2)};
            rho3(ii) = {gpuArray(transitionDipoleMatrix*rho13)};
        elseif flag_t1_gpu && ~flag_t3_gpu
            rho1(ii) = {gather(rho13*transitionDipoleMatrix)};
            rho2(ii) = {gather(transitionDipoleMatrix*rho2t2)};
            rho3(ii) = {gather(transitionDipoleMatrix*rho13)};
        end
        clear rho46 rho5t2
        fprintf('Percent Complete: %0.2f %%\t Elapsed Time: %0.2f s\t Time For Loop: %0.2f s\n', ii/t1Length*100, toc(tElapsed2), toc(tLoop2));
    end
end

fprintf('\n=============================================\n');
fprintf('\t\tTotal Elapsed Time: %0.2f s\n', toc(tElapsedTot));
fprintf('\n=============================================\n\n');
%% Calculate Non-Rephasing Third-Order Response

t3Length = numel(t3TimeProp); % Determine length of t3
% Setup matrices to hold response values
R1 = zeros(t3Length, t1Length);
R2 = zeros(t3Length, t1Length);
R3 = zeros(t3Length, t1Length);

% Transforming the transition dipole matrix to allow for use of GPU or CPU
% depending on options
if ~flag_t3_gpu && flag_t1_gpu
    transitionDipoleMatrix = gather(transitionDipoleMatrix);
elseif flag_t3_gpu && ~flag_t1_gpu
    transitionDipoleMatrix = gpuArray(transitionDipoleMatrix);
end

tElapsed = tic; % Start timer to track how long it takes to calculate t3 time dependence
for ii = 1:t1Length % Loop over each density matrix
    tLoop = tic; % Start time to track how long it takes to run each loop
    r1t2 = rho1{ii};
    r2t2 = rho2{ii};
    r3t2 = rho3{ii};
    if flag_t3_parallel_pool % If t3 evolution is to be calculated in parallel do this
        parfor jj = 1:t3Length % Start parallel loop over t3 time propgators
            
            % Calculate the t3 time dependence of each density matrix
            r1 = transitionDipoleMatrix*t3TimeProp{jj}*r1t2*conj(t3TimeProp{jj});
            r2 = transitionDipoleMatrix*t3TimeProp{jj}*r2t2*conj(t3TimeProp{jj});
            r3 = transitionDipoleMatrix*t3TimeProp{jj}*r3t2*conj(t3TimeProp{jj});
            
            if flag_t3_gpu % If calculations were completed on GPU gather them back into RAM
                % Take the trace of the density matrix and multiply it by the corresponding 2D linshape
                R1(jj,ii) = gather(lineshape2D(jj, ii).*trace(r1));
                R2(jj,ii) = gather(lineshape2D(jj, ii).*trace(r2));
                R3(jj,ii) = gather(lineshape2D(jj, ii).*trace(r3));
            else % Do the same but you don't have to gather output back into RAM
                R1(jj,ii) = lineshape2D(jj, ii).*trace(r1);
                R2(jj,ii) = lineshape2D(jj, ii).*trace(r2);
                R3(jj,ii) = lineshape2D(jj, ii).*trace(r3);
            end
        end        
    else % If not parallel:
        % This loop operates the same way as its parallel counterpart but in a linear fashion
        for jj = 1:t3Length
            r1 = transitionDipoleMatrix*t3TimeProp{jj}*r1t2*conj(t3TimeProp{jj});
            r2 = transitionDipoleMatrix*t3TimeProp{jj}*r2t2*conj(t3TimeProp{jj});
            r3 = transitionDipoleMatrix*t3TimeProp{jj}*r3t2*conj(t3TimeProp{jj});
            
            if flag_t3_gpu
                R1(jj,ii) = gather(lineshape2D(jj, ii).*trace(r1));
                R2(jj,ii) = gather(lineshape2D(jj, ii).*trace(r2));
                R3(jj,ii) = gather(lineshape2D(jj, ii).*trace(r3));
            else
                R1(jj,ii) = lineshape2D(jj, ii).*trace(r1);
                R2(jj,ii) = lineshape2D(jj, ii).*trace(r2);
                R3(jj,ii) = lineshape2D(jj, ii).*trace(r3);
            end
            clear r1 r2 r3
        end
    end
    clear r1t2 r2t2 r3t2
    % Print progress, loop time, and time to calculate t3 dependence
    fprintf('Percent Complete: %0.2f %%\t Elapsed Time: %0.2f s\t Time For Loop: %0.2f s\n', ii/t1Length*100, toc(tElapsed), toc(tLoop));
end

% Print total time of simulation
fprintf('\n=============================================\n');
fprintf('\t\tTotal Elapsed Time: %0.2f s\n', toc(tElapsedTot));
fprintf('\n=============================================\n\n');

Rr = R3-R1-R2; % Combine responses

% Devide the first row and column of response by 2
Rr(:,1) = Rr(:,1)./2;
Rr(1,:) = Rr(1,:)./2;
