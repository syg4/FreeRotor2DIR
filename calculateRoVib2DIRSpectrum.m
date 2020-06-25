params_new.path = fileparts(mfilename('fullpath'));

load(fullfile(params_new.path, 'Data', '2DIR', '076.mat'));         % 1: Load and process 2DIR data
r1 = [2300 2380];
r3 = [data(1).w3(1) data(1).w3(end)];

expData = cropData(data, r1, r3);                                   % crop experimental data to desired frequency window
clear data
%%
baseline_range = [3985 3990];
range1 = [2295 2400];
peakThreshold = 0.2;
params_new.symmetry = 'even';

params_new.dataPathFTIR = fullfile(params_new.path, 'Data', 'FTIR', 'CO2_THF_950um_FlowCell.CSV');      % 2: Load FTIR data and define parameters

[~, ~, params_new.gasConst] = getEmpericalGasPhaseHamiltonian(params_new.dataPathFTIR, range1,...       
    baseline_range, peakThreshold, params_new.symmetry);

clear baseline_range range1 peakThreshold
close([figure(100) figure(200)])
%%

params_new.flag_t1_parallel_pool = false; % Do you want t1 time dependence to be calculated in parallel
params_new.flag_t1_gpu = false; % Do you want t1 time-dependence to be calculated using gpu
params_new.flag_t3_parallel_pool = false; % Do you want t3 time dependence to be calculated in parallel
params_new.flag_t3_gpu = false;% Do you want t3 time-dependence to be calculated using gpu

params_new.maxVibLvl = 2;                           % define maximum vibrational and rotational levels obtainable in model
params_new.maxRotLvl = 5;

params_new.windowType = 'hamming';
params_new.windowParams = struct();

% J = 5     T = 3
% J = 7     T = 12
% J = 15    T = 50
% J = 25    T = 75
% J = 45    T = 273

params_new.temperature = 3;                   % define temperature for simulation

params_new.rotatingWaveShift = 2250;            % apply rotating frame to achieve desired sampling frequency(?)
params_new.w0 = 2350;

params_new.lineshapeForm = '1fast';
params_new.lineshapeParams = struct('T2', 100000);
%%
% params_new.t1s = 0:133.3347:30000;
% params_new.t3s = 0:132.7897:14000;
% params_new.t3s = params_new.t3s';
dw1 = expData(1).w1(2) - expData(1).w1(1);
w1 = 2225:dw1:2474;

dw3 = expData(1).w3(2) - expData(1).w3(1);
w3 = 2225:dw3:2475;
% %% 
params_new.t1s = fftTimeAxis(w1, 'time_units', 'fs');
params_new.t3s = fftTimeAxis(w3, 'time_units', 'fs')';

params_new.t2s = [200]; % 1000 10000 60000 100000];
params_new.t_J = 6000;

%%
if params_new.flag_t1_parallel_pool || params_new.flag_t3_parallel_pool
    if isempty(gcp('nocreate'))
        prpl = parpool;
    else
        fprintf('Parallel Pool Already Running.\nGetting Current Parallel Pool ... ');
        prpl = gcp;
        fprintf('Done.\n\n')
    end
end
% 
%%
try
    load(fullfile(params_new.path, 'Inputs', sprintf('simWorkspace_%iJ_%iK.mat',...     
        [params_new.maxRotLvl, params_new.temperature])));                              % 4: Begin calculating simulation matrices
    
    if ~isequal(params_new, params)
        setupSimWorkspace(params_new);      % 3: Build Hamiltonian (move to setupSimWorkspace.m)
        load(fullfile(params_new.path, 'Inputs', sprintf('simWorkspace_%iJ_%iK.mat',...
        [params_new.maxRotLvl, params_new.temperature])));
    end
catch
    fprintf('\nWorkspace not found\nGenerating workspace ...\n\n')
    setupSimWorkspace(params_new);
    load(fullfile(params_new.path, 'Inputs', sprintf('simWorkspace_%iJ_%iK.mat',...
        [params_new.maxRotLvl, params_new.temperature])));
    fprintf('Generating workspace ... Done\n')
end

clear params_new

n_zp1 = 1*numel(params.t1s); % Set amount of zero-padding in the t1 dimension
n_zp3 = 1*numel(params.t3s); % Set amount of zero-padding in the t3 dimension

%%
if params.flag_t1_gpu % If using gpu for t1 calculations convert these matrices to gpu arrays
    densityMatrix = gpuArray(densityMatrix);
    transitionDipoleMatrix = gpuArray(transitionDipoleMatrix);
    block = gpuArray(block);
end


%%
data = struct('t1s', [], 't3s', [], 't2', 0, 'w1', [], 'w3', [], 'R', [],...        % Store data in a structure so it's nice and organized
    'J', [], 'zeropad1', 0, 'zeropad3', 0, 'R_nr_raw', [], 'R_r_raw', [], ...
    'J_nr_raw', [], 'J_r_raw', [], 'R_nr', [], 'R_r', [], 'J_nr', [], 'J_r', [],...
    'lineshape1D', [], 'lineshape2D', [], 'lineshapeForm', params.lineshapeForm,...
    'lineshapeParams', params.lineshapeParams, 'window1D', [], 'window2D', [],...
    'windowType', params.windowType, 'windowParams', params.windowParams);
%%
for ii = 1:length(params.t2s)           % 5: For i=1:length(t2)
    
%     t2TimeProp = getT2TimePropagation(t2s(ii), t_J, maxVibLvl, maxRotLvl, temperature, 'even');

    t2TimeProp = 1;         % 6: Calculate t2 time propogation

    %%
    fprintf('Calculating Lineshape ... ');
    tic
    [lineshape2D, lineshape1D] = getLineshape(params.t1s, params.t2s(ii), params.t3s, ...       % 7: Calculate lineshape functions
        params.lineshapeForm, params.lineshapeParams);
    lineshapeTime = toc;
    fprintf('Done.\nTime: %f s\n\n', lineshapeTime);
    clear lineshapeTime
    
    if params.flag_t1_gpu % If using gpu for t1 calculations convert these matrices to gpu arrays
        lineshape1D = gpuArray(lineshape1D);
        t2TimeProp = gpuArray(t2TimeProp);
    end
    if params.flag_t3_gpu % If using gpu for t3 calculations convert these matrices to gpu arrays
        lineshape2D.Rr = gpuArray(lineshape2D.Rr);
        lineshape2D.Rnr = gpuArray(lineshape2D.Rnr);
    end

    %% Calculate The First- and Third-Order Response

    % Calculate the non-rephasing response
    fprintf('\nCalculating Non-Rephasing Response\n\n')             % 8: Calculate non-rephasing response
    [Rnr, Jnr] = getNonRephasingResponse(densityMatrix, transitionDipoleMatrix,...
        t1TimeProp, t3TimeProp, t2TimeProp, lineshape1D, lineshape2D.Rnr, block,...
        'flag_t1_parallel_pool', params.flag_t1_parallel_pool,...
        'flag_t1_gpu', params.flag_t1_gpu, ...
        'flag_t3_parallel_pool', params.flag_t3_parallel_pool,...
        'flag_t3_gpu', params.flag_t3_gpu);
    %%
    fprintf('\nCalculating Rephasing Response\n\n')                 % 9: Calculate the rephasing response
    [Rr, Jr] = getRephasingResponse(densityMatrix, transitionDipoleMatrix,...
        t1TimeProp, t3TimeProp, t2TimeProp, lineshape1D, lineshape2D.Rr, block,...
        'flag_t1_parallel_pool', params.flag_t1_parallel_pool,...
        'flag_t1_gpu', params.flag_t1_gpu, ...
        'flag_t3_parallel_pool', params.flag_t3_parallel_pool,...
        'flag_t3_gpu', params.flag_t3_gpu);
    %%
    data(ii).t1s = params.t1s;              % 10: Set up for and perform Fourier transforms
    data(ii).t3s = params.t3s;
    data(ii).t2 = params.t2s(ii);
    data(ii).zeropad1 = n_zp1;
    data(ii).zeropad3 = n_zp3;
    data(ii).w1 = fftFreqAxis(params.t1s, 'time_units', 'fs', 'zeropad', n_zp1)+params.rotatingWaveShift;
    data(ii).w3 = fftFreqAxis(params.t3s, 'time_units', 'fs', 'zeropad', n_zp3)+params.rotatingWaveShift;
    data(ii).R_nr_raw = Rnr;
    data(ii).R_r_raw = Rr;
    data(ii).J_nr_raw = Jnr;
    data(ii).J_r_raw = Jr;
    data(ii).R_nr = fftshift(ifft2(window2D.*Rnr, n_zp3, n_zp1)); % Perform a 2D inverse FT on the non-rephasing third-order response
    data(ii).R_r = fftshift(fliplr(circshift(ifft2(window2D.*Rr, n_zp3, n_zp1),[0 -1]))); % Perform a 2D inverse FT on the rephasing third-order response and transform spectra into w1 vs w3
    data(ii).R =  real(data(ii).R_nr + data(ii).R_r);
    data(ii).J_nr = fftshift(sgrsifft(window1D.*Jnr, n_zp1)); % Do inverse Fourier transform and shift the transform
    data(ii).J_r = fftshift(fliplr(circshift(sgrsifft(window1D.*Jr, n_zp1),[0 -1]))); % Do inverse Fourier transform
    data(ii).J = real(data(ii).J_nr + data(ii).J_r);
    data(ii).lineshape1D = lineshape1D;
    data(ii).lineshape2D = lineshape2D;
    data(ii).lineshapeForm = params.lineshapeForm;
    data(ii).lineshapeParams = params.lineshapeParams;
    data(ii).window1D = window1D;
    data(ii).window2D = window2D;
    data(ii).windowType = params.windowType;
    data(ii).windowParams = params.windowParams;
end

%%
todaysDir = fullfile(params.path, 'Outputs', sprintf('%s', datestr(now,'yyyy-mm-dd')));     % 11: Save simulation data
if ~exist(todaysDir, 'dir')
    mkdir(todaysDir)
end

filename = fullfile(todaysDir, sprintf('%iJ_%iK_%itJ_out.mat', [params.maxRotLvl, params.temperature, params.t_J]));

fileExists = true;
fileNum = 1;
while fileExists
    if exist(strcat(filename), 'file')
        if fileNum == 1
            filename = strrep(filename, '.mat', sprintf('_%i.mat', fileNum));
        else
            filename = strrep(filename, sprintf('_%i.mat', fileNum-1), sprintf('_%i.mat', fileNum));
        end
        fileNum = fileNum+1;
    else
        fileExists = false;
        save(filename, 'data');
    end
end
%%
% r1 = [2300 2380];
% r3 = [2300 2380];

cpData = cropData(data, r1, r3);           
%%

f = prepareGlobalFitData(cpData);
dim = [1 1];
fig = freeRotorMatrixPlot2DIR(real(f), cpData(1).w1, cpData(1).w3, [cpData(:).t2]./1000, dim, 'fignum', 1, 'zlimit', 1);

%%
simData = data;

for ii = 1:length(simData)
    
    [W1, W3] = meshgrid(expData(ii).w1, expData(ii).w3');

    simData.R = interp2(simData.w1, simData.w3, simData.R, W1, W3);          % 12: Interpolate simulation data to experiment
    simData.w1 = expData(ii).w1;
    simData.w3 = expData(ii).w3;
end

figure(2)                                                                    % 13: Plot simulation data
my2dRoVibPlot(simData(1).w1, simData(1).w3, simData(1).R)