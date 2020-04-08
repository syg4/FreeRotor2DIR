function setupSimWorkspace(params)

flag_t1_parallel_pool = params.flag_t1_parallel_pool; % Do you want t1 time dependence to be calculated in parallel
flag_t1_gpu = params.flag_t1_gpu; % Do you want t1 time-dependence to be calculated using gpu
flag_t3_parallel_pool = params.flag_t3_parallel_pool; % Do you want t3 time dependence to be calculated in parallel
flag_t3_gpu = params.flag_t3_gpu;% Do you want t3 time-dependence to be calculated using gpu

maxVibLvl = params.maxVibLvl;
maxRotLvl = params.maxRotLvl;

temperature = params.temperature;

t1s = params.t1s;
t3s = params.t3s;

symmetry = params.symmetry;
rotatingWaveShift = params.rotatingWaveShift;
gasConst = params.gasConst;

%%
% maxVibLvl = 2;
% maxRotLvl = 45;

[v1, v2] = getVibrationalStateMatrices(maxVibLvl, maxRotLvl);
[j1, j2] = getJStateMatrices(maxVibLvl, maxRotLvl);
[m1, m2] = getMStateMatrices(maxVibLvl, maxRotLvl);

% t1s = 0:100:30000;
% t3s = 0:100:14000;
% t3s = t3s';

% temperature = 273; %K 12
old_workspace = [];
%%
try
    old_workspace = load(sprintf(['C:\\Users\\kaigr\\OneDrive\\Documents\\Graduate Research\\'...
        'MATLAB\\FreeRotorGPUCPUHybrid\\Inputs\\simWorkspace_%iJ_%iK.mat'],...
        [maxRotLvl, temperature]));
    
    if ~isequal(old_workspace.params.maxRotLvl, maxRotLvl) || ~isequal(old_workspace.params.maxVibLvl, maxVibLvl)
        fprintf('Calculating Transition Dipole Matrices ... ');
        tic
        
        [transitionDipoleMatrix] = getRoVibTransDipoleMatrix(v1, v2, j1, j2, m1, m2, 'z');
        
        transitionDipoleMatrix = sparse(transitionDipoleMatrix);
        dipoleTime = toc;
        fprintf('Done.\n')
        fprintf('Time to make %i states for transition dipole matrices: %f s\n\n', ...
            3*(maxRotLvl+1)^2, dipoleTime);
        clear dipoleTime
    else
        transitionDipoleMatrix = old_workspace.transitionDipoleMatrix;
    end
catch
    fprintf('Calculating Transition Dipole Matrices ... ');
    tic
    
    [transitionDipoleMatrix] = getRoVibTransDipoleMatrix(v1, v2, j1, j2, m1, m2, 'z');
    
    transitionDipoleMatrix = sparse(transitionDipoleMatrix);
    dipoleTime = toc;
    fprintf('Done.\n')
    fprintf('Time to make %i states for transition dipole matrices: %f s\n\n', ...
        3*(maxRotLvl+1)^2, dipoleTime);
    clear dipoleTime
end
%%
if isempty(old_workspace)
    fprintf('Calculating Hamiltonain Matrix ... ');
    tic
    HamiltonianMatrix_rotWave = sparse(getRoVibHamiltonianMatrix(v1, j1, gasConst, 'rotating_wave_shift', rotatingWaveShift));
    HamTime = toc;
    fprintf('Done.\n');
    fprintf('Time to make %i states for Hamiltonian Matrix: %f s\n\n', (maxRotLvl+1)^2, HamTime);
else
    if (~isequal(old_workspace.params.maxRotLvl, maxRotLvl)...
            || ~isequal(old_workspace.params.maxVibLvl, maxVibLvl)...
            || ~isequal(old_workspace.params.gasConst, gasConst)...
            || ~isequal(old_workspace.params.rotatingWaveShift, rotatingWaveShift))
        fprintf('Calculating Hamiltonain Matrix ... ');
        tic
        HamiltonianMatrix_rotWave = sparse(getRoVibHamiltonianMatrix(v1, j1, gasConst, 'rotating_wave_shift', rotatingWaveShift));
        HamTime = toc;
        fprintf('Done.\n');
        fprintf('Time to make %i states for Hamiltonian Matrix: %f s\n\n', (maxRotLvl+1)^2, HamTime);
    else
        HamiltonianMatrix_rotWave = old_workspace.HamiltonianMatrix_rotWave;
    end
end

clear HamTime

%%
if isempty(old_workspace)
    roVibHamiltonianMatrix = sparse(getRoVibHamiltonianMatrix(v1, j1, gasConst));

    fprintf('Calculating Density Matrix ... ');
    tic
    densityMatrix = sparse(getRoVibDensityMatrix(v1, v2, j1, j2, m1, m2, ...
        roVibHamiltonianMatrix, temperature, symmetry));
    DenTime = toc;
    fprintf('Done.\n');
    fprintf('Time to make %i states for density matrix: %f s\n\n', (maxRotLvl+1)^2, DenTime);

    densityMatrix = densityMatrix.*1000;
    clear roVibHamiltonianMatrix DenTime
else
    if (~isequal(old_workspace.params.maxRotLvl, maxRotLvl) ...
            || ~isequal(old_workspace.params.maxVibLvl, maxVibLvl) ...
            || ~isequal(old_workspace.params.gasConst, gasConst) ...
            || ~isequal(old_workspace.params.temperature, temperature) ...
            || ~isequal(old_workspace.params.symmetry, symmetry))
        
        roVibHamiltonianMatrix = sparse(getRoVibHamiltonianMatrix(v1, j1, gasConst));
        
        fprintf('Calculating Density Matrix ... ');
        tic
        densityMatrix = sparse(getRoVibDensityMatrix(v1, v2, j1, j2, m1, m2, ...
            roVibHamiltonianMatrix, temperature, symmetry));
        DenTime = toc;
        fprintf('Done.\n');
        fprintf('Time to make %i states for density matrix: %f s\n\n', (maxRotLvl+1)^2, DenTime);
        
        densityMatrix = densityMatrix.*1000;
        clear roVibHamiltonianMatrix DenTime
    else
        densityMatrix = old_workspace.densityMatrix;
    end
end

%%
if isempty(old_workspace)
    
    block = eq(v1,v2);
    block = ones(length(block)).*block;
    block = sparse(block);
    
else
    if (~isequal(old_workspace.params.maxVibLvl, maxVibLvl) ...
            || ~isequal(old_workspace.params.maxRotLvl, maxRotLvl))
        
        block = eq(v1,v2);
        block = ones(length(block)).*block;
        block = sparse(block);
        
    else
        block = old_workspace.block;
    end
end
%%
clear symmetry DenTime v1 v2 j1 j2 m1 m2 gasConst

%%
if isempty(old_workspace)
    fprintf('Calculating t1 Time Propagators ... ');
    tic
    t1TimeProp = getTimePropagationCell(HamiltonianMatrix_rotWave, t1s,...
        'flag_parallel_pool', flag_t1_parallel_pool, 'flag_gpu', flag_t1_gpu);
    timePropTime = toc;
    fprintf('Done.\nTime: %f s\n\n', timePropTime);
else
    if (~isequal(old_workspace.HamiltonianMatrix_rotWave, HamiltonianMatrix_rotWave) ...
            || ~isequal(old_workspace.params.t1s, t1s) ...
            || ~isequal(old_workspace.params.flag_t1_gpu, flag_t1_gpu))
        fprintf('Calculating t1 Time Propagators ... ');
        tic
        t1TimeProp = getTimePropagationCell(HamiltonianMatrix_rotWave, t1s,...
            'flag_parallel_pool', flag_t1_parallel_pool, 'flag_gpu', flag_t1_gpu);
        timePropTime = toc;
        fprintf('Done.\nTime: %f s\n\n', timePropTime);
    else
        t1TimeProp = old_workspace.t1TimeProp;
    end
end

%%
if isempty(old_workspace)
    fprintf('Calculating t1 Time Propagators ... ');
    tic
    t3TimeProp = getTimePropagationCell(HamiltonianMatrix_rotWave, t3s,...
        'flag_parallel_pool', flag_t3_parallel_pool, 'flag_gpu', flag_t3_gpu);
    timePropTime = toc;
    fprintf('Done.\nTime: %f s\n\n', timePropTime);
else
    if (~isequal(old_workspace.Hamiltonian_rotWave, HamiltonianMatrix_rotWave) ...
            || ~isequal(old_workspace.params.t3s, t3s) ...
            || ~isequal(old_workspace.params.flag_t3_gpu, flag_t3_gpu))
        
        fprintf('Calculating t1 Time Propagators ... ');
        tic
        t3TimeProp = getTimePropagationCell(HamiltonianMatrix_rotWave, t3s,...
            'flag_parallel_pool', flag_t3_parallel_pool, 'flag_gpu', flag_t3_gpu);
        timePropTime = toc;
        fprintf('Done.\nTime: %f s\n\n', timePropTime);
        
    else
        
        t3TimeProp = old_workspace.t3TimeProp;
    end
end
%%

if isempty(old_workspace)
    [windowFunction2D, windowFunction1D] = getWindowFunction(params.windowType,...
        params.windowParams);

    window1D = windowFunction1D(t1s, t1s(end).*ones(1,length(t1s)));
    [T1s, T3s] = meshgrid(t1s, t3s);
    window2D = windowFunction2D(T1s, T3s, t1s(end).*ones(length(t3s),length(t1s)),...
        t3s(end).*ones(length(t3s), length(t1s)));
else
    if (~isequal(old_workspace.params.windowType, params.windowType) ...
            || ~isequal(old_workspace.params.windowParams, params.windowParams) ...
            || ~isequal(old_workspace.params.t1s, t1s) ...
            || ~isequal(old_workspace.params.t3s, t3s))
        [windowFunction2D, windowFunction1D] = getWindowFunction(params.windowType,...
            params.windowParams);
        
        window1D = windowFunction1D(t1s, t1s(end).*ones(1,length(t1s)));
        [T1s, T3s] = meshgrid(t1s, t3s);
        window2D = windowFunction2D(T1s, T3s, t1s(end).*ones(length(t3s),length(t1s)),...
            t3s(end).*ones(length(t3s), length(t1s)));
    else
        windowFunction1D = old_workspace.windowFunction1D;
        windowFunction2D = old_workspace.windowFunction2D;
        window1D = old_workspace.window1D;
        window2D = old_workspace.window2D;
    end
end
        
%%
clear T1s T3s
clear timePropTime old_workspace
clear flag_t1_parallel_pool flag_t1_gpu flag_t3_parallel_pool flag_t3_gpu
clear maxVibLvl maxRotLvl temperature t1s t3s symmetry rotatingWaveShift gasConst

save(sprintf('C:\\Users\\kaigr\\OneDrive\\Documents\\Graduate Research\\MATLAB\\FreeRotorGPUCPUHybrid\\Inputs\\simWorkspace_%iJ_%iK.mat', [params.maxRotLvl, params.temperature]))

