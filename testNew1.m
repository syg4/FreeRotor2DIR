% global h c

maxVibLvl = 2;
maxRotLvl = 7; %7

[v1, v2] = getVibrationalStateMatrices(maxVibLvl, maxRotLvl);
[j1, j2] = getJStateMatrices(maxVibLvl, maxRotLvl);
[m1, m2] = getMStateMatrices(maxVibLvl, maxRotLvl);

t1s = 0:200:55000;
t3s = 0:200:55000;
t3s = t3s';
t2 = 1;
t_J = 20000;

temperature = 12; %K 12


lineshapeForm = '1fast';
lineshapeParams = struct('T2',20000);

parapool_flag = true;
%%
fprintf('Calculating Transition Dipole Matrices ... ');
tic
[transDipoleMatrixZ] = getRoVibTransDipoleMatrix_Vectorized(v1, v2, j1, j2, m1, m2, 'z');

transDipoleMatrixZ = sparse(transDipoleMatrixZ);
dipoleTime = toc;
fprintf('Done.\n')
fprintf('Time to make %i states for transition dipole matrices: %f s\n\n', ...
    3*(maxRotLvl+1)^2, dipoleTime);
clear dipoleTime

%%

baseline_range = [3985 3990];
range1 = [2295 2400];
peakThreshold = 0.2;
symmetry = 'even';

dataPathFTIR = ['C:\Users\kaigr\OneDrive\Documents\Graduate Research\MATLAB'...
    '\FreeRotor2DIR\new\Data\CO2_THF_950um_FlowCell.CSV'];

[~, ~, gasConst] = getEmpericalGasPhaseHamiltonian(dataPathFTIR, range1,...
    baseline_range, peakThreshold, symmetry);

% momentOfInertia = h/(8*pi^2*c*gasConst.B_e);

fprintf('Calculating Hamiltonain Matrix ... ');
tic
roVibHamiltonianMatrix_rotWave = sparse(getRoVibHamiltonianMatrix_RotWave(v1, v2, j1, j2, m1, m2, gasConst));
HamTime = toc;
fprintf('Done.\n');
fprintf('Time to make %i states for Hamiltonian Matrix: %f s\n\n', (maxRotLvl+1)^2, HamTime);
% fprintf('Do the peak identifications and polynomial fits look good?\n\tIf so press any key to continue ... ')
% pause;
% fprintf('Continuing with computations!\n\n');
pause(1);
close([figure(100) figure(200)])
clear baseline_range HamTime dataPathFTIR peakThreshold range1 

%%
roVibHamiltonianMatrix = sparse(getRoVibHamiltonianMatrix_Vectorized(v1, v2, j1, j2, m1, m2, gasConst));

fprintf('Calculating Density Matrix ... ');
tic
roVibDensityMatrix = sparse(getRoVibDensityMatrix(v1, v2, j1, j2, m1, m2, ...
    roVibHamiltonianMatrix, temperature, symmetry));
DenTime = toc;
fprintf('Done.\n');
fprintf('Time to make %i states for density matrix: %f s\n\n', (maxRotLvl+1)^2, DenTime);

roVibDensityMatrix = roVibDensityMatrix.*1000;

%%
% momentOfInertia = ;
% temperature = 3;
t2TimeProp = getT2TimePropagation(t2, t_J, maxVibLvl, maxRotLvl, temperature, 'even');

t2TimeProp.Norm = 1;
t2TimeProp.Conj = 1;
%%
clear symmetry temperature DenTime v1 v2 j1 j2 m1 m2 maxRotLvl maxVibLvl gasConst

%%
if parapool_flag
    if isempty(gcp('nocreate'))
        prpl = parpool;
    else
        fprintf('Parallel Pool Already Running.\nGetting Current Parallel Pool ... ');
        prpl = gcp;
        fprintf('Done.\n\n')
    end
end
%%
fprintf('Calculating Time Propagators ... ');
tic
t1TimeProp = getTimePropagationCells(roVibHamiltonianMatrix_rotWave, t1s, parapool_flag);
t3TimeProp = getTimePropagationCells(roVibHamiltonianMatrix_rotWave, t3s, parapool_flag);
timePropTime = toc;
fprintf('Done.\nTime: %f s\n\n', timePropTime);
clear timePropTime
%%
fprintf('Calculating Lineshape ... ');
tic
[lineshape2D, lineshape1D] = getLineshape(t1s, t2, t3s, ...
    lineshapeForm, lineshapeParams);
lineshapeTime = toc;
fprintf('Done.\nTime: %f s\n\n', lineshapeTime);
clear lineshapeTime


%%
R_r_fun = {@getR1Cell, @getR2Cell, @getR3Cell};
R_nr_fun = {@getR4Cell, @getR5Cell, @getR6Cell};
R_r = cell(1, length(R_r_fun));
R_nr = cell(1, length(R_nr_fun));
J_r = cell(1, length(R_r_fun));
J_nr = cell(1, length(R_nr_fun));
R_r_time = zeros(1, length(R_r_fun));
R_nr_time = zeros(1, length(R_nr_fun));

if parapool_flag
    fprintf('Calculating Rephasing Response ... \n');
    tic
    parfor ii = 1:length(R_r)
        fprintf('Calculating R%i ... ', ii);
        [R_r_i, J_r_i] = R_r_fun{ii}(roVibDensityMatrix, transDipoleMatrixZ,...
            t1TimeProp, t3TimeProp, t2TimeProp, lineshape1D, lineshape2D);
        fprintf('Done.\n\n');
        R_r(ii) = {R_r_i};
        J_r(ii) = {J_r_i};    
    end
    R_rtime = toc;
    fprintf('Rephasing Response Caluclations Complete.\nTime: %0.2f s\n\n', R_rtime);
    clear ii R_rtime

    fprintf('Calculating Non-Rephasing Response ... \n');
    tic
    parfor ii = 1:length(R_nr)
        fprintf('Calculating R%i ... ', ii+3);
        [R_nr_i, J_nr_i] = R_nr_fun{ii}(roVibDensityMatrix, transDipoleMatrixZ,...
            t1TimeProp, t3TimeProp, t2TimeProp, lineshape1D, lineshape2D);
        fprintf('Done.\n\n');
        R_nr(ii) = {R_nr_i};
        J_nr(ii) = {J_nr_i};
    end
    R_nrtime = toc;
    fprintf('Non-Rephasing Response Caluclations Complete.\nTime: %0.2f s\n\n', R_nrtime);
    clear ii R_nrtime
    
else
    
    fprintf('Calculating Rephasing Response ... \n');
    tic
    for ii = 1:length(R_r)
        fprintf('Calculating R%i ... ', ii);
        tic
        [R_r_i, J_r_i] = R_r_fun{ii}(roVibDensityMatrix, transDipoleMatrixZ,...
            t1TimeProp, t3TimeProp, t2TimeProp, lineshape1D, lineshape2D);
        R_r_time(ii) = toc;
        fprintf('Done.\nTime: %0.2f s\n\n', R_r_time(ii));
        R_r(ii) = {R_r_i};
        J_r(ii) = {J_r_i};    
    end
    R_rtime = sum(R_r_time);
    fprintf('Rephasing Response Caluclations Complete.\nTime: %0.2f s\n\n', R_rtime);
    clear ii R_rtime
    
    fprintf('Calculating Non-Rephasing Response ... \n');
    tic
    for ii = 1:length(R_nr)
        fprintf('Calculating R%i ... ', ii+3);
        tic
        [R_nr_i, J_nr_i] = R_nr_fun{ii}(roVibDensityMatrix, transDipoleMatrixZ,...
            t1TimeProp, t3TimeProp, t2TimeProp, lineshape1D, lineshape2D);
        R_nr_time(ii) = toc;
        fprintf('Done.\nTime: %0.2f s\n\n', R_nr_time(ii));
        R_nr(ii) = {R_nr_i};
        J_nr(ii) = {J_nr_i};
    end
    R_nrtime = sum(R_nr_time);
    fprintf('Non-Rephasing Response Caluclations Complete.\nTime: %0.2f s\n\n', R_nrtime);
    clear ii R_nrtime
end
%%
R_r_tot = R_r{3}-R_r{1}-R_r{2};
R_nr_tot = R_nr{3}-R_nr{1}-R_nr{2};

fig1 = figure(1);
clf
contour(t1s, t3s, real(R_r_tot));
xlabel('t_1 (fs)')
ylabel('t_3 (fs)')
title('R_r')

fig2 = figure(2);
clf
contour(t1s, t3s, real(R_nr_tot));
xlabel('t_1 (fs)')
ylabel('t_3 (fs)')
title('R_{nr}')
%%

n_zp1 = 4*length(t1s);
n_zp3 = 4*length(t3s);

R_r_fft = fftshift(fliplr(circshift(ifft2(R_r_tot, n_zp3, n_zp1),[0 -1])));
R_nr_fft = fftshift(ifft2(R_nr_tot, n_zp3, n_zp1));

R_tot = R_r_fft + R_nr_fft;

freq1 = fftFreqAxis(t1s, 'time_units', 'fs', 'zeropad', n_zp1)+2300;
freq3 = fftFreqAxis(t3s, 'time_units', 'fs', 'zeropad', n_zp3)+2300;

range1 = [2340 2360];
range3 = [2315 2360];

Rcrop = R_tot(freq3 >= range3(1) & freq3 <= range3(2), freq1 >= range1(1) & freq1 <= range1(2));

freq1 = freq1(freq1 >= range1(1) & freq1 <= range1(2));
freq3 = freq3(freq3 >= range3(1) & freq3 <= range3(2));

fig3 = figure(3);
clf
my2dPlot(freq1, freq3, real(Rcrop),'pumpprobe',false,'n_contours',6, 'zlimit', 0.75)

