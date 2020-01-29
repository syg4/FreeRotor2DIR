maxVibLvl = 2;
maxRotLvl = 4;

[v1, v2] = getVibrationalStateMatrices(maxVibLvl, maxRotLvl);
[j1, j2] = getJStateMatrices(maxVibLvl, maxRotLvl);
[m1, m2] = getMStateMatrices(maxVibLvl, maxRotLvl);
%%

[transDipoleMatrixZ, transDipoleMatrixX, transDipoleMatrixY, Cz, Az, Cx,...
    Ax, Cy, Ay] = getRoVibTransDipoleMatrix_Vectorized(v1, v2, j1, j2, m1, m2);

%%

baseline_range = [3985 3990];
range1 = [2295 2400];
peakThreshold = 0.2;
symmetry = 'even';

dataPath = ['C:\Users\kaigr\OneDrive\Documents\Graduate Research\MATLAB'...
    '\FreeRotor2DIR\New\Data\CO2_THF_950um_FlowCell.CSV'];

[PFit, RFit, gasConst] = getEmpericalGasPhaseHamiltonian(dataPath, range1,...
    baseline_range, peakThreshold, symmetry);

roVibHamiltonianMatrix = getRoVibHamiltonianMatrix_Vectorized(v1, v2, j1, j2, m1, m2, gasConst);

fprintf('Do the peak identifications and polynomial fits look good?\n\tIf so press any key to continue ... ')
pause;
fprintf('Continuing with computations!\n\n');

close all
%%
temperature = 300; %K

roVibDensityMatrix = getRoVibDensityMatrix(v1, v2, j1, j2, m1, m2, ...
    roVibHamiltonianMatrix, temperature, 'even')







