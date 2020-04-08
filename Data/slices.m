dataPath = 'C:\Users\kaigr\OneDrive\Documents\Graduate Research\MATLAB\FreeRotorGPUCPUHybrid\Data';
temp = load2DIRdata(dataPath, [58 63 70 73 76]);
temp = sort2DIRdata(temp, 'sortby', 't2');
expData = calibrate2DIRdata(temp,[0.875,297]);

r1 = [2300 2380];
r3 = [2300 2380];

cpExpData = cropData(expData,r1,r3);
clear temp dataPath
%%

simData = load('C:\Users\kaigr\OneDrive\Documents\Graduate Research\MATLAB\FreeRotorGPUCPUHybrid\Outputs\2020-30-03\25J_75K_6000tJ_out_1.mat');
simData = simData.data;

cpSimData = cropData(simData, r1, r3);
%%
figure(2);
my2dRoVibPlot(cpExpData(1).w1, cpExpData(1).w3, cpExpData(1).R, 'pumpprobe',false,'n_contours',6, 'zlimit', 0.75)

figure(3)
my2dRoVibPlot(cpSimData(1).w1, cpSimData(1).w3, cpSimData(1).R, 'pumpprobe',false,'n_contours',6, 'zlimit', 0.25)
%%
figure(1);
clf
f2 = prepareGlobalFitData(cpExpData);
dm = [1 length(cpExpData)];
freeRotorMatrixPlot2DIR(f2,cpExpData(1).w1, cpExpData(1).w3, [cpExpData(:).t2]./1000, dm, 1);


%%
w1_ind = 113;

legend_list = cell(1, length(cpExpData));

fig1 = figure(100);
clf
hold on

w3NormRange = [2366 2376];

for ii = 1:length(cpExpData)
    plot(cpExpData(ii).w3, cpExpData(ii).R(:, w1_ind)'./abs(min(cpExpData(ii).R(cpExpData(ii).w3 >= w3NormRange(1) & cpExpData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpExpData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpExpData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);


%%
w1_ind = 399;

legend_list = cell(1, length(cpSimData));

fig1 = figure(101);
clf
hold on

w3NormRange = [2366 2376];

for ii = 1:length(cpSimData)
    plot(cpSimData(ii).w3, cpSimData(ii).R(:, w1_ind)'./abs(min(cpSimData(ii).R(cpSimData(ii).w3 >= w3NormRange(1) & cpSimData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpSimData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpSimData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);

%%

w1_ind = 82;

legend_list = cell(1, length(cpExpData));

fig2 = figure(200);
clf
hold on

w3NormRange = [2348 2361];

for ii = 1:length(cpExpData)
    plot(cpExpData(ii).w3, cpExpData(ii).R(:, w1_ind)./abs(min(cpExpData(ii).R(cpExpData(ii).w3 >= w3NormRange(1) & cpExpData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpExpData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpExpData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);

%%

w1_ind = 288;

legend_list = cell(1, length(cpExpData));

fig2 = figure(201);
clf
hold on

w3NormRange = [2348 2361];

for ii = 1:length(cpSimData)
    plot(cpSimData(ii).w3, cpSimData(ii).R(:, w1_ind)'./abs(min(cpSimData(ii).R(cpSimData(ii).w3 >= w3NormRange(1) & cpSimData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpSimData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpSimData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);
%%
w1_ind = 94;

legend_list = cell(1, length(cpExpData));

fig3 = figure(300);
clf
hold on

w3NormRange = [2353 2363];

for ii = 1:length(cpExpData)
    plot(cpExpData(ii).w3, cpExpData(ii).R(:, w1_ind)./abs(min(cpExpData(ii).R(cpExpData(ii).w3 >= w3NormRange(1) & cpExpData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpExpData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpExpData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);

%%

w1_ind = 331;

legend_list = cell(1, length(cpExpData));

fig2 = figure(301);
clf
hold on

w3NormRange = [2353 2363];

for ii = 1:length(cpSimData)
    plot(cpSimData(ii).w3, cpSimData(ii).R(:, w1_ind)'./abs(min(cpSimData(ii).R(cpSimData(ii).w3 >= w3NormRange(1) & cpSimData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpSimData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpSimData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);

%%
w1_ind = 108;

legend_list = cell(1, length(cpExpData));

fig4 = figure(400);
clf
hold on

w3NormRange = [2363 2373];

for ii = 1:length(cpExpData)
    plot(cpExpData(ii).w3, cpExpData(ii).R(:, w1_ind)./abs(min(cpExpData(ii).R(cpExpData(ii).w3 >= w3NormRange(1) & cpExpData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpExpData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpExpData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);

%%

w1_ind = 380;

legend_list = cell(1, length(cpExpData));

fig2 = figure(401);
clf
hold on

w3NormRange = [2363 2373];

for ii = 1:length(cpSimData)
    plot(cpSimData(ii).w3, cpSimData(ii).R(:, w1_ind)'./abs(min(cpSimData(ii).R(cpSimData(ii).w3 >= w3NormRange(1) & cpSimData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpSimData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpSimData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);

%%
w1_ind = 64;

legend_list = cell(1, length(cpExpData));

fig5 = figure(500);
clf
hold on

w3NormRange = [2348 2363];

for ii = 1:length(cpExpData)
    plot(cpExpData(ii).w3, cpExpData(ii).R(:, w1_ind)./abs(min(cpExpData(ii).R(cpExpData(ii).w3 >= w3NormRange(1) & cpExpData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpExpData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpExpData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);

%%

w1_ind = 224;

legend_list = cell(1, length(cpExpData));

fig2 = figure(501);
clf
hold on

w3NormRange = [2348 2363];

for ii = 1:length(cpSimData)
    plot(cpSimData(ii).w3, cpSimData(ii).R(:, w1_ind)'./abs(min(cpSimData(ii).R(cpSimData(ii).w3 >= w3NormRange(1) & cpSimData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpSimData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpSimData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);
