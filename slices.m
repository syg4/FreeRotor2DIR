dataPath = 'C:\Users\kaigr\OneDrive\Documents\Graduate Research\MATLAB\FreeRotor2DIR\Data';
temp = load2DIRdata(dataPath, [58 63 70 73 76]);
temp = sort2DIRdata(temp, 'sortby', 't2');
data = calibrate2DIRdata(temp,[0.875,297]);

r1 = [2300 2380];
r3 = [2300 2380];

cpData = cropData(data,r1,r3);
%%
figure(2);
my2dRoVibPlot(cpData(1).w1, cpData(1).w3, cpData(1).R, 'pumpprobe',false,'n_contours',6, 'zlimit', 0.75)

%%
figure(1);
clf
f2 = prepareGlobalFitData(cpData);
dm = [1 length(cpData)];
freeRotorMatrixPlot2DIR(f2,cpData(1).w1, cpData(1).w3, [cpData(:).t2]./1000, dm, 1);
%%
w1_ind = 113;

legend_list = cell(1, length(cpData));

fig1 = figure(100);
clf
hold on

w3NormRange = [2366 2376];

for ii = 1:length(cpData)
    plot(cpData(ii).w3, cpData(ii).R(:, w1_ind)'./abs(min(cpData(ii).R(cpData(ii).w3 >= w3NormRange(1) & cpData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);

%%

w1_ind = 82;

legend_list = cell(1, length(cpData));

fig2 = figure(200);
clf
hold on

w3NormRange = [2348 2361];

for ii = 1:length(cpData)
    plot(cpData(ii).w3, cpData(ii).R(:, w1_ind)./abs(min(cpData(ii).R(cpData(ii).w3 >= w3NormRange(1) & cpData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);

%%
w1_ind = 94;

legend_list = cell(1, length(cpData));

fig3 = figure(300);
clf
hold on

w3NormRange = [2353 2363];

for ii = 1:length(cpData)
    plot(cpData(ii).w3, cpData(ii).R(:, w1_ind)./abs(min(cpData(ii).R(cpData(ii).w3 >= w3NormRange(1) & cpData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);

%%
w1_ind = 108;

legend_list = cell(1, length(cpData));

fig4 = figure(400);
clf
hold on

w3NormRange = [2363 2373];

for ii = 1:length(cpData)
    plot(cpData(ii).w3, cpData(ii).R(:, w1_ind)./abs(min(cpData(ii).R(cpData(ii).w3 >= w3NormRange(1) & cpData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);

%%
w1_ind = 64;

legend_list = cell(1, length(cpData));

fig5 = figure(500);
clf
hold on

w3NormRange = [2348 2363];

for ii = 1:length(cpData)
    plot(cpData(ii).w3, cpData(ii).R(:, w1_ind)./abs(min(cpData(ii).R(cpData(ii).w3 >= w3NormRange(1) & cpData(ii).w3 <= w3NormRange(2), w1_ind))), 'LineWidth', 1.5);
    legend_list(ii) = {sprintf('%i ps', cpData(ii).t2)};
end

titleStr = strcat(strcat('\omega_1/2\pic = ', sprintf(' %0.1f', cpData(1).w1(w1_ind))), ' cm^{-1}');
title(titleStr)

% title('\omega_1/2\pic = %f cm^{-1}', cpData(1).w1(w1_ind))
xlabel('\omega_3/2\pic');
ylabel('\DeltamOD');
legend(legend_list);