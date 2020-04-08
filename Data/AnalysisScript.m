dataPath = 'C:\Users\kaigr\OneDrive\Documents\Graduate Research\MATLAB\FreeRotorGPUCPUHybrid\Data\';
load(strcat(dataPath, 'nonRephasingResults.mat'));
load(strcat(dataPath, 'rephasingResults.mat'));
load(strcat(dataPath, 'simWorkspace.mat'));
%%
fig1 = figure(1);
clf
plot(t1s, real(Jnr), 'LineWidth', 1.5);
box off
xlabel('t1 (fs)')
ylabel('J_{nr}')

fig2 = figure(2);
clf
plot(t1s, real(Jr), 'LineWidth', 1.5);
box off
xlabel('t1 (fs)')
ylabel('J_{r}')
%%
n_zp1 = 2*length(t1s);
freq1 = fftFreqAxis(t1s, 'time_units', 'fs', 'zeropad', n_zp1)+2300;
Jnr_fft = fftshift(sgrsifft(Jnr, n_zp1));
Jr_fft = fftshift(fliplr(circshift(sgrsifft(Jr, n_zp1),[0 -1])));
Jtot = Jnr_fft + Jr_fft;
fig3 = figure(3);
clf
plot(freq1, real(Jtot), 'LineWidth', 1);
box off
xlim([2300 2400]);
xlabel('\omega_1 / 2\pic')
ylabel('Abs')
%%

fig4 = figure(4);
clf
contour(t1s, t3s, real(Rr));
xlabel('t_1 (fs)')
ylabel('t_3 (fs)')
title('R_r')

fig5 = figure(5);
clf
contour(t1s, t3s, real(Rnr));
xlabel('t_1 (fs)')
ylabel('t_3 (fs)')
title('R_{nr}')
%%
n_zp1 = 2*length(t1s);
n_zp3 = 1*length(t3s);

R_r_fft = fftshift(fliplr(circshift(ifft2(Rr, n_zp3, n_zp1),[0 -1])));
R_nr_fft = fftshift(ifft2(Rnr, n_zp3, n_zp1));

R_tot = R_r_fft + R_nr_fft;

freq1 = fftFreqAxis(t1s, 'time_units', 'fs', 'zeropad', n_zp1)+2300;
freq3 = fftFreqAxis(t3s, 'time_units', 'fs', 'zeropad', n_zp3)+2300;

range1 = [2300 2380];
range3 = [2300 2380];

Rcrop = R_tot(freq3 >= range3(1) & freq3 <= range3(2), freq1 >= range1(1) & freq1 <= range1(2));

freq1 = freq1(freq1 >= range1(1) & freq1 <= range1(2));
freq3 = freq3(freq3 >= range3(1) & freq3 <= range3(2));

fig6 = figure(6);
clf
my2dRoVibPlot(freq1, freq3, real(Rcrop),'pumpprobe',false,'n_contours',8, 'zlimit', .25)

