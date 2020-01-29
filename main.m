%% Intro
% This is the executable for calculating the 2DIR spectum of a linear free
% rotor.
%%
% Defining molecular variables
maxVibLvl = 2;
maxRotLvl = 44;
%% Load the dipole matrix and Hamiltonian for CO2
load('dipoleMatrices.mat');
load('HamiltonianMatrix.mat');

[transDipoleMatrixZ, transDipoleMatrixX, transDipoleMatrixY, Cz, Az, Cx, Ax, Cy, Ay] ...
    = getCO2RoVibTransDipoleMatrix(maxVibLvl, maxRotLvl);

%%
% Setting up the transition dipole matrices and hamitonian

pmodes.H = sparse(roVibCO2HamiltonianMatrix);
pmodes.H_ = pmodes.H;
pmodes.MUX = sparse(transDipoleMatrixX);
pmodes.MUY = sparse(transDipoleMatrixY);
pmodes.MUZ = sparse(transDipoleMatrixZ);
pmodes.A = sparse(Az); %take upper triangular part (???)
pmodes.C = sparse(Cz); %take lower triangular part (???)

% some general things
pmodes.NSTATES = length(pmodes.H);
pmodes.IDENTITY = eye(pmodes.NSTATES);



%% setup response calcs

% experimental parameters
BW = 300; %bandwidth
w_laser = 2350; %central frequency

% setup response parameters
roptions.order = 3;
roptions.w0 = 2350;
roptions.n_t = 2048;%1024;
%roptions.n_t = 128;
roptions.n_zp = 2*roptions.n_t; % the zeropadded length
%roptions.n_zp = 4*roptions.n_t; % the zeropadded length
roptions.dt = 0.025*6;
roptions.t2 = 0; %the population time
roptions.flag_plot = false;
roptions.polarizations = {[0; 0; 1],[0; 0; 1],[0;0;1],[0;0;1] };
roptions.w_laser = w_laser;
roptions.BW = BW;
roptions.c2form = 'voigt'; %'1fast'; 
roptions.c2params = struct('T2',500, 'Delta1_cm', 0.035);
roptions.T = 300; 
roptions.thermal_cutoff = 0.0033;%1e-4;%1e-6;%???
roptions.verbose = 1;

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

%[V,E]=analyzeEnergyLevels(lmodes,pmodes,'roptions',roptions);

%
% calculate response
tic
out = freeRotorResponse(pmodes,roptions);
toc


%% plot
%
%range = [(w_laser-BW/2) (w_laser+BW/2)];
% range1 = [2250 2450];
% range3 = [2100 2450];
range1 = [2290 2390];
range3 = [2290 2390];
ind1 = (out.w1 >= range1(1) & out.w1 <= range1(2));
ind3 = (out.w3 >= range3(1) & out.w3 <= range3(2));
w1 = out.w1(ind1);
w3 = out.w3(ind3);
R = out.R(ind3,ind1);
J = out.J(ind1);
fig = figure(103);
[h1 h2 h3]= my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20,'zlimit',0.3);
% [h1, h2, h3] = my2dPlot(w1,w3,R,'pumpprobe',false,'n_contours',20);
plot(h2,w1,J)
set(h2,'XTickLabel',[],'XLim',[w1(1) w1(end)])
%%
fig.Color = 'w';
% fig.Children.TickDir = 'out'; % sets your ticks to point outwards
% fig.Children.TickLength = 1.5.*fig.Children.TickLength; % bigger ticks (useful as you shrink it)
[fig,~] = resizefig(fig,3.37,'inches'); % makes it 3.37" wide
print(fig,'C:\Users\Kai\OneDrive\Documents\Graduate Research\PreComp\FigureScripts\Figures\2DIR\SimFig.eps','-depsc');