%% This properly handles the symmetry requirements of the wavefunction
% allowing only the even numbers of rotational states
% to be populated in the ground state allowing for the correct values of
% the gas constants to be calculated
function [PFit, RFit, gasConst] = getEmpericalGasPhaseHamiltonian(dataPath, freq_range, baseline_range, peakThreshold, symmetry)

data = struct('temp', [], 'name', [], 'frequency', [], 'absorbance', [], 'abs_baselined', []);

files = dir(dataPath);

data(1).temp = load(files(1).name); % loads the data
data(1).name = files(1).name; % assigns the filename
data(1).frequency = data(1).temp(:,1); % unpacks the x variable
data(1).absorbance = data(1).temp(:,2); % unpacks the y variable

%%

z = mean(data(1).absorbance(data(1).frequency>baseline_range(1) & data(1).frequency<baseline_range(2))); % finds the average value of absorbance over that range
data(1).abs_baselined = data(1).absorbance - z; % defines a new variable that has subtracted the mean value found above
clear index z

data(1).name = strrep(data(1).name, '_', ' ');
data(1).name = strrep(data(1).name, '.CSV', '');

%%

dataCropped(1).name = data(1).name;
dataCropped(1).frequency = data(1).frequency(data(1).frequency >= freq_range(1) & data(1).frequency <= freq_range(2));
dataCropped(1).absorbance = data(1).abs_baselined(data(1).frequency >= freq_range(1) & data(1).frequency <= freq_range(2));

%%

[~, locs] = findpeaks(dataCropped(1).absorbance);

peakAbs = dataCropped(1).absorbance(locs);
peakFreq = dataCropped(1).frequency(locs);
peakFreq = peakFreq(peakAbs >= peakThreshold);
peakAbs = peakAbs(peakAbs >= peakThreshold);

fig1 = figure(100);
clf
plot(dataCropped(1).frequency, dataCropped(1).absorbance, peakFreq, peakAbs, 'o');
xlim([freq_range(1) 2385]);
ylim([-0.05 1.4]);
xlabel('\omega_1/2\pic');
ylabel('A');
box off
fig1.Children.TickDir = 'out'; % sets your ticks to point outwards
fig1.Children.TickLength = 1.5.*fig1.Children.TickLength;
fig1.Color = 'w';




%%
switch symmetry
    case 'even'
        PFreqs = zeros(idivide(length(peakAbs), int32(2), 'floor'), 1);
        PJNums = zeros(idivide(length(peakAbs), int32(2), 'floor'), 1);
        PJNum = idivide(length(peakAbs), int32(2), 'floor') * 2;
        iP = 1;
        
        RFreqs = zeros(idivide(length(peakAbs), int32(2), 'ceil'), 1);
        RJNums = zeros(idivide(length(peakAbs), int32(2), 'ceil'), 1);
        RJNum = 0;
        iR = 1;
        
        for jj = 1:length(peakAbs)
            
            if peakFreq(jj) >= 2349
                RFreqs(iR) = peakFreq(jj);
                RJNums(iR) = RJNum;
                RJNum = RJNum + 2;
                iR = iR + 1;
            else
                PFreqs(iP) = peakFreq(jj);
                PJNums(iP) = PJNum;
                PJNum = PJNum - 2;
                iP = iP + 1;
            end
        end
    case 'odd'
        % TODO: Implement case for only 'odd' rotational ground states
        % being populated
        
    case 'all'
        % TODO: Implement case for 'all' rotational ground states being
        % populated
        
    otherwise
        error('Symmetry must be either "even", "odd", or "all".')
end

%%
PFit = fit(PJNums, PFreqs, 'poly3');
RFit = fit(RJNums, RFreqs, 'poly3');

%%
fig2 = figure(200);
clf
box off
fig2.Children.TickDir = 'out'; % sets your ticks to point outwards
fig2.Children.TickLength = 1.5.*fig2.Children.TickLength;
fig2.Color = 'w';

hold on
plot(PJNums, PFreqs, 'co', RJNums, RFreqs, 'ro');
plot(PFit, 'c-')
plot(RFit, 'r-')
xlabel('J"')
ylabel('\omega_1/2\pic')
legend({'P Data', 'R Data', 'P Fit', 'R Fit'}, 'Location', 'northwest')
hold off



%% Determing Gas Constants

gasConst = struct('D', 0, 'a_e', 0, 'B_e', 0, 'I', 0);

c = 2.997900000000000e+10;
h = 6.626000000000000e-34;
%     k_B = 1.380000000000000e-23;
% Equations used to fit

% v_P(J") = (v_e - 2*v_e*x_e) - (2*B_e - 2*a_e)*J" - a_e*J"^2 + 4*D*J"^3

% v_R(J") = (v_e - 2*v_e*x_e + 2*B_e - 3*a_e - 4*D) + (2*B_e - 4*a_e -
% 12*D)*J" - (a_e + 12*D)*J"^2 - 4*D*J"^3

PD = PFit.p1 ./ 4;
Pa_e = -1*PFit.p2;
PB_e = Pa_e - PFit.p3/2;

RD = -1*RFit.p1./4;
Ra_e = (-1*RFit.p2)-12*RD;
RB_e = (RFit.p3 + 4*Ra_e + 12*RD)/2;


DAvg = (PD+RD)/2;
a_eAvg = (Pa_e + Ra_e)/2;
B_eAvg = (PB_e + RB_e)/2;

I = h/(800*pi^2*c*B_eAvg);

gasConst.D = DAvg;
gasConst.a_e = a_eAvg - 0.000347855;
gasConst.B_e = B_eAvg;
gasConst.I = I;
gasConst.v_ex_e = 12; %12;
gasConst.v_e = PFit.p4 + 2*gasConst.v_ex_e;