%%

data = struct('temp', [], 'name', [], 'frequency', [], 'absorbance', [], 'abs_baselined', []);

if isunix
    path = '/Users/kaigronborg/OneDrive/Documents/Graduate Research/Projects/THFFreeRotor/Simulations/Data/CO2_THF_950um_FlowCell.CSV';
else
    path = 'C:\Users\Kai\OneDrive\Documents\Graduate Research\Projects\THFFreeRotor\Simulations\Data\CO2_THF_950um_FlowCell.CSV';
end
    
data(1).temp = load(path); % loads the data
data(1).name = 'CO_2 in THF'; % assigns the filename
data(1).frequency = data(1).temp(:,1); % unpacks the x variable
data(1).absorbance = data(1).temp(:,2); % unpacks the y variable

%%
range = [3985 3990];

z = mean(data(1).absorbance(data(1).frequency>range(1) & data(1).frequency<range(2))); % finds the average value of absorbance over that range
data(1).abs_baselined = data(1).absorbance - z; % defines a new variable that has subtracted the mean value found above
clear index z

data(1).name = strrep(data(1).name, '_', ' ');
data(1).name = strrep(data(1).name, '.CSV', '');

%%

range2 = [2295 2400];

dataCropped(1).name = data(1).name;
dataCropped(1).frequency = data(1).frequency(data(1).frequency >= range2(1) & data(1).frequency <= range2(2));
dataCropped(1).absorbance = data(1).abs_baselined(data(1).frequency >= range2(1) & data(1).frequency <= range2(2));

%%
load('dipoleMatrices.mat');
load('HamiltonianMatrix3.mat');
%%
%dipoleMat = sparse(transDipoleMatrix);
%HamMat = sparse(roVibCO2HamiltonianMatrix);

% I am confused about what A and C should be in this case. Probably
% something here is not right. For now, I will let MUY and MUZ = 0, set MUX
% to transDipoleMatrix and H and H_ to roVibCO2HamiltonianMatrix and see
% what I get.

pmodes.H = sparse(roVibCO2HamiltonianMatrix);
pmodes.H_ = pmodes.H;
pmodes.MUX = sparse(transDipoleMatrixX);
pmodes.MUY = sparse(transDipoleMatrixY); %sparse(length(pmodes.H),length(pmodes.H));
pmodes.MUZ = sparse(transDipoleMatrixZ); %sparse(length(pmodes.H),length(pmodes.H));
pmodes.A = triu(pmodes.MUZ,1); %take upper triangular part (???)
pmodes.C = tril(pmodes.MUZ,1); %take lower triangular part (???)

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
roptions.n_t = 2048;%1024;%
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
roptions.c2params = struct('T2',5000, 'Delta1_cm', 0.03);
roptions.T = 300;
roptions.thermal_cutoff = 1e-4;%0.001;%1e-6;%???
roptions.verbose = 0;

% set up lineshape function
lineshape = chooseLineshapeFunction(roptions.c2form,roptions.c2params);
roptions.g = lineshape.g;

%[V,E]=analyzeEnergyLevels(lmodes,pmodes,'roptions',roptions);



%% Generagting First Order Response

wavenumbersToInvPs = 0.0300;
k_B_SI = 1.3806E-23;
h = 6.6260e-34;
c_SI = 299790000;

k_B_cm_K = k_B_SI/h/c_SI/100;%k_B in cm-1/K
thermal_cutoff = 0.01;
T = 0;
verbose = 0;

out = [];
n_sparse_states = estimateNSparseStates(pmodes,roptions);
n_sparse_states = min(n_sparse_states,pmodes.NSTATES-2);
order = roptions.order;

flag_print = false;

%n_excitons = 2;
%n_exciton_sig_figs = 1;
n_exciton_tol = 0.0001;

if isfield(roptions,'w0')
    if isempty(roptions.w0)
        w0=0;
    else
        w0 = roptions.w0;
    end
end
if isfield(roptions,'n_sparse_states')
    if isempty(roptions.n_sparse_states)
        %do nothing
    else
        n_sparse_states = roptions.n_sparse_states;
    end
end
if isfield(roptions,'flag_plot')
    flag_plot = roptions.flag_plot;
end
if flag_print
    fid = 1;
else
    if isunix
        fid = fopen('/dev/null');
%         fid = fopen('out.txt', 'w');
    else
        warning('this windows "dev/null" is not tested yet');
        fid = fopen('out.txt', 'w'); %fopen('nul');
    end
end
if isfield(roptions,'T')
    if isempty(roptions.T)
        %do nothing (see default at top)
    else
        T = roptions.T;
    end
end
if isfield(roptions,'thermal_cutoff')
    if isempty(roptions.thermal_cutoff)
        %do nothing (see default at top)
    else
        thermal_cutoff = roptions.thermal_cutoff;
    end
end
if isfield(roptions,'verbose')
    if isempty(roptions.verbose) || roptions.verbose == 0 
        verbose=0;
    else
        verbose = 1;
    end
end

% simulation parameters
if verbose>=1,disp('initialize variables'),end
n_t = roptions.n_t;
n_zp = roptions.n_zp;
dt = roptions.dt; %time step
t2 = roptions.t2; %population time (ps)
pol = roptions.polarizations;
e_1 = pol{1};
e_2 = pol{2};
e_3 = pol{3};
e_4 = pol{4};
w_laser = roptions.w_laser;
BW = roptions.BW;

% set up time and response functions
J = zeros(1,n_t);
J_accum = zeros(1,n_t);
R_r = zeros(n_t,n_t);
R_nr = zeros(n_t,n_t);
R_r_accum = zeros(n_t,n_t);
R_nr_accum = zeros(n_t,n_t);
t=0:dt:(n_t-1)*dt;
[T1,T3] = meshgrid(t,t);

% set up thermal weight
kT = k_B_cm_K*T;


%upack the results
if verbose>=1,disp('unpack pmodes'),end
f=fieldnames(pmodes);
for ii=1:length(f)
    eval(strcat(f{ii},'=pmodes.',f{ii},';'))
end

flag_diagonal_input = false;
if issparse(H_)
    % calculate the eigenvectors of the coupled system
    if verbose>=1,disp('calculate sparse eigenvalues'),end
    
    if isdiag(H_)
        disp('-------- RUNNING --------');
        disp('H_ already diagonal');
        flag_diagonal_input = true;
        
        E = full(diag(H_));
        V = speye(length(E),length(E));
        VV = speye(length(E),length(E)); %eigenvectors in eigenstate basis
        
    else
        
        [V,E]=eigs(H_,n_sparse_states,'SM');
        E = diag(E);
        if verbose>=1,disp('sort energies and vectors'),end
        [E,ordering] = sort(E);
        E = E - E(1); %remove zero point energy
        V = V(:,ordering); %eigenvectors in input basis
        VV = speye(length(E),length(E)); %eigenvectors in eigenstate basis
        
    end
else
    % calculate the eigenvectors of the coupled system
    if verbose>=1,disp('calculate full eigenvalues'),end
    [V,E]=eig(H_,'vector');
    if verbose>=1,disp('sort energies and vectors'),end
    [E,ordering] = sort(E);
    E = E - E(1); %remove zero point energy
    V = V(:,ordering); %eigenvectors in input basis
    VV = eye(size(V)); %eigenvectors in eigenstate basis
end

%
% set up operators
%
if verbose>=1,disp('set up operators'),end

if flag_diagonal_input
    
    A0 = A;
    C0 = C;
    
else
    
    % rotate creation annihil operators to the eigenstate basis
    A0 = V'*A*V;
    C0 = V'*C*V;
    
    % rotate dipole operators to the eigenstate basis
    MUX = V'*MUX*V;
    MUY = V'*MUY*V;
    MUZ = V'*MUZ*V;
    
end

% could add a loop over possible initial states here. The idea would be to
% look at the thermal density matrix elements relative to some cutoff. then
% loop through the response function calculation for each state with the
% appropriate thermal weight.

flag_finished_thermal_loop = false;

i_thermal = 0;
J_rot = 0;
m_rot = 0;
while ~flag_finished_thermal_loop
    if mod(J_rot, 2) ~= 0
        i_thermal = i_thermal + 1;
        if m_rot == J_rot
            J_rot = J_rot+1;
            m_rot = -J_rot;
        else
            m_rot = m_rot + 1;
        end
    else
        i_thermal = i_thermal+1;
        
        if T == 0
            flag_finished_thermal_loop = true;
            thermal_weight = 1;
        else
            thermal_weight = full(exp(-E(i_thermal)/kT));
            if thermal_weight < thermal_cutoff
                flag_finished_thermal_loop = true;
                continue;
            end
        end
        if verbose>=1,fprintf(fid,'loop over thermal states, i_thermal = %i; J = %i; m = %i\n',i_thermal, J_rot, m_rot);end
        if verbose>=1,fprintf(fid,'    thermal weight = %f\n',thermal_weight);end
        
        if verbose>=1,fprintf(1,'loop over thermal states, i_thermal = %i; J = %i; m = %i\n',i_thermal, J_rot, m_rot);end
        if verbose>=1,fprintf(1,'    thermal weight = %f\n',thermal_weight);end
        
        if m_rot == J_rot
            J_rot = J_rot + 1;
            m_rot = -J_rot;
        else
            m_rot = m_rot + 1;
        end
        
        if verbose>=1,disp('determine active states...'),end
        
        % density matrix
        PSIi = VV(:,i_thermal); %take first eigenstate for the time being
        %rho = PSIi*PSIi'; % could do thermal density here!
        
        %one way to go would be to define mui muj etc from inputs and then add
        %invariants function (see thoughts below)
        %tests:
        %reproduce Rhcomplex
        %make sure in strong mixing dipoles are orthogonal
        % check amplitudes of parallel and perp polarizations
        
        % calc mus and omegas from inputs (ultimately want to refactor this)
        % calculate the one and two exciton manifolds. Might need to be modified if
        % thermal states are allowed. not sure.
        
        % %find all one and two exciton states
        % [ind_1ex ind_2ex] =  findNExcitonStates(PSIi,C0,n_exciton_sig_figs);
        %
        % %keep only the ones in the laser bandwidth
        % [ind_1ex ind_2ex] = filterExcitons(w_laser,BW,E,i_thermal,ind_1ex,ind_2ex);
        
        
        lower_limit_energy = roptions.w_laser-roptions.BW/2;
        upper_limit_energy = roptions.w_laser+roptions.BW/2;
        
        % calculate dipole matrix elements
        
        % can probably speed this up by looking only for nonzero elements
        if flag_diagonal_input
            
            ind_1ex = find(sum([MUX(i_thermal,:).^2; MUY(i_thermal,:).^2; MUZ(i_thermal,:).^2],1)>n_exciton_tol);
            ind_2ex = find(sum([MUX(ind_1ex,:).^2; MUY(ind_1ex,:).^2; MUZ(ind_1ex,:).^2],1)>n_exciton_tol);
            
            n = length(ind_1ex);
            n2 = length(ind_2ex);
            
            mu1 = full([MUX(i_thermal,ind_1ex)' MUY(i_thermal,ind_1ex)' MUZ(i_thermal,ind_1ex)']);
            mu2 = zeros(n2,3,n);
            for ii = 1:n
                for jj = 1:n2
                    mu2(jj,:,ii) = [MUX(ind_1ex(ii),ind_2ex(jj)) MUY(ind_1ex(ii),ind_2ex(jj)) MUZ(ind_1ex(ii),ind_2ex(jj))];
                end
            end
            
            mu0 = mu2(1:length(mu2)/2,:,:);
            mu2 = mu2(length(mu2)/2+1:end,:,:);
            
            ind_1ex = ind_1ex(:);
            jj = 1;
            kk = 1;
            ind_2ex_loop = ind_2ex(:);
            ind_2ex = [];
            ind_0ex = [];
            for ii = 1:length(ind_2ex_loop)
                if ind_2ex_loop(ii) < 2*(length(H_)/3)
                    ind_0ex(jj) = ind_2ex_loop(ii);
                    jj = jj+1;
                else
                    ind_2ex(kk) = ind_2ex_loop(ii);
                    kk = kk+1;
                end
            end
            
            n0 = length(ind_0ex);
            n1 = length(ind_1ex);
            n2 = length(ind_2ex);
            
            
            % I don't know if I need to put the energy filtering back in here...
            % SGR
        else
            
            fprintf('This One')
            n = length(E);
            mu1 = zeros(n,3);
            mu2 = zeros(n,3,n);
            
            for ii = i_thermal+1:n-1
                PSIf = VV(:,ii);
                Ef = E(ii)-E(i_thermal);
                
                %filter out energies outside our window
                if (Ef < lower_limit_energy)||(Ef > upper_limit_energy)
                    continue;
                end
                
                mu1(ii,:) = [PSIf'*MUX*PSIi PSIf'*MUY*PSIi PSIf'*MUZ*PSIi];
                for jj =  (ii+1):n
                    PSIf2 = VV(:,jj);
                    
                    Ef2 = E(jj)-Ef-E(i_thermal);
                    
                    %filter out energies outside our window
                    if (Ef2 < lower_limit_energy)||(Ef2 > upper_limit_energy)
                        continue;
                    end
                    
                    mu2(jj,:,ii) = [PSIf2'*MUX*PSIf PSIf2'*MUY*PSIf PSIf2'*MUZ*PSIf];
                end
            end
            
            ind_1ex = find(sum(mu1.^2,2)>n_exciton_tol);
            ind_2ex = [];
            for ii = 1:length(ind_1ex)
                ind_2ex = [ind_2ex; find(sum(mu2(:,:,ind_1ex(ii)).^2,2)>n_exciton_tol)];
            end
            ind_2ex = unique(ind_2ex);
            
            ind_2ex_loop = ind_2ex(:);
            ind_2ex = [];
            ind_0ex = [];
            jj = 0;
            kk = 0;
            
            for ii = 1:length(ind_2ex_loop)
                if ind_2ex_loop(ii) < 2*(length(H_)/3)
                    jj = jj + 1;
                    ind_0ex(jj) = ind_2ex_loop(ii);
                else
                    kk = kk + 1;
                    ind_2ex(kk) = ind_2ex_loop(ii);
                end
            end
            
            
            %reduce dipole matrix elements to only the needed size
            mu1 = mu1(ind_1ex,:);
            mu2 = mu2(ind_2ex,:,ind_1ex);
            
            mu0 = mu2(1:length(mu2)/2,:,:);
            mu2 = mu2(length(mu2)/2+1:end, :, :);
            
            n0 = length(ind_0ex);
            n1 = length(ind_1ex);
            n2 = length(ind_2ex);
            
        end
        
        % ind_1ex = find(abs(round(C0*PSIi,1))>0);
        % ind_2ex = find(abs(round(C0*C0*PSIi,1))>0);
        
        % energies -- subtract zero point energy
        w0s = E(ind_0ex) - E(i_thermal);
        w1s = E(ind_1ex) - E(i_thermal);
        w2s = E(ind_2ex) - E(i_thermal);
        
        fprintf(fid,'zero exciton state energies\n');
        fprintf(fid,'%d\t%8.1f\n',[ind_0ex',w0s]');
        fprintf(fid,'one exciton state energies\n');
        fprintf(fid,'%d\t%8.1f\n',[ind_1ex,w1s]');
        fprintf(fid,'two exciton state energies\n');
        fprintf(fid,'%d\t%8.1f\n',[ind_2ex',w2s]');
        fprintf(fid,'\n');
        
        % subtract rotating frame frequency and convert to rad/ps
        w0s = w0s * 2 * pi * wavenumbersToInvPs;
        w1s = (w1s - w0)*2*pi*wavenumbersToInvPs;
        w2s = (w2s - 2*w0)*2*pi*wavenumbersToInvPs;
        
        if verbose>=1,disp('determine matrix elements...'),end
        
        % % calculate dipole matrix elements
        % mu = zeros(n,3);
        % mu2 = zeros(n2,3,n);
        % for ii = 1:n
        %     PSIf = VV(:,ind_1ex(ii));
        %     mu(ii,:) = [PSIf'*MUX*PSIi PSIf'*MUY*PSIi PSIf'*MUZ*PSIi];
        %     for jj =  1:n2
        %         PSIf2 = VV(:,ind_2ex(jj));
        %         mu2(jj,:,ii) = [PSIf2'*MUX*PSIf PSIf2'*MUY*PSIf PSIf2'*MUZ*PSIf];
        %     end
        % end
        
        g = @(t) roptions.g(t,roptions.c2params);
        
        if verbose>=1,disp('start linear spectroscopy...'),end
        %linear spectroscopy
        for j = 1:n
            [~,mu0j] = unit_vector(mu1(j,:));
            J = J + mu0j^2.*exp(-1i*w1s(j).*t);
        end
        %add lineshape (same for all peaks for now)
        J = J.*exp(-g(t));
    end
    
    J_accum = J_accum + J*thermal_weight;
end

J = J_accum;
figure(2)
plot(t, real(J));

% calculate 1D spectrum (freq domain)
J = real(fftshift(sgrsifft(J,n_zp)));

freq  = fftFreqAxis(t,'time_units','ps','zeropad',n_zp);
freq = freq+w0;

disp('-------- FINISHED --------')
%%
figure(1)
plot(freq, J./max(J),'.-', dataCropped(1).frequency, dataCropped(1).absorbance./max(dataCropped(1).absorbance), '.-');
%%
mu = 2349;
sig = 3;
gamma = 25;
gaussian = @(x) 1./(sig.*sqrt(2*pi)) .* exp(-1/2 .* ((x - mu) ./ sig).^2) .* (1/pi) .* (1/2 .* gamma)./((x - mu).^2 + (1/2 .* gamma)^2);


fig = figure(3);
hold on
plot(dataCropped(1).frequency(dataCropped(1).frequency >= 2349.5), dataCropped(1).absorbance(dataCropped(1).frequency >= 2349.5)./max(dataCropped(1).absorbance), 'r', 'LineWidth',1);
plot(dataCropped(1).frequency(dataCropped(1).frequency <= 2349.6), dataCropped(1).absorbance(dataCropped(1).frequency <= 2349.6)./max(dataCropped(1).absorbance), 'color', rgb('dark purple'), 'LineWidth',1);
plot(dataCropped(1).frequency, gaussian(dataCropped(1).frequency)./max(gaussian(dataCropped(1).frequency)), 'b', 'LineWidth',2)
hold off
xlim([dataCropped(1).frequency(1) 2385])
ylim([-0.05 1.05])
xlabel('\omega_1/2\pic')
ylabel('A (a.u.)')
fig.Color = 'w';
% set(findall(gca, 'Type', 'Line'),'LineWidth',1);
%%
function [out,n] = unit_vector(in)
n  = norm(in,2);
if n>1e-6
    out = in./n;
else
    n=0;
    out = 0*in;
end
out = out(:); %convert to a column matrix
end

function n_sparse_states = estimateNSparseStates(pmodes,options)
    original_energies = diag(pmodes.H);
    original_energies = original_energies-original_energies(1);
    max_e = 2*(options.w_laser+options.BW);
    n_sparse_states = sum(original_energies<=max_e);
end
    