function out = freeRotorResponse(pmodes,roptions)
% global wavenumbersToInvPs k_B_SI h c_SI
%% This fixes the handling of the dipoles

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

%canonical results
% Rh.n = 2;
% Rh.w = [1889.184 1947.621];
% Rh.mu = [sqrt(0.908) 0 0 ; 0 sqrt(0.671) 0];
% Rh.w_off = 1915;
% Rh.w = Rh.w - Rh.w_off;
% Rh.n2 = 3;
% Rh.w2 = [3767.517 3883.911 3813.698];
% Rh.w2 = Rh.w2 - 2*Rh.w_off;
% Rh.mu2 = zeros(Rh.n2,3,Rh.n);
% Rh.mu2(:,:,1) = [sqrt(2).*Rh.mu(1,:) ; 0,0,0 ; Rh.mu(2,:)];
% Rh.mu2(:,:,2) = [ 0,0,0 ; sqrt(2).*Rh.mu(2,:) ; Rh.mu(1,:)];

% default values
flag_plot = false;
flag_print = false;
w0 = 0;
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
%         fid = fopen('/dev/null');
        fid = fopen('out.txt', 'w');
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
    if isempty(roptions.verbose)
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
        
        if verbose>=1,disp('start third order spectroscopy...'),end
        % first calculate all rephasing diagrams
        fprintf(fid,'Rephasing transitions\n');
%         fprintf(fid,'i\tj\tk\tw_1\t(w_2)\tw_3\tu\n');
        fprintf(fid,'Rn\t%-8s%-8s%-8s%-9s%-9s%-9s%-11s\n', 'i', 'j', 'k', 'w_1', 'w_2', 'w_3', 'u');
        
        for j = 1:n1
            for i  = 1:n1
                
%                 [aa,mui] = unit_vector(mu1(i,:));
%                 [bb,muj] = unit_vector(mu1(j,:));
%                 dipole = mui^2*muj^2;
%                 angle = polarizationInvariant(e_1,e_2,e_3,e_4,...
%                     aa,bb,aa,bb);
%                 
% %                             n0_start = 1;
% %                             if length(w0s) > 2
% %                                 n0_start = i;
% %                             end
%                 
%                 ind_1ex_sub = ind_1ex - length(H_)./3;
%                 ind_dif_i = ind_0ex - ind_1ex_sub(i);
%                 ind_0ex_neg_i = [];
%                 ind_0ex_pos_i = [];
%                 ii = 0;
%                 jj = 0;
%                 for m = 1:length(ind_dif_i)
%                     if ind_dif_i(m) < 0
%                         ii = ii + 1;
%                         ind_0ex_neg_i(ii) = ind_0ex(m);
%                     elseif ind_dif_i(m) > 0
%                         jj = jj + 1;
%                         ind_0ex_pos_i(jj) = ind_0ex(m);
%                     end
%                 end
%                 
%                 if isempty(ind_0ex_neg_i)
%                     ind_0ex_access_i = [ind_0ex(find(ind_0ex == min(ind_0ex_pos_i)))];
%                     w0s_access_i = [w0s(find(ind_0ex == min(ind_0ex_pos_i)))];
%                 elseif isempty(ind_0ex_pos_i)
%                     ind_0ex_access_i = [ind_0ex(find(ind_0ex == max(ind_0ex_neg_i)))];
%                     w0s_access_i = [w0s(find(ind_0ex == max(ind_0ex_neg_i)))];
%                 else
%                     ind_0ex_access_i = [ind_0ex(find(ind_0ex == max(ind_0ex_neg_i))) ind_0ex(find(ind_0ex == min(ind_0ex_pos_i)))];
%                     w0s_access_i = [w0s(find(ind_0ex == max(ind_0ex_neg_i))) w0s(find(ind_0ex == min(ind_0ex_pos_i)))];
%                 end
%                 
%                 ind_dif_j = ind_0ex - ind_1ex_sub(j);
%                 ind_0ex_neg_j = [];
%                 ind_0ex_pos_j = [];
%                 ii = 0;
%                 jj = 0;
%                 for m = 1:length(ind_dif_j)
%                     if ind_dif_j(m) < 0
%                         ii = ii + 1;
%                         ind_0ex_neg_j(ii) = ind_0ex(m);
%                     elseif ind_dif_j(m) > 0
%                         jj = jj + 1;
%                         ind_0ex_pos_j(jj) = ind_0ex(m);
%                     end
%                 end
%                 
%                 if isempty(ind_0ex_neg_j)
%                     ind_0ex_access_j = [ind_0ex(find(ind_0ex == min(ind_0ex_pos_j)))];
%                     w0s_access_j = [w0s(find(ind_0ex == min(ind_0ex_pos_j)))];
%                 elseif isempty(ind_0ex_pos_j)
%                     ind_0ex_access_j = [ind_0ex(find(ind_0ex == max(ind_0ex_neg_j)))];
%                     w0s_access_j = [w0s(find(ind_0ex == max(ind_0ex_neg_j)))];
%                 else
%                     ind_0ex_access_j = [ind_0ex(find(ind_0ex == max(ind_0ex_neg_j))) ind_0ex(find(ind_0ex == min(ind_0ex_pos_j)))];
%                     w0s_access_j = [w0s(find(ind_0ex == max(ind_0ex_neg_j))) w0s(find(ind_0ex == min(ind_0ex_pos_j)))];
%                 end
%                 
%                 ind_0ex_shared = intersect(ind_0ex_access_i, ind_0ex_access_j);
%                 w0s_shared = intersect(w0s_access_i, w0s_access_j);
                
                for l = 1:n0
                    
                    [aa,mu0i] = unit_vector(mu1(i,:));
                    [bb,mu0j] = unit_vector(mu1(j,:));
                    [cc,mujl] = unit_vector(mu0(l,:,j));
                    [dd,muil] = unit_vector(mu0(l,:,i));
                    dipole = mu0i*mu0j*mujl*muil;
                    angle = polarizationInvariant(e_1,e_2,e_3,e_4,...
                        aa,bb,cc,dd);
                    
                    % rephasing diagram R1
                    fprintf(fid,'R1\t%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(i),ind_0ex(l),w1s(j)./(2*pi*wavenumbersToInvPs)+w0,(w1s(j)-w1s(i))./(2*pi*wavenumbersToInvPs),(w1s(i)-w0s(l))./(2*pi*wavenumbersToInvPs)+w0,-dipole*angle);
                    R_r = R_r - dipole*angle*exp(+ 1i*w1s(j).*T1 ...
                        + 1i*(w1s(j)-w1s(i))*t2 ...
                        - 1i*(w1s(i)-w0s(l)).*T3);
                    
                    % rephasing diagram R2
                    fprintf(fid,'R2\t%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_0ex(l),ind_1ex(i),w1s(j)./(2*pi*wavenumbersToInvPs)+w0,w0s(l)./(2*pi*wavenumbersToInvPs),(w1s(i)-w0s(l))./(2*pi*wavenumbersToInvPs)+w0,-dipole*angle);
                    R_r = R_r - dipole*angle*exp(+ 1i*w1s(j).*T1 ...
                        + 1i*w0s(l).*t2 ...
                        - 1i*(w1s(i)-w0s(l)).*T3);
                end
                
                %             n2_start = 1;
                %             if length(w1s) > 1;
                %                 n2_start = i;
                %             end
                for k = 1:n2
                    %molecular dipoles?
                    [cc,muik_] = unit_vector(mu2(k,:,i));
                    [dd,mujk_] = unit_vector(mu2(k,:,j));
                    dipole = mu0i*mu0j*muik_*mujk_;
                    angle = polarizationInvariant(e_1,e_2,e_3,e_4,...
                        aa,bb,cc,dd);
                    
                    %rephasing diagram R3
                    fprintf(fid,'R3\t%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(i),ind_2ex(k),w1s(j)./(2*pi*wavenumbersToInvPs)+w0,(w1s(j)-w1s(i))./(2*pi*wavenumbersToInvPs),(w2s(k)-w1s(j))./(2*pi*wavenumbersToInvPs)+w0,dipole*angle);
                    R_r = R_r + dipole*angle*exp(+ 1i*w1s(j).*T1 ...
                        + 1i*(w1s(j)-w1s(i)).*t2 ...
                        - 1i*(w2s(k)-w1s(j)).*T3);
                end
            end
        end
        % add lineshape (same for all peaks for now)
        R_r = exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*R_r;
        fprintf(fid,'\n\n');
        
        % now non-rephasing diagrams
        fprintf(fid,'Non-rephasing transitions\n');
%         fprintf(fid,'i\tj\tk\tw_1\t(w_2)\tw_3\tu\n');
        fprintf(fid,'Rn\t%-8s%-8s%-8s%-9s%-9s%-9s%-11s\n', 'i', 'j', 'k', 'w_1', 'w_2', 'w_3', 'u');

        for j = 1:n1
            for i  = 1:n1
%                 [aa,mui] = unit_vector(mu1(i,:));
%                 [bb,muj] = unit_vector(mu1(j,:));
%                 dipole = mui^2*muj^2;
%                 angle = polarizationInvariant(e_1,e_2,e_3,e_4,...
%                     aa,bb,aa,bb);
%                 
%                 ind_1ex_sub = ind_1ex - length(H_)./3;
%                 ind_dif_i = ind_0ex - ind_1ex_sub(i);
%                 ind_0ex_neg_i = [];
%                 ind_0ex_pos_i = [];
%                 ii = 0;
%                 jj = 0;
%                 for m = 1:length(ind_dif_i)
%                     if ind_dif_i(m) < 0
%                         ii = ii + 1;
%                         ind_0ex_neg_i(ii) = ind_0ex(m);
%                     elseif ind_dif_i(m) > 0
%                         jj = jj + 1;
%                         ind_0ex_pos_i(jj) = ind_0ex(m);
%                     end
%                 end
%                 
%                 if isempty(ind_0ex_neg_i)
%                     ind_0ex_access_i = [ind_0ex(find(ind_0ex == min(ind_0ex_pos_i)))];
%                     w0s_access_i = [w0s(find(ind_0ex == min(ind_0ex_pos_i)))];
%                 elseif isempty(ind_0ex_pos_i)
%                     ind_0ex_access_i = [ind_0ex(find(ind_0ex == max(ind_0ex_neg_i)))];
%                     w0s_access_i = [w0s(find(ind_0ex == max(ind_0ex_neg_i)))];
%                 else
%                     ind_0ex_access_i = [ind_0ex(find(ind_0ex == max(ind_0ex_neg_i))) ind_0ex(find(ind_0ex == min(ind_0ex_pos_i)))];
%                     w0s_access_i = [w0s(find(ind_0ex == max(ind_0ex_neg_i))) w0s(find(ind_0ex == min(ind_0ex_pos_i)))];
%                 end
%                 
%                 ind_1ex_sub = ind_1ex - length(H_)./3;
%                 ind_dif_j = ind_0ex - ind_1ex_sub(j);
%                 ind_0ex_neg_j = [];
%                 ind_0ex_pos_j = [];
%                 ii = 0;
%                 jj = 0;
%                 for m = 1:length(ind_dif_j)
%                     if ind_dif_j(m) < 0
%                         ii = ii + 1;
%                         ind_0ex_neg_j(ii) = ind_0ex(m);
%                     elseif ind_dif_j(m) > 0
%                         jj = jj + 1;
%                         ind_0ex_pos_j(jj) = ind_0ex(m);
%                     end
%                 end
%                 
%                 if isempty(ind_0ex_neg_j)
%                     ind_0ex_access_j = [ind_0ex(find(ind_0ex == min(ind_0ex_pos_j)))];
%                     w0s_access_j = [w0s(find(ind_0ex == min(ind_0ex_pos_j)))];
%                 elseif isempty(ind_0ex_pos_j)
%                     ind_0ex_access_j = [ind_0ex(find(ind_0ex == max(ind_0ex_neg_j)))];
%                     w0s_access_j = [w0s(find(ind_0ex == max(ind_0ex_neg_j)))];
%                 else
%                     ind_0ex_access_j = [ind_0ex(find(ind_0ex == max(ind_0ex_neg_j))) ind_0ex(find(ind_0ex == min(ind_0ex_pos_j)))];
%                     w0s_access_j = [w0s(find(ind_0ex == max(ind_0ex_neg_j))) w0s(find(ind_0ex == min(ind_0ex_pos_j)))];
%                 end
%                 
%                 ind_0ex_shared = intersect(ind_0ex_access_i, ind_0ex_access_j);
%                 w0s_shared = intersect(w0s_access_i, w0s_access_j);
                
                for l = 1:n0
                    
                    [aa,mu0i] = unit_vector(mu1(i,:));
                    [bb,mu0j] = unit_vector(mu1(j,:));
                    [cc,mujl] = unit_vector(mu0(l,:,j));
                    [dd,muil] = unit_vector(mu0(l,:,i));
                    dipole = mu0i*mu0j*mujl*muil;
                    angle = polarizationInvariant(e_1,e_2,e_3,e_4,...
                        aa,bb,cc,dd);
                    
                    % non-rephasing diagram R4
                    fprintf(fid,'R4\t%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(i),ind_0ex(l),w1s(j)./(2*pi*wavenumbersToInvPs)+w0,(w1s(j)-w1s(i))./(2*pi*wavenumbersToInvPs),(w1s(j)-w0s(l))./(2*pi*wavenumbersToInvPs)+w0,-dipole*angle);
                    R_nr = R_nr - dipole*angle*exp(- 1i*w1s(j).*T1 ...
                        - 1i*(w1s(j)-w1s(i))*t2 ...
                        - 1i*(w1s(j)-w0s(l)).*T3);
                    
                    % non-rephasing diagram R5
                    fprintf(fid,'R5\t%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_0ex(l),ind_1ex(i),w1s(j)./(2*pi*wavenumbersToInvPs)+w0,w0s(l)./(2*pi*wavenumbersToInvPs),w1s(i)./(2*pi*wavenumbersToInvPs)+w0,-dipole*angle);
                    R_nr = R_nr - dipole*angle*exp(- 1i*w1s(j).*T1 ...
                        - 1i*w0s(l).*t2 ...
                        - 1i*w1s(i).*T3);
                end
                
                %             n2_start = 1;
                %             if length(w1s) > 1;
                %                 n2_start = i;
                %             end
                
                for k = 1:n2
                    %molecular dipoles
                    [cc,muik_] = unit_vector(mu2(k,:,i));
                    [dd,mujk_] = unit_vector(mu2(k,:,j));
                    dipole = mu0i*mu0j*muik_*mujk_;
                    angle = polarizationInvariant(e_1,e_2,e_3,e_4,...
                        aa,bb,cc,dd);
                    
                    %non-rephasing diagram R6
                    fprintf(fid,'R6\t%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(i),ind_2ex(k),w1s(j)./(2*pi*wavenumbersToInvPs)+w0,(w1s(j)-w1s(i))./(2*pi*wavenumbersToInvPs),(w2s(k)-w1s(i))./(2*pi*wavenumbersToInvPs)+w0,dipole*angle);
                    R_nr = R_nr + dipole*angle*exp(- 1i*w1s(j).*T1 ...
                        - 1i*(w1s(j)-w1s(i)).*t2 ...
                        - 1i*(w2s(k)-w1s(i)).*T3);
                end
            end
        end
    end
    % add lineshape (same for all peaks for now)
    R_nr = exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*R_nr;
    
    R_r_accum  = R_r_accum  + R_r *thermal_weight;
    R_nr_accum = R_nr_accum + R_nr*thermal_weight;
    J_accum = J_accum + J*thermal_weight;
    
    fprintf(fid,'\n\n');
end

J = J_accum;
R_r = R_r_accum;
R_nr = R_nr_accum;

if verbose>=1,disp('calculate fourier transforms...'),end

% calculate 1D spectrum (freq domain)
J = real(fftshift(sgrsifft(J,n_zp)));

% divide first points (by row and column) by 2
R_r(:,1) = R_r(:,1)./2;
R_r(1,:) = R_r(1,:)./2;
R_nr(:,1) = R_nr(:,1)./2;
R_nr(1,:) = R_nr(1,:)./2;

if flag_plot
    %what we have so far in the time domain
    figure(1),clf
    subplot(1,2,1)
    contourf(real(R_r'),10);
    axis equal tight
    subplot(1,2,2)
    contourf(real(R_nr'),10);
    axis equal tight
end

% do the fft
R_r = ifft2(R_r,n_zp,n_zp); %given the frequency definitions used
%above, use the ifft to get the
%frequencies right (Mathematica has the
%opposite definition of the fft by default)
R_nr = ifft2(R_nr,n_zp,n_zp);

%this is the frequency not the energy of the transitions
freq  = fftFreqAxis(t,'time_units','ps','zeropad',n_zp);
freq = freq+w0;

if flag_plot
    %now frequency domain
    figure(2),clf
    subplot(1,2,1)
    contourf(freq,freq,fftshift(real(R_r')),20); %pump-probe axis convention
    %contourf(fftshift(real(R_r)),20; % the (omega_1, omega_3) axis convention
    axis equal tight
    subplot(1,2,2)
    contourf(freq,freq,fftshift(real(R_nr')),20)
    %contourf(fftshift(real(R_nr)),20)
    axis equal tight
end

% flip R_r (being careful to keep zero frequency as the first time
% point), add the response functions, take the real part, and
% finally reorganize so that the 0 frequency is in the center
R = fftshift(real(fliplr(circshift(R_r,[0 -1]))+R_nr));


if flag_plot
    figure(3),clf
    n_contours = 40;
    MAX = max(abs(R(:)));
    level_list = linspace(-MAX,MAX,n_contours+2);
    dl = level_list(2)-level_list(1);
    cmin =level_list(1)-dl/2;
    cmax =level_list(end);
    
    %contourf(freq,freq,R',level_list) %use R' to display the pump-probe axis convention
    contourf(freq,freq,R,level_list) %use R to display the (omega_1, omega_3) axis convention
    caxis([cmin cmax]);
    axis equal tight
end

if verbose>=1,disp('package output...'),end

%package output
out.w1 = freq;
out.w3 = freq;
out.J = J;
out.R_r = R_r;
out.R_nr = R_nr;
out.R = R;
out.E = E;
out.V = V;
out.ind_0ex = ind_0ex;
out.ind_1ex = ind_1ex;
out.ind_2ex = ind_2ex;
out.energy_gap1 = w1s./(2*pi*wavenumbersToInvPs)+w0;
out.energy_gap2 = w2s./(2*pi*wavenumbersToInvPs)+2*w0;
out.mu1 = mu1;
out.mu2 = mu2;

end

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