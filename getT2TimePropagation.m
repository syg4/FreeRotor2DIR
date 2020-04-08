function t2TimeProp = getT2TimePropagation(t2, tau_J, maxVibLvl, maxRotLvl, temperature, symmetry)

k_B_SI = 1.3865e-23;
h = 6.6260e-34;
c = 2.9979e+10;

k = k_B_SI/(h*c);

[ji, jf] = getJStateMatrices(0, maxRotLvl);

j_even1 = zeros(length(ji));
j_even2 = zeros(length(jf));
j_even1(mod(ji,2)==0) = 1;
j_even2(mod(jf,2)==0) = 1;
J_even = j_even1.*j_even2;

j_odd1 = zeros(length(ji));
j_odd2 = zeros(length(jf));
j_odd1(mod(ji,2)~=0) = 1;
j_odd2(mod(jf,2)~=0) = 1;
J_odd = j_odd1.*j_odd2;

switch symmetry
    case 'even'
        J = J_even;
        if maxVibLvl > 0
            for ii = 1:maxVibLvl
                if mod(ii,2)==0
                    J = blkdiag(J, J_even);
                else
                    J = blkdiag(J, J_odd);
                end
            end
        end
    case 'odd'
        J = J_odd;
        if maxVibLvl > 0
            for ii = 1:maxVibLvl
                if mod(ii,2)==0
                    J = blkdiag(J, J_odd);
                else
                    J = blkdiag(J, J_even);
                end
            end
        end
    otherwise
        J = ones(maxRotLvl);
        if maxVibLvl > 0
            for ii = 1:maxVibLvl
                J = blkdiag(J,ones(maxRotLvl));
            end
        end
end


[j1, j2] = getJStateMatrices(maxVibLvl, maxRotLvl);

g = j1.*(j1>=j2)+j2.*(j1<j2);
g = 2.*g+1;

phi = exp((-1.*(0.4 .* j2 .* (j2+1)))./(k.*temperature));

phi = phi./sum(phi);

% vdiag = eq(v1, v2);

t2TimeProp = J.*(eq(j2,j1).*exp(-t2./tau_J) + phi.*(1-exp(-t2./tau_J)))./g;

t2TimeProp = sparse(t2TimeProp);
% J.*(2.*pi.*momentOfInertia.*k.*temperature.*(1-exp(-2.*t2./t_j))).*...
%     exp(-1.*(j2-j1.*exp(-1.*t2./t_j))^2./(2.*momentOfInertia.*k.*temperature.*(1-exp(2.*t2./t_j))));

% t2TimeProp.Conj = conj(t2TimeProp.Norm);