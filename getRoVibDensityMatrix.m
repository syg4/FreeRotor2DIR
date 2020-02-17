function roVibDensityMatrix = getRoVibDensityMatrix(v1, v2, j1, j2, m1, m2, ...
    roVibHamiltonianMatrix, temperature, symmetry)

k_B_SI = 1.3865e-23;
h = 6.6260e-34;
c = 2.9979e+10;

k = k_B_SI/(h*c);

v = v1 == v2;
J = j1 == j2;
M = m1 == m2;
J1 = j1.*J.*M.*v;

switch symmetry
    case 'even'
        symMatrix = diag(diag(mod(J1, 2) == 0));
        
    case 'odd'
        symMatrix = diag(diag(mod(J1, 2) == 1));
        
    case 'all'
        symMatrix = v.*J.*M;
    otherwise
        error('Symmetry must be either "even", "odd", or "all".')
end


roVibDensityMatrix = symMatrix.*exp((-1.*(roVibHamiltonianMatrix))./(k.*temperature));

densityTrace = trace(roVibDensityMatrix);

roVibDensityMatrix = roVibDensityMatrix./densityTrace;

