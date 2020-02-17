function roVibHamiltonianMatrix = ...
    getRoVibHamiltonianMatrix_Vectorized(v1, v2, j1, j2, m1, m2, gasConst)

B_e = gasConst.B_e;
D = gasConst.D;
a_e = gasConst.a_e;
% v_e = sqrt(4*B_e^3/abs(D));
v_ex_e = gasConst.v_ex_e;
v_e = gasConst.v_e;

% finding where the vibrational and rotational levels are the same
v = v1 == v2;
J = j1 == j2;
M = m1 == m2;

% Creating matrix to make non-diagonal terms go to zero
overlap = v.*J.*M;

% Creating hamiltonian matrix
roVibHamiltonianMatrix = overlap.*(v_e .* (v1 + 1/2) ...
    - v_ex_e .* (v1 + 1/2).^2 ...
    + B_e .* j1 .* (j1+1) ...
    - a_e .* (v1+1/2).*j1.*(j1+1) ...
    - D .* j1.^2 * (j1+1).^2);
