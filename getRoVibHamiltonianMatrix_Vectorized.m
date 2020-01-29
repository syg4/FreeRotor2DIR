function roVibHamiltonianMatrix = getRoVibHamiltonianMatrix_Vectorized(maxVibLvl, maxRotLvl, gasConsts)
    
    
    
    
    function energy = H(v, J, gasConst)
        
        B_e = gasConst.B_e;
        D = gasConst.D;
        a_e = gasConst.a_e;
%         v_e = sqrt(4*B_e^3/abs(D));
        v_ex_e = gasConst.v_ex_e;
        v_e = gasConst.v_e;
        
        
        energy = v_e * (v + 1/2) - v_ex_e * (v + 1/2)^2 + B_e * J * (J+1) - a_e * (v+1/2)*J*(J+1) - D * J^2 * (J+1)^2;
    end
end
        