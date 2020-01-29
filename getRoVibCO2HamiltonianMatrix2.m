function roVibCO2HamiltonianMatrix = getRoVibCO2HamiltonianMatrix2(maxVibLvl, maxRotLvl)

%     numMLvls = 0;
%     for ii = 0:maxRotLvl
%         numMLvls = numMLvls + 2*ii+1;
%     end
%     clear ii

%     if isunix
%         path = '/Users/kaigronborg/OneDrive/Documents/Graduate Research/Projects/THFFreeRotor/Simulations/Data/CO2_THF_950um_FlowCell.CSV';
%     else
%         path = 'C:\Users\Kai\OneDrive\Documents\Graduate Research\Projects\THFFreeRotor\Simulations\Data\CO2_THF_950um_FlowCell.CSV';
%     end
    path = 'C:\Users\kaigr\OneDrive\Documents\Graduate Research\MATLAB\FreeRotor2DIR\Data\CO2_THF_950um_FlowCell.CSV';
    
    numMLvls = (maxRotLvl+1)^2;
    
    roVibCO2HamiltonianMatrix = zeros((maxVibLvl+1)*numMLvls);
    
    [~, ~, gasConst] = getEmpericalCO2GasPhaseHamiltonian2(path);

    close all
    
%     H(2,0, gasConst)
%     
%     H(1, 0, gasConst) - H(0, 0, gasConst)
%     
%     H(2, 0, gasConst) - H(1, 0, gasConst)
    
    index = 0;
    for jj = 0:maxVibLvl
        for ii = 0:maxRotLvl
            
            for kk = 1:(2*(ii)+1)
                ind = kk + index;
                roVibCO2HamiltonianMatrix(ind, ind) = H(jj, ii, gasConst);
            end
            index = index + 2*(ii)+1;
            clear kk
        end
    end
    
    
    
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
        