% getCO2RoVibTransitionDipoleMatrix2   This is a function that creates a rovibrational transition dipole matrix
%   for CO_2. This function can also return a creation matrix and annihilation
%   matrix.

    
function [transDipoleMatrixZ, transDipoleMatrixX, transDipoleMatrixY, Cz, Az, Cx, Ax, Cy, Ay] ...
    = getCO2RoVibTransDipoleMatrix_Vectorized(maxVibLvl, maxRotLvl)

%     numMLvls = (maxRotLvl + 1)^2;

    
    [v1, v2] = getVibrationalStateMatrices(maxVibLvl, maxRotLvl);

    [j1, j2] = getJStateMatrices(maxVibLvl, maxRotLvl);
    
    [m1, m2] = getMStateMatrices(maxVibLvl, maxRotLvl);
    
    transDipoleMatrixZ = arrayfun(@getRoVibTransDipoleMomentZ, v1, v2, j1, j2, m1, m2);
    transDipoleMatrixX = arrayfun(@getRoVibTransDipoleMomentX, v1, v2, j1, j2, m1, m2);
    transDipoleMatrixY = arrayfun(@getRoVibTransDipoleMomentY, v1, v2, j1, j2, m1, m2);
    
    
    
%     transDipoleMatrixZ = zeros((maxVibLvl+1)*numMLvls);
% %     Cz = zeros((maxVibLvl+1)*numMLvls);
% %     Az = zeros((maxVibLvl+1)*numMLvls);
%     
%         
%     transDipoleMatrixX = zeros((maxVibLvl+1)*numMLvls);
% %     Cx = zeros((maxVibLvl+1)*numMLvls);
% %     Ax = zeros((maxVibLvl+1)*numMLvls);
%         
%     transDipoleMatrixY = complex(zeros((maxVibLvl+1)*numMLvls));
% %     Cy = complex(zeros((maxVibLvl+1)*numMLvls));
% %     Ay = complex(zeros((maxVibLvl+1)*numMLvls));
    
    
%     rowInd = 1; % row index
%     colInd = 1; % column index
%     
%     for v_final = 0:maxVibLvl
% %         fprintf('----------   Starting v = %i Transitions   ----------\n', v_final);
%         for J_final = 0:maxRotLvl
%             for m_final = -J_final:J_final
%                 for v_init = 0:maxVibLvl
%                     for J_init = 0:maxRotLvl
%                         for m_init = -J_init:J_init
%                             
%                             transDipoleMatrixZ(rowInd, colInd) = getRoVibTransDipoleMomentZ(v_init, v_final, J_init, J_final, m_init, m_final);
%                             transDipoleMatrixX(rowInd, colInd) = getRoVibTransDipoleMomentX(v_init, v_final, J_init, J_final, m_init, m_final);
%                             transDipoleMatrixY(rowInd, colInd) = getRoVibTransDipoleMomentY(v_init, v_final, J_init, J_final, m_init, m_final);
%                             
%                             colInd = colInd+1;
%                         end
%                     end
%                 end
%                 rowInd = rowInd+1;
%                 colInd = 1;
%             end
%         end
%     end
    
    Cz = triu(transDipoleMatrixZ);
    Cx = triu(transDipoleMatrixX);
    Cy = triu(transDipoleMatrixY);
%     for ii = 1:length(transDipoleMatrixZ)
%         for jj = ii:length(transDipoleMatrixZ)
%             Cz(ii, jj) = transDipoleMatrixZ(ii, jj);
%             Cx(ii, jj) = transDipoleMatrixX(ii, jj);
%             Cy(ii, jj) = transDipoleMatrixY(ii, jj);
%         end
%     end
%     clear ii jj

    Az = tril(transDipoleMatrixZ);
    Ax = tril(transDipoleMatrixX);
    Ay = tril(transDipoleMatrixY);
    
%     for ii = 1:length(transDipoleMatrixZ)
%         for jj = ii:length(transDipoleMatrixZ)
%             Az(jj, ii) = transDipoleMatrixZ(jj, ii);
%             Ax(jj, ii) = transDipoleMatrixX(jj, ii);
%             Ay(jj, ii) = transDipoleMatrixY(jj, ii);
%         end
%     end    
end
