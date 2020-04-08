% getCO2RoVibTransitionDipoleMatrix2   This is a function that creates a rovibrational transition dipole matrix
%   for CO_2. This function can also return a creation matrix and annihilation
%   matrix.


function varargout = getRoVibTransDipoleMatrix(v1, v2, j1, j2, m1, m2, whichDims)

nOutputs = 3*length(whichDims);
varargout = cell(1,nOutputs);

ii = 1;

for jj = 1:length(whichDims)
    % generate transition dipole matrices
    if strcmp(whichDims(jj), 'z')
        transDipoleMatrixZ = arrayfun(@getRoVibTransDipoleMomentZ, v1, v2, j1, j2, m1, m2);
        % generate creation matrices
        Cz = triu(transDipoleMatrixZ);
        % generate annihilation matrices
        Az = tril(transDipoleMatrixZ);

        varargout{ii} = transDipoleMatrixZ;
        varargout{ii+1} = Cz;
        varargout{ii+2} = Az;
        ii = ii+3;
    elseif strcmp(whichDims(jj), 'x')
        transDipoleMatrixX = arrayfun(@getRoVibTransDipoleMomentX, v1, v2, j1, j2, m1, m2);
        Cx = triu(transDipoleMatrixX);
        Ax = tril(transDipoleMatrixX);
        
        varargout{ii} = transDipoleMatrixX;
        varargout{ii+1} = Cx;
        varargout{ii+2} = Ax;
        ii = ii+3;
    elseif strcmp(whichDims(jj), 'y')
        transDipoleMatrixY = arrayfun(@getRoVibTransDipoleMomentY, v1, v2, j1, j2, m1, m2);
        Cy = triu(transDipoleMatrixY);
        Ay = tril(transDipoleMatrixY);
        
        varargout{ii} = transDipoleMatrixY;
        varargout{ii+1} = Cy;
        varargout{ii+2} = Ay;
        ii = ii+3;
    end
end

end
