function timePropagationStruct = getTimePropagationCells(HamiltonianMatrix, ts, parpool_flag) %transDipoleMatrix,

h = 6.626E-34*5.034E22*10^15; %[cm^-1*fs]
hbar = h/(2*pi);
i = sqrt(-1);

timePropagationCell = cell(1, length(ts));
timePropagationCellConj = cell(1, length(ts));

if parpool_flag
    parfor ii = 1:length(ts)
            timePropagationMatrix = expm(-1.*i.*HamiltonianMatrix.*ts(ii)./hbar);
            timePropagationMatrix = sparse(timePropagationMatrix);

            timePropagationCell(ii) = {timePropagationMatrix};
            timePropagationCellConj(ii) = {conj(timePropagationMatrix)};
    end
else
    for ii = 1:length(ts)
            timePropagationMatrix = expm(-1.*i.*HamiltonianMatrix.*ts(ii)./hbar);
            timePropagationMatrix = sparse(timePropagationMatrix);

            timePropagationCell(ii) = {timePropagationMatrix};
            timePropagationCellConj(ii) = {conj(timePropagationMatrix)};
    end
end

timePropagationStruct.Norm = timePropagationCell;
timePropagationStruct.Conj = timePropagationCellConj;