function [v1, v2] = getVibrationalStateMatrices(maxVibLvl, maxRotLvl)

totRoLvls = (maxRotLvl+1)^2;
totVibLvls = maxVibLvl + 1;

v1 = zeros(totRoLvls, totVibLvls*totRoLvls);

for ii = 1:maxVibLvl
    v1 = cat(1, v1, ii.*ones(totRoLvls, totVibLvls*totRoLvls));
end
v2 = v1';
