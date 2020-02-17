function [j1, j2] = getJStateMatrices(maxVibLvl, maxRotLvl)

totVibLvls = maxVibLvl + 1;
totRotLvls = (maxRotLvl + 1)^2;


j = zeros(1, totRotLvls*totVibLvls);

for ii = 1:maxRotLvl
    j = cat(1, j, ii.*ones(2*ii+1, totRotLvls*totVibLvls));
end

j1 = j;

for ii = 1:maxVibLvl
    j1 = cat(1, j1, j);
end

j2 = j1';