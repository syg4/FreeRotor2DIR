function [m1, m2] = getMStateMatrices(maxVibLvl, maxRotLvl)

totVibLvls = maxVibLvl + 1;
totRotLvls = (maxRotLvl + 1)^2;
m0 = zeros(1, totRotLvls*totVibLvls);


for ii = 1:maxRotLvl
    ri = meshgrid(-ii:ii, 1:(totRotLvls*totVibLvls))';
    m0 = cat(1, m0, ri);
end

m1 = m0;

for ii = 1:maxVibLvl
    m1 = cat(1, m1, m0);
end

m2 = m1';