function J4 = getJ4(densityMatrix, transitionDipoleMatrix, t1TimeProp)

rho0 = cell(1, length(t1TimeProp.Norm)); 
rho0(:) = {transitionDipoleMatrix*densityMatrix};

transitionDipoleCell1 = cell(1, length(t1TimeProp.Norm));
transitionDipoleCell1(:) = {transitionDipoleMatrix};


rho1 = cellfun(@mtimes, rho0, t1TimeProp.Conj, 'UniformOutput', false);
rho1 = cellfun(@mtimes, t1TimeProp.Norm, rho1, 'UniformOutput', false);
rho1 = cellfun(@mtimes, rho1, transitionDipoleCell1, 'UniformOutput', false);

J4 = cellfun(@trace, rho1);