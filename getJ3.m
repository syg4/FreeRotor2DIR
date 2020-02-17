function J3 = getJ3(densityMatrix, transitionDipoleMatrix, t1TimeProp)

rho0 = cell(1, length(t1TimeProp.Norm)); 
rho0(:) = {densityMatrix*transitionDipoleMatrix};

transitionDipoleCell1 = cell(1, length(t1TimeProp.Norm));
transitionDipoleCell1(:) = {transitionDipoleMatrix};


rho1 = cellfun(@mtimes, rho0, t1TimeProp.Conj, 'UniformOutput', false);
rho1 = cellfun(@mtimes, t1TimeProp.Norm, rho1, 'UniformOutput', false);
rho1 = cellfun(@mtimes, transitionDipoleCell1, rho1, 'UniformOutput', false);

J3 = cellfun(@trace, rho1);