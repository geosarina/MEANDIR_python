% define matricies and cells to hold the results of the inversion
smprescell                   = cell(s,1);                 % will hold results for the fraction contribution of end-members to the normalization variable
InvRivDatacell               = cell(s,1);                 % will hold results for the river chemistry inverted on each simulation
RiverColumnSucessMatrixcell  = cell(s,numberiterations);  % will hold results for the normalized river chemistry inverted on each simulation
EMSucessIter0Matrixcell      = cell(s,numberiterations);  % will hold results for the end-member chemistry of each simulation
EMSucessIter0Matrixcell_frac = cell(s,numberiterations);  % will hold results for the end-member chemistry of each simulation containing fractionating
iterationcountercell         = cell(s,numberiterations);  % cell to count number of iterations attempted
successcountercell           = cell(s,numberiterations);  % cell to count the number of successfull simulations
badmbcountercell             = cell(s,numberiterations);  % cell to count number of simulations out of mass balance
smpresmisfitcell             = cell(s,numberiterations);  % misfit of observations and model results 
smpfunmisfitcell             = cell(s,numberiterations);  % misfit of cost function

iterationcounts              = NaN(s,numberiterations);   % vector to count number of iterations attempted
successfullsimulations       = NaN(s,numberiterations);   % vector to count the number of successfull simulations
badmbcount                   = NaN(s,numberiterations);   % vector to count number of simulations out of mass balance

fractionation                = NaN;                       % to get defined in the workspace