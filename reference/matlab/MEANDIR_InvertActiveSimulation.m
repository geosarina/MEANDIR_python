function [X, EMIterInst0_updated, functioncost] = MEANDIR_InvertActiveSimulation(Solver, EMIterInst0, EMIterInst_r, RiverColumn0, RiverColumn_r, EMList0, Xdirect, MinFractionalContribution_r, MaxFractionalContribution_r, SolveCFList_r, abspos_r, relpos_r, WeightingList_r, fractionation, sources)

        % this function sets up the actual inversion of each simulation, which is 
        % performed as an optimization of "MEANDIR_CostFunction" using fmincon

        % (1) if the river data has a NaN or is empty, set the solution to be all NaNs
        if     sum(isnan(RiverColumn_r))>0 
               X = NaN(size(EMList0'));
               EMIterInst0_updated = EMIterInst0;
               functioncost = NaN;
               disp('unable to perform the inversion! "RiverColumn_r" contains a NaN value')
        elseif isempty(EMIterInst_r) 
               X = NaN(size(EMList0'));
               EMIterInst0_updated = EMIterInst0;
               functioncost = NaN;
               disp('unable to perform the inversion! "EMIterInst_r" is empty')
        else

            % (2) find the initial conditions. if there are fractionations, those values are  are recorded within EMIterInst_r,
            % so define initial_EMIterInst_r to have the end-member express half the fractionation offset from the sample value
            initial_EMIterInst_r = EMIterInst_r;                                                                                            % copy the reduced inversion matrix
            for ii=1:fractionation.red.n_r                                                                                                  % iterate over the requested number of fractionations
                ionvalue     = initial_EMIterInst_r(fractionation.red.ionpos_r(ii),fractionation.red.empos_r(ii));                          % find the normalized concentration ratio for this sample
                fracvalue    = initial_EMIterInst_r(fractionation.red.isopos_r(ii),fractionation.red.empos_r(ii))/ionvalue;                 % calculate the current fractionation factor
                riverisotope = RiverColumn_r(fractionation.red.isopos_r(ii))/RiverColumn_r(fractionation.red.ionpos_r(ii));                 % calculate the isotopic composition of the river sample
                initial_EMIterInst_r(fractionation.red.isopos_r(ii),fractionation.red.empos_r(ii)) = ionvalue*(riverisotope+0.5*fracvalue); % set the inital isotopic composition of the end-member to be the sample value plus half the fractionation
            end            
            if     isequal(Solver,'mldivide');           X0 = mldivide(initial_EMIterInst_r,RiverColumn_r);                        functioncost = sqrt(sum((RiverColumn_r - initial_EMIterInst_r*X0).^2)); % use the mldivide solver, solve for absolute error
            elseif isequal(Solver,'mldivide_optimize');  X0 = mldivide(initial_EMIterInst_r,RiverColumn_r);                        functioncost = sqrt(sum((RiverColumn_r - initial_EMIterInst_r*X0).^2)); % use the mldivide solver, solve for absolute error
            elseif isequal(Solver,'lsqnonneg');          X0 = lsqnonneg(initial_EMIterInst_r,RiverColumn_r);                       functioncost = sqrt(sum((RiverColumn_r - initial_EMIterInst_r*X0).^2)); % use the lsqnonneg solver, solve for absolute error
            elseif isequal(Solver,'lsqnonneg_optimize'); X0 = lsqnonneg(initial_EMIterInst_r,RiverColumn_r);                       functioncost = sqrt(sum((RiverColumn_r - initial_EMIterInst_r*X0).^2)); % use the lsqnonneg solver, solve for absolute error
            else   isequal(Solver,'optimize');           X0 = ones(size(initial_EMIterInst_r,2),1)./size(initial_EMIterInst_r,2);  functioncost = sqrt(sum((RiverColumn_r - initial_EMIterInst_r*X0).^2)); % make a vector of equal initial conditions, solve for absolute error
            end
            
            % if the functioncost is NaN, set functioncost equal to 99999
            if isnan(functioncost)
               functioncost = 999999;
            end
                       
            % (3) find an optimized solution by defining a cost function to be minimized, which is MEANDIR_CostFunction using the current EMIterInst_r and RiverColumn_r    
            if    (isequal(Solver,'optimize') | isequal(Solver,'mldivide_optimize') | isequal(Solver,'lsqnonneg_optimize')) && sum(isnan(X0))==0
                  options_fmincon = optimset('fmincon');                                                                                                                                       % set optimization parameters
                  options_fmincon = optimset('MaxFunEvals',1e3*length(EMList0),'Display','none');                                                                                              % set optimization parameters
                  f = @(X)MEANDIR_CostFunction(X,EMIterInst0, EMIterInst_r, RiverColumn0, RiverColumn_r, SolveCFList_r, abspos_r, relpos_r, WeightingList_r, fractionation, Xdirect, sources); % define a function
                  [Xtemp,functioncost] = fmincon(f,X0,[],[],[],[],MinFractionalContribution_r,MaxFractionalContribution_r,[],options_fmincon);                                                 % optimize that function
            else  
                  % do not run an optimization, set X to be the initial solution  
                  Xtemp = X0; 
            end

            % (4) given the inversion results, re-construct the matrix of fractional contributions to the normalization variable
            X = NaN(size(EMList0'));                           % define a vector of the correct final size
            if  sum(isnan(Xdirect))~=length(Xdirect)           % evaluate if Xdirect is entirely NaN values
                X(~isnan(Xdirect)) = Xdirect(~isnan(Xdirect)); % fill X with the directly chosen values
                X(isnan(Xdirect))  = Xtemp;                    % fill X with the values found through inversion
            else
                X = Xtemp; 
            end       

            % (5) account for any change to the inversion matrix due to fractionation
            EMIterInst0_updated = EMIterInst0;
            for i = 1:fractionation.all.n0
                if  isnan(fractionation.all.isopos_r(i)) | isnan(fractionation.all.empos_r(i))
                    % if a fractionation was requested but the relevant term was not included in the inversion, define the isotopic information as -1001
                    EMIterInst0_updated(fractionation.all.isopos0(i),fractionation.all.empos0(i)) = EMIterInst0_updated(fractionation.all.ionpos0(i),fractionation.all.empos0(i)).*(-1001);         % -1001 (which is impossible) indicates that the parameter was not used
                else
                    Xactive                = EMIterInst0(fractionation.all.ionpos0(i),:).*X'./RiverColumn0(fractionation.all.ionpos0(i));                                                           % convert the X values to fraction of the active ion (equation 14)
                    endmemberisotopevalues = EMIterInst0(fractionation.all.isopos0(i),:)./EMIterInst0(fractionation.all.ionpos0(i),:);                                                              % calculate the isotopic information of the end-members
                    sourceindx             = sources & Xactive~=0;                                                                                                                                  % find the position of sources with non-zero contributions
                    sourceisotopevalue     = sum(Xactive(sourceindx).*endmemberisotopevalues(sourceindx))./sum(Xactive(sourceindx));                                                                % calculate the isotopic information value of the inputs
                    activefractionation    = EMIterInst0(fractionation.all.isopos0(i),fractionation.all.empos0(i))./EMIterInst0(fractionation.all.ionpos0(i),fractionation.all.empos0(i));          % calculate the current fractionation
                    newisotopevalue        = (sourceisotopevalue + activefractionation) + activefractionation.*Xactive(fractionation.all.empos0(i))/sum(Xactive(sourceindx));                       % calculate the isotope information of the secondary phase
                    EMIterInst0_updated(fractionation.all.isopos0(i),fractionation.all.empos0(i)) = EMIterInst0_updated(fractionation.all.ionpos0(i),fractionation.all.empos0(i)).*newisotopevalue; % replace the value of the secondary phase in the inversion matrix
                end % end asking if the fractionation are NaN        
            end     % end iteration on fractionation.all.n0
        end         % end 'if' statement on whether there are NaNs in the river data or if the end-member matrix is empty
end                 % end the function