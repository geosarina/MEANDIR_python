%% MIXING ELEMENTS AND DISSOLVED ISOTOPES IN RIVERS (MEANDIR)
% version 1.3: last edited September 6th, 2023
% written by Preston Cosslett Kemeny, University of Chicago
% pkemeny@uchicago.edu, preston.kemeny@gmail.com
%
% this code is extremely similar to version 1.1, which underlies the manuscript:
% "PRESENTATION AND APPLICATIONS OF MIXING ELEMENTS AND DISSOLVED ISOTOPES IN RIVERS (MEANDIR), 
% A CUSTOMIZABLE MATLAB MODEL FOR MONTE CARLO INVERSION OF DISSOLVED RIVER CHEMISTRY"
%  Preston Cosslett Kemeny and Mark Albert Torres, published in the American Journal of Science (2021)
%
% (1) define the name of a scenario in ScenarioList
% (2) define the parameters of that scenario in MEANDIR_FindScenarioParameters
% (3) add river observations into a spreadsheet, read data with MEANDIR_ReadRiverDataIntoMatlab
% (4) define end-member in spreadsheet, read them with MEANDIR_ReadEndMembersIntoMatlab
% (5) inversion results are saved as fields of the river structure

%% SECTION 1: DEFINE INVERSION
clearvars; close all; clc; % clear all variables, close all figure, and clear the command line
overwritesimnumbers = 1;   % if 1, the number of simulations defined in MEANDIR_FindScenarioParameters will be overwritten
                           % to a value specified below. In general, overwritesimnumbers should equal 0 unless troubleshooting

% default: define a scenario called "MyDefaultScenario" in ScenarioList. To see the parameters associated with this
% scenario, open "MEANDIR_FindScenarioParameters". if only running one scenario, it is easiest to leave this variable as is
%ScenarioList = {'AK_scenario1_carbonate_slctindi'};
%ScenarioList = {'AK_scenario2_carbonate_slctindi_degas_2p5'};
%ScenarioList = {'AK_scenario3_carbonate_slctindi_degas_25'};
%ScenarioList = {'AK_scenario4_carbonate_slctindi_degas_2p5_highfrac'};
%ScenarioList = {'AK_scenario5_carbonate_slctindi_degas_2p5_lowfrac'};
ScenarioList = {'AK_scenario1_carbonate_slctindi','AK_scenario2_carbonate_slctindi_degas_2p5','AK_scenario3_carbonate_slctindi_degas_25','AK_scenario4_carbonate_slctindi_degas_2p5_highfrac','AK_scenario5_carbonate_slctindi_degas_2p5_lowfrac'};


%% SECTION 2: CHECK THAT INVERSION PARAMETERS DO NOT PRODUCE ERRORS
disp('MEANDIR initialized! evaluating requested inversion(s) for potential problems');
if exist('ScenarioList','var'); else; disp('no inversion(s) requested! add an entry to "ScenarioList" to begin using MEANDIR'); return; end
if isempty(ScenarioList);             disp('"ScenarioList" is empty, no inversion(s) requested! add an entry to "ScenarioList" to begin using MEANDIR'); end
for ScenarioListIndex = 1:length(ScenarioList) % iteration over the entries of ScenarioList
    printtxt = 0;                              % do not print the details of each iteration    
    MEANDIR_FindScenarioParameters;            % test that this scenario can be found and does not return errors
    if FailCase==1; disp(sprintf('catastrophic failure! MEANDIR shutting down due to problem with inversion: "%s"',ScenarioList{ScenarioListIndex}));
                    disp('if unable to troubleshoot, contact Preston Cosslett Kemeny at preston.kemeny@gmail.com for assistance'); return; end; end     
if     length(ScenarioList)==1; disp('scenario located! now beginning the inversion');  % tell user if all inversions are found
elseif length(ScenarioList)>1;  disp('scenarios located! now beginning inversions');    % tell user if all inversions are found
end

% the following code attempts to open resources for parallel processing
ncores = feature('numcores');
if ncores > 1
    setenv('TZ','UTC');  % set the time zone to avoid a warning
    v = ver; if  sum(ismember(cellstr(char(v.Name)),'Parallel Computing Toolbox'))>0; p = gcp('nocreate'); else; p = NaN; end % open a pool without activating cores. To delete parpool, enter: 'delete(gcp)'
    if  isempty(p)
        parpool(ncores-1);
    end
end

%% SECTION 3: ITERATE OVER EACH SCENARIO
savenames = {};                                  % declare an empty structure for saving names
for ScenarioListIndex = 1:length(ScenarioList)   % iterate over the entries of ScenarioList. 
    
    tic; % start a timer      
    % clear all variables in the workspace except those defined for the simulation
    clearvars -except ScenarioList ScenarioListIndex overwritesimnumbers savenames  
    
    % display the active scenario name
    disp(sprintf('\nSCENARIO (%i/%i): %s ',ScenarioListIndex,length(ScenarioList),ScenarioList{ScenarioListIndex}));
           
    % Load the data for this scenario and identify basic information
    printtxt = 1;                                       % print the information for this scenario
    MEANDIR_FindScenarioParameters;                     % find the parameters for this scenario
    MEANDIR_ReadRiverDataIntoMatlab;                    % read river observations based on the value of Riverdatasource
    EM              = MEANDIR_ReadEndMembersIntoMatlab; % read end-members chemical compositions
    ConvertDelta2R  = MEANDIR_ReadConvertDelta2R;       % read information on converting delta notation to isotopic ratios
    go              = MEANDIR_Go(ObsList);              % identify which variables are in the inversion
    nEM             = length(EMList0);                  % identify the number of end-members
    nOL             = length(ObsList);                  % udentify the number of observations
    ObsList         = ObsList';                         % transpose the observations vector
    WeightingList   = WeightingList';                   % transpose the weighting vector
    CostFunType     = CostFunType';                     % transpose the cost function vector    
    abspos          = ismember(CostFunType,'abs');      % find the index vector for absolute cost function evaluation
    relpos          = ismember(CostFunType,'rel');      % find the index vector for relative cost function evaluation
    isotopeposition = (ismember(ObsList,'d34S')   | ismember(ObsList,'d18O')   | ismember(ObsList,'Sr8786') | ...
                       ismember(ObsList,'d7Li')   | ismember(ObsList,'d26Mg')  | ismember(ObsList,'d30Si')  | ...
                       ismember(ObsList,'d42Ca')  | ismember(ObsList,'d44Ca')  | ismember(ObsList,'d56Fe')  | ...
                       ismember(ObsList,'d98Mo')  | ismember(ObsList,'Os8788') | ismember(ObsList,'d13C')  | ismember(ObsList,'Fmod'));
    indxisopos      = find(isotopeposition);           % Find the index vector for isotope data        
    
    % If overwritesimnumbers is 1, MEANDIR will overwrite the number of simulations 
    % when IterateOver equals 'Samples', numberiterations should equal 1
    % when IterateOver equals 'End-members', maxsuccess and maxiterations should equal 1    
    if isequal(IterateOver,'Samples')
        numberiterations = 1; 
        if overwritesimnumbers == 1          
           maxiterations = 400000; % attempt up to 100 simulations
           maxsuccess    = 250;    % target finding 10 successful simulations           
           maxzerohits   = 10000;  % limit attempts to 100 without success           
           fprintf('overwriting the maximum number of iterations to be %se%s and the maximum numbers of successes to be %s\n',num2str(floor(maxiterations/(10^floor(log10(maxiterations))))),num2str(floor(log10(maxiterations))),num2str(maxsuccess));
        end
        resultssize = maxsuccess; 
    end
    if isequal(IterateOver,'End-members')
        maxsuccess    = 1; 
        maxiterations = 1; 
        if overwritesimnumbers == 1
           numberiterations = 60; % run 60 simulations
           MisfitCuts       = 25; % keep the 25% of simulations with lowest misfit
           fprintf('overwriting the number of iterations to be %s and misfit cuts fraction to be %s%s\n',num2str(numberiterations),num2str(MisfitCuts),'%'); 
        end
        resultssize = numberiterations;
    end        
    
    % check if the optimization toolbox is required for the chosen solver.
    % if so, test if the user has the optimization toolbox
    if  isequal(Solver,'optimize') | isequal(Solver,'mldivide_optimize') | isequal(Solver,'lsqnonneg_optimize')
         v = ver;      
         if  sum(ismember(cellstr(char(v.Name)),'Optimization Toolbox'))==0
             disp(sprintf('warning! you are attempting to use the solver "%s", but MEANDIR cannot find the Optimization Toolbox installed in this MATLAB environment. install the toolbox or use "mldivide" or "lsqnonneg" instead',Solver))
         end
    else
         options_fmincon = NaN; 
    end    
    % if more than 1 core is found, check whether the parallel computing toolbox is installed in this environment
    ncores = feature('numcores');
    if ncores > 1
        v = ver;
        if  sum(ismember(cellstr(char(v.Name)),'Parallel Computing Toolbox'))==0
             disp('warning! MEANDIR finds multiple cores but cannot find the Parallel Computing Toolbox installed in this MATLAB environment. ignore this warning if not intending to run in parallel')
        end    
    end
    
    % generate matrix of river observations
    RiverMatrix0 = MEANDIR_GenerateRiverMatrix(river, nOL, ObsList, indxisopos, carbonisotopematch, ConvertDelta2RList, ConvertDelta2R);
    [l, s]       = size(RiverMatrix0);
    
    % generate distributions of end-member chemistry
    [OrEMdist, TrEMdist, distEMdist0, distEMdist0min, distEMdist0max, smpadddist_mean, smpadddist_std, smpadddist_min, smpadddist_max, lgu_sign] =  ...
    MEANDIR_makeEMdistributions(nOL, ObsList, nEM, EMList0, EM, EMdatasource, EndMembersWithNegativeRatios, CoupleFeS2SO4intoEM, NormalizationType, ObsInNormalization, IterateOver, indxisopos, ConvertDelta2RList, ConvertDelta2R);

    % generate structures to hold results and declare several variables as NaN 
    MEANDIR_DeclareHoldingVariables; % Pre-allocate a few variables used to hold results
    if isequal(IterateOver,'End-members')        
          MEANDIR_DeclareExternalLoopIterEMStructures 
        % MEANDIR_DeclareExternalLoopIterEMStructures defines IterEM_smpfracStore* w/ size (nEM, numberiterations,s)
        % MEANDIR_DeclareExternalLoopIterEMStructures defines IterEM_InvRivDataStore_* w/ size (s,numberiterations)
    end
 

    %% SECTION 4: PERFORM EACH INVERSION           
    % provide a warning to the user that an error may appear
    disp('when running in parallel you may see two warnings for the temporary variables "IterEM_smpfrac" and "IterEM_InvRivData". these can be ignored.'); 
    
    % only one of the two lines 176 or 177 should be active, and the other commented
        % when IterateOver equals 'Samples', this loop should be a normal loop (uncomment line 176)
        % when IterateOver equals 'End-members', this can be a parallel loop (uncomment line 177), but does not have to be
    for iterindx = 1:numberiterations
    %parfor iterindx = 1:numberiterations

        % if IterateOver equals 'End-members', display 20 progress markers during the run
        % to display inversion progress more/less regularly, increase/decrease the value of 'numberdivisions'
        if isequal(IterateOver,'End-members')
            numberdivisions = 20;
            step = floor(numberiterations/numberdivisions);
            if step==0; step = 1; end 
            ldisp = 1:step:numberiterations;
            if ldisp(end)~=numberiterations;    ldisp(end+1) = numberiterations; end
            if sum(ismember(ldisp,iterindx))>0; fprintf('reached simulation %i/%i\n',iterindx,numberiterations); end
        end     

        % for each simulation, find an instance of end-member chemistry
        % if IterateOver equals 'End-members', this is the matrix of end-member chemistry used in the inversion.
        % if IterateOver equals 'Samples', here define the matrix as NaN because it will be re-defined independently for each sample below.
        EMIterInst0_p = NaN;
        if isequal(IterateOver,'End-members')
           iterationcounter = 0;
           RiverColumn0     = 0;
           [EMIterInst0_p, FailCase] = MEANDIR_PullEndMemberRatios(EMList0, EM, EMdatasource, EndMembersWithNegativeRatios, ObsList, IonCharges, isotopeposition,indxisopos, NaN, iterationcounter, maxiterations, ...
                                       RiverColumn0, NormalizationType, ObsInNormalization, CationsInNormalization, AnionsInNormalization, NeutralInNormalization, ...
                                       ListNormClosure, ListChargeClosure, CoupleFeS2d34SintoEM, AllIonsExplicitlyResolved, BalanceEvaporite, ...
                                       carbonisotopematch, TrEMdist, distEMdist0, distEMdist0min, distEMdist0max, smpadddist_mean, smpadddist_std, smpadddist_min, smpadddist_max, ConvertDelta2RList, ConvertDelta2R, lgu_sign);
           if FailCase == 1; disp('MEANDIR was unable to find end-members!'); end  
        end

        % when IterateOver equals 'End-members', the structures IterEM_smpfrac and IterEM_InvRivData are used to move results from iterating over
        % samples within a for loop to the results of the overall parallel loop. as a result, they are re-set on each value of iterindx.       
        if isequal(IterateOver,'End-members')
           [IterEM_smpfrac, IterEM_InvRivData] = MEANDIR_DeclareInternalLoopIterEMStructures(go, nEM, s);
                % IterEM_smpfrac subfields have size (nEM,s)
                % IterEM_InvRivData subfields have size (s,1)
        end
       
        % only one of the two lines 215 or 216 should be active, and the other commented
            % when IterateOver equals 'End-members', this should be a normal for loop  (uncomment line 215)
            % when IterateOver equals 'Samples', this can be a parallel loop (uncomment line 216), but does not have to be
        %for i=1:s
        parfor i=1:s

            % if IterateOver equals 'Samples', display the sample currently being inverted. 
            if isequal(IterateOver,'Samples')
               fprintf('analyzing sample %i/%i\n',i,s);
            end

            % define variables for inversion and for holding inversion results
            Xtemp                        = NaN;                     % Xtemp is passed to the solver functions and must be defined
            workingresult                = NaN;                     % workingresult is passed to the solver functions and must be defined
            successcounter               = 0;                       % counter for successfull simulations
            badmbcounter                 = 0;                       % counter for the number of physcially realizable solutions out of mass balance
            iterationcounter             = 0;                       % counter for the number of iterations attempted
            EMSucessIter0Matrix          = NaN(nOL,nEM,maxsuccess); % end-member composition in successful simulations (will hold final inversion information)
            EMSucessIter0Matrix_frac     = NaN(nOL,nEM,maxsuccess); % end-member composition in successful simulations (will hold inverison information with fractionation factors)
            RiverColumnSucessMatrix      = NaN(maxsuccess,nOL);     % river chemistry in successful simulations
            smpresmisfit                 = NaN(1,maxsuccess);       % smpresmisfit has size (1,maxsuccess)
            smpfunmisfit                 = NaN(1,maxsuccess);       % smpfunmisfit has size (1,maxsuccess)
                      
            % define matricies for storing simulation results. 
            [InvRivDatsave, InvRivDat, smpres] = MEANDIR_DeclareSimulationStorageMatricies(maxsuccess, nEM, go);
                % InvRivDatsave.* has size (maxsuccess, 1)
                % InvRivDat.* has size (1,1)
                % smpres.smpresults* has size (nEM, maxsuccess)
                
            % set the while loop to only begin if all the data are present for this sample
            if   functionalsamplelist(i)==1
                 stopiterations = 0;
            else
                 stopiterations = 1;
            end
            
            % the following while statement is the actual Monte Carlo inversion. the while statement
            % is satisfied when the iterationcounter >= maxiterations or when the successcounter == maxsuccess
            while stopiterations == 0
                  iterationcounter = iterationcounter+1;  % increase the iteration counter by 1
                  Xdirect = NaN(nEM,1);                   % define Xdirect, which which hold some of the inversion results

                  % (4.1) update river data with error reflecting measurement uncertainty
                  RiverColumn0 = RiverMatrix0(:,i); % start the river column as the entry of RiverMatrix0
                  [InvRivDat, RiverColumn0] = MEANDIR_AdjustRiverDataForAnalyticalError(i,river,ObsList,RiverColumn0,AdjustRiverObs,InvRivDat,ObsInNormalization, carbonisotopematch, ConvertDelta2RList, ConvertDelta2R, go);
                  
                  % (4.2) generate the matrix of end-members for this simulation
                  EMIterInst0 = EMIterInst0_p;       % update EMIterInst0 to equal its preliminary value 
                                                     % if IterateOver equals 'Samples', at this point EMIterInst0 is NaN and will be set in the following lines
                                                     % if IterateOver equals 'End-members', at this point EMIterInst0 is a full end-member matrix
                  if  isequal(IterateOver,'Samples') % if iterating over samples, re-draw the end-member matrix
                      [EMIterInst0, FailCase] = MEANDIR_PullEndMemberRatios(EMList0, EM, EMdatasource, EndMembersWithNegativeRatios, ObsList, IonCharges, isotopeposition,indxisopos, i, iterationcounter, maxiterations, ...
                                                RiverColumn0, NormalizationType, ObsInNormalization, CationsInNormalization, AnionsInNormalization, NeutralInNormalization, ...
                                                ListNormClosure, ListChargeClosure, CoupleFeS2d34SintoEM, AllIonsExplicitlyResolved, BalanceEvaporite, ...
                                                carbonisotopematch, TrEMdist, distEMdist0, distEMdist0min, distEMdist0max, smpadddist_mean, smpadddist_std, smpadddist_min, smpadddist_max, ConvertDelta2RList, ConvertDelta2R, lgu_sign);
                                                if FailCase == 1; disp('MEANDIR was unable to find end-members! error in "MEANDIR_Master" for "MEANDIR_PullEndMemberRatios".'); end
                  end
                  
                  % (4.3) adjust river observations for Cl Critical values
                  if   isequal(PrecProcessing,'ClCrit') && sum(ismember(EMList0,'prec'))>0 && sum(ismember(ObsList,'Cl'))>0
                       [RiverColumn, EMIterInst, EMIterInst0, Xdirect, EMList, distEMdist, iterationcounter] = ...
                       MEANDIR_ClCriticalCorrection(river, InvRivDat, EMIterInst0, Xdirect, IterateOver, ClCriticalValuesGiven, nOL, EMList0, EM, EMdatasource, EndMembersWithNegativeRatios,...
                       ObsList, IonCharges, isotopeposition,indxisopos, i, iterationcounter, maxiterations, RiverColumn0, NormalizationType, ObsInNormalization, CationsInNormalization, AnionsInNormalization, NeutralInNormalization,...
                       ListNormClosure, ListChargeClosure, CoupleFeS2d34SintoEM, AllIonsExplicitlyResolved, BalanceEvaporite, carbonisotopematch, TrEMdist, distEMdist0, distEMdist0min, distEMdist0max, ...
                       smpadddist_mean, smpadddist_std, smpadddist_min, smpadddist_max, ConvertDelta2RList, ConvertDelta2R, lgu_sign, MinFractionalContribution);
                  else % if the PrecProcessing is not ClCrit or if precipitation is not included in the list of end-members
                       RiverColumn = RiverColumn0; % RiverColumn0 has been adjusted for analytical error
                       EMIterInst  = EMIterInst0;  % the EM chemistry matrix
                       EMList      = EMList0;      % the EMlist to invert
                       distEMdist  = distEMdist0;  % the distEMdist to invert
                       fNorm_rain  = NaN;          % the fraction of the normalization due to rainwater
                  end

                  % at this point in the code:
                  % RiverColumn contains the river data to be inverted (differs from RiverColumn0 when accounting for ClCritical values, where RiverColumn0 has been updated to reflect analytical error)
                  % EMIterInst contains the end-member data to be inverted (differs from EMIterInst0 when accounting for ClCritical values, where EMIterInst0 has been updated to make the ClCritical removal physically possible)
                  % EMList contains the list of end-members to be accounted for (differs from EMList0 when accounting for ClCritical values)

                  % (4.4) remove zeros from river data and EM matrix when evaluating equations for relative misfit 
                  [RiverColumn_r, EMIterInst_r, EMList_r, ObsList_r, Xdirect, distEMdist_r] = MEANDIR_PrepareUpdatedDataForInversion(Xdirect,RiverColumn,ObsList,EMList,EMList0,relpos,EMIterInst,Solver,distEMdist);
                  % RiverColumn_r and EMIterInst_r are the actual matricies that will be inverted. however, by this point in the code
                  % the contribution of precipitation may already be defined in Xdirect. similarly, based on the relative structure of
                  % the river observations and end-member chemistry, end-members known to have no fractional contribution may also have zero values in Xdirect. 

                  % (4.5) in the case of fractionation, find the pairs of isotope, ions, and end-members. the structure fractionation is given
                  % to the cost function so that these same calculations only need to be done here and not on each loop of the optimization          
                  fractionation = MEANDIR_FindFractionationPairs(distEMdist0, distEMdist_r, ObsList, ObsList_r, EMList0, EMList_r, carbonisotopematch, indxisopos);
    

                 % If requested, constrain the amount of degassing to be a fraction of the amount of DIC
                  if  ResetDegasDICContribution == 1 && sum(ismember(EMList0,'degas'))>0                    
                      [MinFractionalContribution, MaxFractionalContribution] = MEANDIR_resetDICcont(ObsList, EMList0, RiverColumn, DegasDICContributionMin, DegasDICContributionMax, MinFractionalContribution0, MaxFractionalContribution0);
                  else
                      MinFractionalContribution = MinFractionalContribution0; 
                      MaxFractionalContribution = MaxFractionalContribution0; 
                  end

                  % (4.6) prepare and perform the actual inversion of the data 
                  MinFractionalContribution_r = MinFractionalContribution(isnan(Xdirect));  % sample the original vetor of minimum fractional contributions
                  MaxFractionalContribution_r = MaxFractionalContribution(isnan(Xdirect));  % sample the original vetor of maximum fractional contributions

                  SolveCFList                 = ~ismember(ObsList,nCFList);                 % find the vector of positions to be included in the cost function evaluation
                  SolveCFList_r               = ~ismember(ObsList_r,nCFList);               % find the vector of positions to be included in the cost function evaluation
                  WeightingList_r             = WeightingList(ismember(ObsList,ObsList_r)); % find the reduced weighting vector
                  abspos_r                    = abspos(ismember(ObsList,ObsList_r));        % find the reduced vector of absolute positions
                  relpos_r                    = relpos(ismember(ObsList,ObsList_r));        % find the reduced vector of relative positions
                  sources                     = ismember(EMList0,EMsources);                % make a vector for which end-members count as sources
                  [X, EMIterInst0_updated, functioncost] = MEANDIR_InvertActiveSimulation(Solver, EMIterInst0, EMIterInst_r, RiverColumn0, RiverColumn_r, EMList0, Xdirect,...
                                                           MinFractionalContribution_r, MaxFractionalContribution_r, SolveCFList_r, abspos_r, relpos_r, WeightingList_r, fractionation, sources);

                  % (4.7) now evaluate the result of the inversion. if the inversion is acceptable, update the variables that store results
                  [successcounter, smpres, EMSucessIter0Matrix, EMSucessIter0Matrix_frac, RiverColumnSucessMatrix, smpresmisfit, smpfunmisfit, badmbcounter, InvRivDatsave] = ...
                  MEANDIR_EvaluateInversionInstance(i,X,nEM,nOL,EMIterInst0, EMIterInst0_updated,RiverColumn0,ObsList,ErrorCutMinMB,ErrorCutMaxMB,...
                  successcounter,InvRivDat,smpres,NormalizationType, badmbcounter, EMSucessIter0Matrix, EMSucessIter0Matrix_frac, RiverColumnSucessMatrix,...
                  smpresmisfit, smpfunmisfit, river, EM, EMList0, EMdatasource, InvRivDatsave, EMHaveSO4,...
                  MinFractionalContribution, MaxFractionalContribution, Solver, isotopeposition, carbonisotopematch,...
                  SolveCFList, abspos, relpos, WeightingList, ConvertDelta2RList, ConvertDelta2R, ImposeNormalizationCheck, EndMembersWithNegativeRatios, go, functioncost);
                              
                  if iterationcounter >= maxiterations; stopiterations = 1; end  % if the maximum number of attempts has been made, cease simulations
                  if successcounter   == maxsuccess;    stopiterations = 1; end  % if the desired number of successful simulations has been found, cease simulations
                  if isequal(IterateOver,'Samples') && successcounter==0 && iterationcounter>=maxzerohits; stopiterations = 1; fprintf('Hit the maxzerohits threshold, currently set to: %s\n',num2str(maxzerohits)); end

            end % END OF WHILE LOOP on 'stopiterations'

            % save the relevant simulation results 
            RiverColumnSucessMatrixcell{i,iterindx}  = RiverColumnSucessMatrix;  % save information on river chemistry
            EMSucessIter0Matrixcell{i,iterindx}      = EMSucessIter0Matrix;      % save information on end-member composition
            EMSucessIter0Matrixcell_frac{i,iterindx} = EMSucessIter0Matrix_frac; % save information on end-member composition
            iterationcountercell{i,iterindx}         = iterationcounter;         % save information on number of iterations attempted
            successcountercell{i,iterindx}           = successcounter;           % save information on number of successfull iterations
            badmbcountercell{i,iterindx}             = badmbcounter;             % save information on mass balance results
            smpresmisfitcell{i,iterindx}             = smpresmisfit;             % save information on model-observation misfit
            smpfunmisfitcell{i,iterindx}             = smpfunmisfit;             % save information on cost function misfit

            if  isequal(IterateOver,'Samples')                
                smprescell{i,iterindx}     = smpres;                             % save the results for fractional contributions to smprescell (smpres has a subfield for each variable, each with dimension nEM x maxsuccess)
                InvRivDatacell{i,iterindx} = InvRivDatsave;                      % save the results for actual inverted chemistry to InvRivDatacell (InvRivDatsave has a subfield for each variable, each with dimension maxsuccess x 1)
            end          
            
            if  isequal(IterateOver,'End-members')
                % save the results for fractional contributions to IterEM_smpfrac and for actual inverted chemistry to IterEM_InvRivData 
                IterEM_smpfrac    = MEANDIR_updateIterEM_smpfrac(smpres,IterEM_smpfrac,i,go);              % smpres has a subfield for each variable, each with dimension nEM x 1, IterEM_smpfrac has subfield for each variable, with dimension nEM x s
                IterEM_InvRivData = MEANDIR_updateIterEM_InvRivData(InvRivDatsave,IterEM_InvRivData,i,go); % InvRivDatsave has a subfield for each variable, each with dimension 1 x 1, IterEM_InvRivData has subfield for each variable, each with dimension s x 1
            end
            
        end % END OF INNER INVERSION LOOP (for/parfor i = 1:s)
        
        if  isequal(IterateOver,'End-members')
            % the IterEM_smpfracStore* have size of (nEM, numberiterations, s), IterEM_smpfrac subfields have dimension nEM x s
            % the IterEM_InvRivDataStore_* have size of (s, numberiteration), IterEM_InvRivData subfield have dimension s x 1
                               IterEM_smpfracStorenorm(:,iterindx,:)   = reshape(IterEM_smpfrac.smpresultsnorm,nEM,1,s);   IterEM_InvRivDataStore_norm(:,iterindx)   = IterEM_InvRivData.norm;       % norm
                               IterEM_smpfracStoreSO4(:,iterindx,:)    = reshape(IterEM_smpfrac.smpresultsSO4,nEM,1,s);    IterEM_InvRivDataStore_SO4(:,iterindx)    = IterEM_InvRivData.SO4;        % SO4
                               IterEM_smpfracStored34S(:,iterindx,:)   = reshape(IterEM_smpfrac.smpresultsd34S,nEM,1,s);   IterEM_InvRivDataStore_d34S(:,iterindx)   = IterEM_InvRivData.d34S;       % d34S
            if go.ALK    == 1; IterEM_smpfracStoreALK(:,iterindx,:)    = reshape(IterEM_smpfrac.smpresultsALK,nEM,1,s);    IterEM_InvRivDataStore_ALK(:,iterindx)    = IterEM_InvRivData.ALK;    end % ALK
            if go.DIC    == 1; IterEM_smpfracStoreDIC(:,iterindx,:)    = reshape(IterEM_smpfrac.smpresultsDIC,nEM,1,s);    IterEM_InvRivDataStore_DIC(:,iterindx)    = IterEM_InvRivData.DIC;    end % DIC
            if go.Ca     == 1; IterEM_smpfracStoreCa(:,iterindx,:)     = reshape(IterEM_smpfrac.smpresultsCa,nEM,1,s);     IterEM_InvRivDataStore_Ca(:,iterindx)     = IterEM_InvRivData.Ca;     end % Ca
            if go.Mg     == 1; IterEM_smpfracStoreMg(:,iterindx,:)     = reshape(IterEM_smpfrac.smpresultsMg,nEM,1,s);     IterEM_InvRivDataStore_Mg(:,iterindx)     = IterEM_InvRivData.Mg;     end % Mg
            if go.Na     == 1; IterEM_smpfracStoreNa(:,iterindx,:)     = reshape(IterEM_smpfrac.smpresultsNa,nEM,1,s);     IterEM_InvRivDataStore_Na(:,iterindx)     = IterEM_InvRivData.Na;     end % Na
            if go.K      == 1; IterEM_smpfracStoreK(:,iterindx,:)      = reshape(IterEM_smpfrac.smpresultsK,nEM,1,s);      IterEM_InvRivDataStore_K(:,iterindx)      = IterEM_InvRivData.K;      end % K
            if go.Sr     == 1; IterEM_smpfracStoreSr(:,iterindx,:)     = reshape(IterEM_smpfrac.smpresultsSr,nEM,1,s);     IterEM_InvRivDataStore_Sr(:,iterindx)     = IterEM_InvRivData.Sr;     end % Sr
            if go.Fe     == 1; IterEM_smpfracStoreFe(:,iterindx,:)     = reshape(IterEM_smpfrac.smpresultsFe,nEM,1,s);     IterEM_InvRivDataStore_Fe(:,iterindx)     = IterEM_InvRivData.Fe;     end % Fe
            if go.Cl     == 1; IterEM_smpfracStoreCl(:,iterindx,:)     = reshape(IterEM_smpfrac.smpresultsCl,nEM,1,s);     IterEM_InvRivDataStore_Cl(:,iterindx)     = IterEM_InvRivData.Cl;     end % Cl
            if go.NO3    == 1; IterEM_smpfracStoreNO3(:,iterindx,:)    = reshape(IterEM_smpfrac.smpresultsNO3,nEM,1,s);    IterEM_InvRivDataStore_NO3(:,iterindx)    = IterEM_InvRivData.NO3;    end % NO3
            if go.PO4    == 1; IterEM_smpfracStorePO4(:,iterindx,:)    = reshape(IterEM_smpfrac.smpresultsPO4,nEM,1,s);    IterEM_InvRivDataStore_PO4(:,iterindx)    = IterEM_InvRivData.PO4;    end % PO4
            if go.Si     == 1; IterEM_smpfracStoreSi(:,iterindx,:)     = reshape(IterEM_smpfrac.smpresultsSi,nEM,1,s);     IterEM_InvRivDataStore_Si(:,iterindx)     = IterEM_InvRivData.Si;     end % Si
            if go.Ge     == 1; IterEM_smpfracStoreGe(:,iterindx,:)     = reshape(IterEM_smpfrac.smpresultsGe,nEM,1,s);     IterEM_InvRivDataStore_Ge(:,iterindx)     = IterEM_InvRivData.Ge;     end % Ge
            if go.Li     == 1; IterEM_smpfracStoreLi(:,iterindx,:)     = reshape(IterEM_smpfrac.smpresultsLi,nEM,1,s);     IterEM_InvRivDataStore_Li(:,iterindx)     = IterEM_InvRivData.Li;     end % Li
            if go.F      == 1; IterEM_smpfracStoreF(:,iterindx,:)      = reshape(IterEM_smpfrac.smpresultsF,nEM,1,s);      IterEM_InvRivDataStore_F(:,iterindx)      = IterEM_InvRivData.F;      end % F
            if go.B      == 1; IterEM_smpfracStoreB(:,iterindx,:)      = reshape(IterEM_smpfrac.smpresultsB,nEM,1,s);      IterEM_InvRivDataStore_B(:,iterindx)      = IterEM_InvRivData.B;      end % B
            if go.Re     == 1; IterEM_smpfracStoreRe(:,iterindx,:)     = reshape(IterEM_smpfrac.smpresultsRe,nEM,1,s);     IterEM_InvRivDataStore_Re(:,iterindx)     = IterEM_InvRivData.Re;     end % Re
            if go.Mo     == 1; IterEM_smpfracStoreMo(:,iterindx,:)     = reshape(IterEM_smpfrac.smpresultsMo,nEM,1,s);     IterEM_InvRivDataStore_Mo(:,iterindx)     = IterEM_InvRivData.Mo;     end % Mo
            if go.Os     == 1; IterEM_smpfracStoreOs(:,iterindx,:)     = reshape(IterEM_smpfrac.smpresultsOs,nEM,1,s);     IterEM_InvRivDataStore_Os(:,iterindx)     = IterEM_InvRivData.Os;     end % Os
            if go.HCO3   == 1; IterEM_smpfracStoreHCO3(:,iterindx,:)   = reshape(IterEM_smpfrac.smpresultsHCO3,nEM,1,s);   IterEM_InvRivDataStore_HCO3(:,iterindx)   = IterEM_InvRivData.HCO3;   end % HCO3
            if go.d7Li   == 1; IterEM_smpfracStored7Li(:,iterindx,:)   = reshape(IterEM_smpfrac.smpresultsd7Li,nEM,1,s);   IterEM_InvRivDataStore_d7Li(:,iterindx)   = IterEM_InvRivData.d7Li;   end % d7Li
            if go.d13C   == 1; IterEM_smpfracStored13C(:,iterindx,:)   = reshape(IterEM_smpfrac.smpresultsd13C,nEM,1,s);   IterEM_InvRivDataStore_d13C(:,iterindx)   = IterEM_InvRivData.d13C;   end % d13C
            if go.d18O   == 1; IterEM_smpfracStored18O(:,iterindx,:)   = reshape(IterEM_smpfrac.smpresultsd18O,nEM,1,s);   IterEM_InvRivDataStore_d18O(:,iterindx)   = IterEM_InvRivData.d18O;   end % d18O
            if go.d26Mg  == 1; IterEM_smpfracStored26Mg(:,iterindx,:)  = reshape(IterEM_smpfrac.smpresultsd26Mg,nEM,1,s);  IterEM_InvRivDataStore_d26Mg(:,iterindx)  = IterEM_InvRivData.d26Mg;  end % d26Mg
            if go.d30Si  == 1; IterEM_smpfracStored30Si(:,iterindx,:)  = reshape(IterEM_smpfrac.smpresultsd30Si,nEM,1,s);  IterEM_InvRivDataStore_d30Si(:,iterindx)  = IterEM_InvRivData.d30Si;  end % d30Si
            if go.d42Ca  == 1; IterEM_smpfracStored42Ca(:,iterindx,:)  = reshape(IterEM_smpfrac.smpresultsd42Ca,nEM,1,s);  IterEM_InvRivDataStore_d42Ca(:,iterindx)  = IterEM_InvRivData.d42Ca;  end % d42Ca
            if go.d44Ca  == 1; IterEM_smpfracStored44Ca(:,iterindx,:)  = reshape(IterEM_smpfrac.smpresultsd44Ca,nEM,1,s);  IterEM_InvRivDataStore_d44Ca(:,iterindx)  = IterEM_InvRivData.d44Ca;  end % d44Ca
            if go.d56Fe  == 1; IterEM_smpfracStored56Fe(:,iterindx,:)  = reshape(IterEM_smpfrac.smpresultsd56Fe,nEM,1,s);  IterEM_InvRivDataStore_d56Fe(:,iterindx)  = IterEM_InvRivData.d56Fe;  end % d56Fe
            if go.Sr8786 == 1; IterEM_smpfracStoreSr8786(:,iterindx,:) = reshape(IterEM_smpfrac.smpresultsSr8786,nEM,1,s); IterEM_InvRivDataStore_Sr8786(:,iterindx) = IterEM_InvRivData.Sr8786; end % Sr8786
            if go.d98Mo  == 1; IterEM_smpfracStored98Mo(:,iterindx,:)  = reshape(IterEM_smpfrac.smpresultsd98Mo,nEM,1,s);  IterEM_InvRivDataStore_d98Mo(:,iterindx)  = IterEM_InvRivData.d98Mo;  end % d98Mo
            if go.Os8788 == 1; IterEM_smpfracStoreOs8788(:,iterindx,:) = reshape(IterEM_smpfrac.smpresultsOs8788,nEM,1,s); IterEM_InvRivDataStore_Os8788(:,iterindx) = IterEM_InvRivData.Os8788; end % Os8788
            if go.Fmod   == 1; IterEM_smpfracStoreFmod(:,iterindx,:)   = reshape(IterEM_smpfrac.smpresultsFmod,nEM,1,s);   IterEM_InvRivDataStore_Fmod(:,iterindx)   = IterEM_InvRivData.Fmod;   end % Fmod
        end
        
    end % END OF OUTER INVERSION LOOP (for/parfor iterindx = 1:numberiterations)

    disp('aggregating the results');
    if isequal(IterateOver,'End-members')
       % take the output of the scenario and un-pack the results for MEANDIR_UnpackInversionResults
       MEANDIR_PrepareIterEMOutputForUnpacking
    end          
    
    % unpack the results of the inversion model
    MEANDIR_UnpackInversionResults 
        
    %% SECTION 5: CALCULATIONS
    disp('performing calculations on results');   
    MEANDIR_CalculateCullDatabyMisfit          % when iterating over end-members, cull data by misfit
    MEANDIR_CalculateMassBalance               % calculate and save the fraction of reconstructed river observations
    MEANDIR_SaveFractionalContributions        % save fractional contributions of end-members to river observations
    MEANDIR_SaveInvertedRiverValues            % save the InvRivDat date into the data structure
    MEANDIR_SaveFractionations                 % save the fractionation factors used for each sample and inversion
    MEANDIR_CalculateReconstructedObservations % calculate and save inversion-constrained river chemistry
    MEANDIR_CalculateExcessSO4                 % calculate and save the excess SO4 from the simulation
    MEANDIR_CalculateOther_d34S                % calculate and save expected value for the pyrite/SO4 excess d34S
    MEANDIR_CalculateEndMemberValues           % calculate and save the composition of the end-members
    MEANDIR_CalculateRZCWY                     % calculate and save R,Z,C,Y, W for each sample
    clear results;
    
    %% SECTION 6: REORDER AND SAVE RESULTS 
    MEANDIR_SaveInversionSettings    
    % save the matrix of inversion end-members
    river.inversionmatrix = EMSucessIterMatrixAll; 
    clear EMSucessIterMatrixAll; 
    clear EMSucessIterMatrixAll_frac; 

    % re-order the structures of river 
    clear reorder_river;
    reorder_river = struct;
    fn = fieldnames(river);
    if sum(ismember(fn,'settings'))               >0; reorder_river.settings               = 1;  end  %  settings
    if sum(ismember(fn,'info'))                   >0; reorder_river.info                   = 1;  end  %  info
    if sum(ismember(fn,'observations'))           >0; reorder_river.observations           = 1;  end  %  observations
    if sum(ismember(fn,'model_variable'))         >0; reorder_river.model_variable         = 1;  end  %  model_variable
    if sum(ismember(fn,'RiverMatrix'))            >0; reorder_river.RiverMatrix            = 1;  end  %  RiverMatrix
    if sum(ismember(fn,'InvertedRiverValues'))    >0; reorder_river.InvertedRiverValues    = 1;  end  %  InvertedRiverValues
    if sum(ismember(fn,'fractionation'))          >0; reorder_river.fractionation          = 1;  end  %  fractionation
    if sum(ismember(fn,'inversionmatrix'))        >0; reorder_river.inversionmatrix        = 1;  end  %  inversionmatrix
    if sum(ismember(fn,'reconstructed'))          >0; reorder_river.reconstructed          = 1;  end  %  reconstructed
    if sum(ismember(fn,'fraction'))               >0; reorder_river.fraction               = 1;  end  %  fraction
    if sum(ismember(fn,'misfit'))                 >0; reorder_river.misfit                 = 1;  end  %  misfit
    if sum(ismember(fn,'massbalance'))            >0; reorder_river.massbalance            = 1;  end  %  massbalance
    if sum(ismember(fn,'EndMembers'))             >0; reorder_river.EndMembers             = 1;  end  %  EndMembers
    if sum(ismember(fn,'RZCWY'))                  >0; reorder_river.RZCWY                  = 1;  end  %  RZCWY
    if sum(ismember(fn,'calculated_other_d34S'))  >0; reorder_river.calculated_other_d34S  = 1;  end  %  calculated_other_d34S
    if sum(ismember(fn,'ExcessSO4'))              >0; reorder_river.ExcessSO4              = 1;  end  %  ExcessSO4
    river = orderfields(river,reorder_river);    
    clear reorder_river;
    
    % determine the filename for saving inversion results
    d = date; d = d(~ismember(d,'-')); % find the date today (saved in the filename of output)
    if     isequal(river.settings.IterateOver,'Samples');     savename = ['MEANDIR_' d '_' ScenarioList{ScenarioListIndex} '_max_' num2str(river.settings.maxsuccess) '_try_' num2str(floor(river.settings.maxiterations/(10^floor(log10(river.settings.maxiterations))))) 'e' num2str(floor(log10(river.settings.maxiterations)))];
    elseif isequal(river.settings.IterateOver,'End-members'); savename = ['MEANDIR_' d '_' ScenarioList{ScenarioListIndex} '_n' num2str(river.settings.numberiterations) '_cut' num2str(river.settings.MisfitCuts)]; end

    % end the timer and record elapsed duration
    t = toc; river.info.time_elapsed = t;
    river.info.functionalsamplelist  = functionalsamplelist; 
    
    % define the variables to save and then save them. display the name of the file, the 
    % location where it was saved, and inform the user that the inversion is complete
    % cd /Users/prestonkemeny/Documents/Alaska/MEANDIR_Alaska
    save(savename,'river','-v7.3');
    savenames{ScenarioListIndex,1} = savename;
    fprintf('saved the data file: %s\n',savename)
    fprintf('saved the data file in folder: %s\n',cd)
    disp(sprintf('completed simulation "%s".',ScenarioList{ScenarioListIndex}));


    
end % END OF LOOP OVER ScenarioList
disp('all requested scenarios have been completed. thank you for using MEANDIR!'); % END OF MEANDIR