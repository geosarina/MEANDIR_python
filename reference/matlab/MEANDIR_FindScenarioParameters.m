% this script identifies the parameters associated with each entry in "ScenarioList".
% the function of each variable is briefly summarized and described in more detail in the user guide.

% set FailCase equal to 0
 FailCase = 0;
 % FailCase is a variable that will be set to 1 when a requested simulation does not meet requirements. for example, if a variable should take on integer values but the 
 % user input is 'caterpillar', FailCase will be set equal to 1 rather than its default value of 0. although a series of common possible errors are
 % tested, FailCase will not catch all possible mistakes. if MEANDIR identifies a known problem, it will inform whe user why it has halted the scenario and 
 % how to change the variable values to allow the inversion to continue. MEANDIR_Master loops through all requested inversions in ScenarioList to
 % check that FailCase does not attain a value of 1 prior to beginning the scenarios.  


%% (1) scenario parameters
if isequal(ScenarioList{ScenarioListIndex},'AK_scenario1_carbonate_slctindi')     
  % group 1/6, river observations - 12 variables
    Riverdatasource              = 'PCKAlaskaData';                                           % name for the source of the river observations
    AdjustRiverObs               = 1;                                                         % whether river measurements should be adjusted to reflect analytical uncertainty. options are: 0 or 1   
    ObsList                      = {'Ca', 'Mg', 'Na', 'K',  'Cl', 'SO4','DIC','d34S','d13C'}; % list of elements, ions, and isotopic ratios to include in the inversion
    CostFunType                  = {'rel','rel','rel','rel','rel','rel','rel','rel', 'rel'};  % whether each variable should be evaluated as proportional cost or absolute cost. options are 'rel' and 'abs'.
    WeightingList                = [1      1     1     1     1     1     1      1      1];    % weighting term for the cost function                     
    ErrorCutMinMB                = [95     95    95    95    95    95    95    -1     -1];    % if IterateOver is "Samples", these are the bounds of mass balance that will be accepted for an inversion result  
    ErrorCutMaxMB                = [105    105   105   105   105   105   105   +1     +1];    % if IterateOver is "Samples", these are the bounds of mass balance that will be accepted for an inversion result             
    nCFList                      = {};                                                        % list of dissolved variables (in ObsList) to not include in the cost function evaluation
    ConvertDelta2RList           = {'d13C','d34S'};                                           % list of isotopic variables (in delta notation) to be converted to isotopic ratios for the analysis (results are converted back to delta)
    AllIonsExplicitlyResolved    = 0;                                                         % Indicate if all cations and anions (all sources of charge) are resolved. options are: 0 or 1    
    ObsInNormalization           = {'Ca','Mg','Na','K','SO4','DIC'};                          % list of variables included in the normalization
    ImposeNormalizationCheck     = 1;                                                         % whether the normalization should be reproduced to the minimum and maximum values for non-isotopic variables 
  % group 2/6, general settings - 8 variables
    Solver                       = 'mldivide_optimize';                                       % type of inversion solution, options are: 'mldivide', 'lsqnonneg', 'mldivide_optimize', 'lsqnonneg_optimize', 'optimize'. enter as string, not a cell.
    IterateOver                  = 'Samples';                                                 % whether to apply the same end-members to all sample ('End-members'), or different end-members to each sample ('Sample'). options are: 'End-members', 'Samples'
    maxiterations                = 5e5;                                                       % if IterateOver equals "Samples", this is the maximum number of inversion attempts     
    maxsuccess                   = 200;                                                       % if IterateOver is "Samples", this is the desired number of successes      
    maxzerohits                  = 1e4;                                                       % if IterateOver is "Samples", this is the maximum number of iterations to try without success before moving on
    numberiterations             = NaN;                                                       % if IterateOver equals "End-members", this is the number of inversions to perform
    MisfitCuts                   = NaN;                                                       % if IterateOver equals "End-members", this percentage of simulations with lowest misfit between model result and data will be kept    
    CullOn                       = 'EachSample';                                              % if IterateOver is "End-members", this must be either 'EachSample' or 'AllSample' and describes HOW to cull the data
    saveuncutdata                = 0;                                                         % if 1, saves all the simulation data (results in very large files)    
  % group 3/6, end-members - 14 variables
    EMdatasource                 = 'PCK23_Alaska_SumCatSO4DIC_uniform';                       % the name of the end-member group (the entry on column A of the spreadsheet containing end-member information). 
    EMList0                      = {'prec','carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri','evap','corg'};  % list of the end-members to use in the inversion (row 2 of spreadsheet containing end-member information).     
    MinFractionalContribution    = [0       0       0          0         0        0        0       0      0];    % the minumum fractional contribution of each end-member to the normalization ion (typically 0, unless secondary mineral formation is in inversion)
    MaxFractionalContribution    = [1       1       1          1         1        1        1       1      1];    % the maximum fractional contribution of each end-member to the normalization ion (typically 1, unless secondary mineral formation is in inversion) 
    ListNormClosure              = {'Ca','Ca','Ca','Mg','Na','K','SO4','Ca','DIC'};            % the variable for each end-member that will be calculated by mass balance to ensure interal consistency
    ListChargeClosure            = {};                                                         % the variable for each end-member that will be calculated by charge balance to ensure internal consistency 
    EMUnits                      = 'equi';                                                     % whether the end-member chemistry is given as concentration or charge-equivalent. options are: 'conc', 'equi'    
    EMsources                    = {'prec','carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri','evap','corg'};  % end-members that count as sources of dissolved constituents (typically rocks and other solute sources)
    EMsinks                      = {};                                                         % end-members that count as sinks of dissolved constituents (clays and other secondary phases)         
    EndMembersWithNegativeRatios = {};                                                         % end-members with negative chemical ratios    
    CoupleFeS2SO4intoEM          = {};                                                         % list of end-members where SO4 ratios represent pyrite
    CoupleFeS2d34SintoEM         = {};                                                         % list of end-members where SO4 d34S values represent pyrite
    RecordFullFeS2Distribution   = 0;                                                          % whether to record the full distribution of calculated FeS2 values, when d34S is not included in the inversion
    BalanceEvaporite             = 1;                                                          % whether evaporite SO4 = Ca+Mg+Sr and Cl = Na+K.    
  % group 4/6, Cl correction - 2 variables
    PrecProcessing               = 'EndMember';                                                % whether Cl should be treated as a common ion ('EndMember') or if there are Cl critical values ('ClCrit'). options are: 'EndMember', 'ClCrit'                            
    ClCriticalValuesGiven        = 1;                                                          % whether Cl critical values are given. if PrecProcessing is 'ClCrit' and ClCriticalValuesGiven is 1, those ClCritical values are used. if PrecProcessing is 'ClCrit' and ClCriticalValuesGiven is 0, 100% of river Cl is used.        
  % group 5/6, calculation of RZWYC - 8 variables
    CalculateRZCWY               = 1;                                                          % whether model should attempt to calculate R, Z, C, W, and Y. options are: 0 or 1
    R_Numerator_EMList           = {'carb'};                                                   % end-members that should contribute to the numerator of R
    R_Numerator_IonList          = {'Na','Ca','Mg','K'};                                       % ions that should contribute to the numerator of R        
    Z_NumeratorType              = {'ZfromEM'};                                                % how the numerator of Z (SO4 from FeS2) should be calculated. options are: 'ZfromEM', 'ZfromSO4excess', 'ZfromriverSO4', 'Znotcalculated'
    Z_Numerator_EMList           = {'pyri'};                                                   % end-members that should contribute to the numerator of Z (only used if Z_NumeratorType = {'ZfromEM'}).       
    C_Numerator_EMList           = {'corg'};                                                   % End-member that should contribute to the numerator of C
    RZC_Denominator_EMList       = {'carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri'};     % end-members that should contribute to the denominator of R, Z, and C (weathering end-members)
    RZC_Denominator_IonList      = {'Na','Ca','Mg','K'};                                       % ions that should contribute to the denominator of R, Z, and C. we recommend including Na and K.                                 
    MinFractionalContribution0   = MinFractionalContribution;
    MaxFractionalContribution0   = MaxFractionalContribution; 

% scenario 2: AK_scenario2_carbonate_slctindi_degas_2p5
elseif isequal(ScenarioList{ScenarioListIndex},'AK_scenario2_carbonate_slctindi_degas_2p5')     
  % group 1/6, river observations - 12 variables
    Riverdatasource              = 'PCKAlaskaData';                                            % name for the source of the river observations
    AdjustRiverObs               = 1;                                                          % whether river measurements should be adjusted to reflect analytical uncertainty. options are: 0 or 1   
    ObsList                      = {'Ca', 'Mg', 'Na', 'K',  'Cl', 'SO4','DIC','d34S','d13C'};  % list of elements, ions, and isotopic ratios to include in the inversion
    CostFunType                  = {'rel','rel','rel','rel','rel','rel','rel','rel', 'rel'};   % whether each variable should be evaluated as proportional cost or absolute cost. options are 'rel' and 'abs'.
    WeightingList                = [1      1     1     1     1     1     1      1      1];     % weighting term for the cost function                     
    ErrorCutMinMB                = [95     95    95    95    95    95    95    -1     -1];     % if IterateOver is "Samples", these are the bounds of mass balance that will be accepted for an inversion result  
    ErrorCutMaxMB                = [105    105   105   105   105   105   105   +1     +1];     % if IterateOver is "Samples", these are the bounds of mass balance that will be accepted for an inversion result                
    nCFList                      = {};                                                         % list of dissolved variables (in ObsList) to not include in the cost function evaluation
    ConvertDelta2RList           = {'d13C','d34S'};                                            % list of isotopic variables (in delta notation) to be converted to isotopic ratios for the analysis (results are converted back to delta)
    AllIonsExplicitlyResolved    = 0;                                                          % Indicate if all cations and anions (all sources of charge) are resolved. options are: 0 or 1    
    ObsInNormalization           = {'Ca','Mg','Na','K','SO4','DIC'};                           % list of variables included in the normalization
    ImposeNormalizationCheck     = 1;                                                          % whether the normalization should be reproduced to the minimum and maximum values for non-isotopic variables 
  % group 2/6, general settings - 8 variables
    Solver                       = 'mldivide_optimize';                                        % type of inversion solution, options are: 'mldivide', 'lsqnonneg', 'mldivide_optimize', 'lsqnonneg_optimize', 'optimize'. enter as string, not a cell.
    IterateOver                  = 'Samples';                                                  % whether to apply the same end-members to all sample ('End-members'), or different end-members to each sample ('Sample'). options are: 'End-members', 'Samples'
    maxiterations                = 5e5;                                                       % if IterateOver equals "Samples", this is the maximum number of inversion attempts     
    maxsuccess                   = 200;                                                       % if IterateOver is "Samples", this is the desired number of successes      
    maxzerohits                  = 1e4;                                                       % if IterateOver is "Samples", this is the maximum number of iterations to try without success before moving on
    numberiterations             = NaN;                                                       % if IterateOver equals "End-members", this is the number of inversions to perform
    MisfitCuts                   = NaN;                                                       % if IterateOver equals "End-members", this percentage of simulations with lowest misfit between model result and data will be kept    
    CullOn                       = 'EachSample';                                              % if IterateOver is "End-members", this must be either 'EachSample' or 'AllSample' and describes HOW to cull the data
    saveuncutdata                = 0;                                                         % if 1, saves all the simulation data (results in very large files)    
  % group 3/6, end-members - 14 variables
    EMdatasource                 = 'PCK23_Alaska_SumCatSO4DIC_uniform';                        % the name of the end-member group (the entry on column A of the spreadsheet containing end-member information). 
    EMList0                      = {'prec','carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri','evap','corg','degas'}; % list of the end-members to use in the inversion (row 2 of spreadsheet containing end-member information).     
    MinFractionalContribution    = [0       0         0        0         0        0        0       0     0      -9.99]; % the minumum fractional contribution of each end-member to the normalization ion (typically 0, unless secondary mineral formation is in inversion)
    MaxFractionalContribution    = [1       inf       1        1         1        1        1       1     inf    -9.99]; % the maximum fractional contribution of each end-member to the normalization ion (typically 1, unless secondary mineral formation is in inversion) 
    ListNormClosure              = {'Ca','Ca','Ca','Mg','Na','K','SO4','Ca','DIC','DIC'};      % the variable for each end-member that will be calculated by mass balance to ensure interal consistency
    ListChargeClosure            = {};                                                         % the variable for each end-member that will be calculated by charge balance to ensure internal consistency 
    EMUnits                      = 'equi';                                                     % whether the end-member chemistry is given as concentration or charge-equivalent. options are: 'conc', 'equi'    
    EMsources                    = {'prec','carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri','evap','corg'};  % end-members that count as sources of dissolved constituents (typically rocks and other solute sources)
    EMsinks                      = {};                                                         % end-members that count as sinks of dissolved constituents (clays and other secondary phases)
    EndMembersWithNegativeRatios = {};                                                         % end-members with negative chemical ratios    
    CoupleFeS2SO4intoEM          = {};                                                         % list of end-members where SO4 ratios represent pyrite
    CoupleFeS2d34SintoEM         = {};                                                         % list of end-members where SO4 d34S values represent pyrite
    RecordFullFeS2Distribution   = 0;                                                          % whether to record the full distribution of calculated FeS2 values, when d34S is not included in the inversion
    BalanceEvaporite             = 1;                                                          % whether evaporite SO4 = Ca+Mg+Sr and Cl = Na+K.    
  % group 4/6, Cl correction - 2 variables
    PrecProcessing               = 'EndMember';                                                % whether Cl should be treated as a common ion ('EndMember') or if there are Cl critical values ('ClCrit'). options are: 'EndMember', 'ClCrit'                            
    ClCriticalValuesGiven        = 1;                                                          % whether Cl critical values are given. if PrecProcessing is 'ClCrit' and ClCriticalValuesGiven is 1, those ClCritical values are used. if PrecProcessing is 'ClCrit' and ClCriticalValuesGiven is 0, 100% of river Cl is used.        
  % group 5/6, calculation of RZWYC - 8 variables
    CalculateRZCWY               = 1;                                                          % whether model should attempt to calculate R, Z, C, W, and Y. options are: 0 or 1
    R_Numerator_EMList           = {'carb'};                                                   % end-members that should contribute to the numerator of R
    R_Numerator_IonList          = {'Na','Ca','Mg','K'};                                       % ions that should contribute to the numerator of R        
    Z_NumeratorType              = {'ZfromEM'};                                                % how the numerator of Z (SO4 from FeS2) should be calculated. options are: 'ZfromEM', 'ZfromSO4excess', 'ZfromriverSO4', 'Znotcalculated'
    Z_Numerator_EMList           = {'pyri'};                                                   % end-members that should contribute to the numerator of Z (only used if Z_NumeratorType = {'ZfromEM'}).       
    C_Numerator_EMList           = {'corg'};                                                   % End-member that should contribute to the numerator of C
    RZC_Denominator_EMList       = {'carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri'};     % end-members that should contribute to the denominator of R, Z, and C (weathering end-members)
    RZC_Denominator_IonList      = {'Na','Ca','Mg','K'};                                       % ions that should contribute to the denominator of R, Z, and C. we recommend including Na and K.                              
  % group 6/6, setting specific contributions from organic carbon
    ResetDegasDICContribution    = 1;
    DegasDICContributionMin      = -2.5;  % fraction relative to the DIC in solution, not the normalization
    DegasDICContributionMax      = +0.00; % fraction relative to the DIC in solution, not the normalization
    MinFractionalContribution0   = MinFractionalContribution;
    MaxFractionalContribution0   = MaxFractionalContribution; 
    
% scenario 3: AK_scenario3_carbonate_slctindi_degas_25
elseif isequal(ScenarioList{ScenarioListIndex},'AK_scenario3_carbonate_slctindi_degas_25')     
  % group 1/6, river observations - 12 variables
    Riverdatasource              = 'PCKAlaskaData';                                            % name for the source of the river observations
    AdjustRiverObs               = 1;                                                          % whether river measurements should be adjusted to reflect analytical uncertainty. options are: 0 or 1   
    ObsList                      = {'Ca', 'Mg', 'Na', 'K',  'Cl', 'SO4','DIC','d34S','d13C'};  % list of elements, ions, and isotopic ratios to include in the inversion
    CostFunType                  = {'rel','rel','rel','rel','rel','rel','rel','rel', 'rel'};   % whether each variable should be evaluated as proportional cost or absolute cost. options are 'rel' and 'abs'.
    WeightingList                = [1      1     1     1     1     1     1      1      1];     % weighting term for the cost function                     
    ErrorCutMinMB                = [95     95    95    95    95    95    95    -1     -1];     % if IterateOver is "Samples", these are the bounds of mass balance that will be accepted for an inversion result  
    ErrorCutMaxMB                = [105    105   105   105   105   105   105   +1     +1];     % if IterateOver is "Samples", these are the bounds of mass balance that will be accepted for an inversion result             
    nCFList                      = {};                                                         % list of dissolved variables (in ObsList) to not include in the cost function evaluation
    ConvertDelta2RList           = {'d13C','d34S'};                                            % list of isotopic variables (in delta notation) to be converted to isotopic ratios for the analysis (results are converted back to delta)
    AllIonsExplicitlyResolved    = 0;                                                          % Indicate if all cations and anions (all sources of charge) are resolved. options are: 0 or 1    
    ObsInNormalization           = {'Ca','Mg','Na','K','SO4','DIC'};                           % list of variables included in the normalization
    ImposeNormalizationCheck     = 1;                                                          % whether the normalization should be reproduced to the minimum and maximum values for non-isotopic variables 
  % group 2/6, general settings - 8 variables
    Solver                       = 'mldivide_optimize';                                        % type of inversion solution, options are: 'mldivide', 'lsqnonneg', 'mldivide_optimize', 'lsqnonneg_optimize', 'optimize'. enter as string, not a cell.
    IterateOver                  = 'Samples';                                                  % whether to apply the same end-members to all sample ('End-members'), or different end-members to each sample ('Sample'). options are: 'End-members', 'Samples'
    maxiterations                = 5e5;                                                        % if IterateOver equals "Samples", this is the maximum number of inversion attempts     
    maxsuccess                   = 200;                                                        % if IterateOver is "Samples", this is the desired number of successes      
    maxzerohits                  = 1e4;                                                       % if IterateOver is "Samples", this is the maximum number of iterations to try without success before moving on
    numberiterations             = NaN;                                                       % if IterateOver equals "End-members", this is the number of inversions to perform
    MisfitCuts                   = NaN;                                                       % if IterateOver equals "End-members", this percentage of simulations with lowest misfit between model result and data will be kept    
    CullOn                       = 'EachSample';                                               % if IterateOver is "End-members", this must be either 'EachSample' or 'AllSample' and describes HOW to cull the data
    saveuncutdata                = 0;                                                          % if 1, saves all the simulation data (results in very large files)    
  % group 3/6, end-members - 14 variables
    EMdatasource                 = 'PCK23_Alaska_SumCatSO4DIC_uniform';                        % the name of the end-member group (the entry on column A of the spreadsheet containing end-member information). 
    EMList0                      = {'prec','carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri','evap','corg','degas'}; % list of the end-members to use in the inversion (row 2 of spreadsheet containing end-member information).     
    MinFractionalContribution    = [0       0        0          0         0       0        0       0     0      -9.99]; % the minumum fractional contribution of each end-member to the normalization ion (typically 0, unless secondary mineral formation is in inversion)
    MaxFractionalContribution    = [1       inf      1          1         1       1        1       1     inf    -9.99]; % the maximum fractional contribution of each end-member to the normalization ion (typically 1, unless secondary mineral formation is in inversion) 
    ListNormClosure              = {'Ca','Ca','Ca','Mg','Na','K','SO4','Ca','DIC','DIC'};      % the variable for each end-member that will be calculated by mass balance to ensure interal consistency
    ListChargeClosure            = {};                                                         % the variable for each end-member that will be calculated by charge balance to ensure internal consistency 
    EMUnits                      = 'equi';                                                     % whether the end-member chemistry is given as concentration or charge-equivalent. options are: 'conc', 'equi'    
    EMsources                    = {'prec','carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri','evap','corg'};  % end-members that count as sources of dissolved constituents (typically rocks and other solute sources)
    EMsinks                      = {};                                                         % end-members that count as sinks of dissolved constituents (clays and other secondary phases)
    EndMembersWithNegativeRatios = {};                                                         % end-members with negative chemical ratios    
    CoupleFeS2SO4intoEM          = {};                                                         % list of end-members where SO4 ratios represent pyrite
    CoupleFeS2d34SintoEM         = {};                                                         % list of end-members where SO4 d34S values represent pyrite
    RecordFullFeS2Distribution   = 0;                                                          % whether to record the full distribution of calculated FeS2 values, when d34S is not included in the inversion
    BalanceEvaporite             = 1;                                                          % whether evaporite SO4 = Ca+Mg+Sr and Cl = Na+K.    
  % group 4/6, Cl correction - 2 variables
    PrecProcessing               = 'EndMember';                                                % whether Cl should be treated as a common ion ('EndMember') or if there are Cl critical values ('ClCrit'). options are: 'EndMember', 'ClCrit'                            
    ClCriticalValuesGiven        = 1;                                                          % whether Cl critical values are given. if PrecProcessing is 'ClCrit' and ClCriticalValuesGiven is 1, those ClCritical values are used. if PrecProcessing is 'ClCrit' and ClCriticalValuesGiven is 0, 100% of river Cl is used.        
  % group 5/6, calculation of RZWYC - 8 variables
    CalculateRZCWY               = 1;                                                          % whether model should attempt to calculate R, Z, C, W, and Y. options are: 0 or 1
    R_Numerator_EMList           = {'carb'};                                                   % end-members that should contribute to the numerator of R
    R_Numerator_IonList          = {'Na','Ca','Mg','K'};                                       % ions that should contribute to the numerator of R        
    Z_NumeratorType              = {'ZfromEM'};                                                % how the numerator of Z (SO4 from FeS2) should be calculated. options are: 'ZfromEM', 'ZfromSO4excess', 'ZfromriverSO4', 'Znotcalculated'
    Z_Numerator_EMList           = {'pyri'};                                                   % end-members that should contribute to the numerator of Z (only used if Z_NumeratorType = {'ZfromEM'}).       
    C_Numerator_EMList           = {'corg'};                                                   % End-member that should contribute to the numerator of C
    RZC_Denominator_EMList       = {'carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri'};     % end-members that should contribute to the denominator of R, Z, and C (weathering end-members)
    RZC_Denominator_IonList      = {'Na','Ca','Mg','K'};                                       % ions that should contribute to the denominator of R, Z, and C. we recommend including Na and K.                              
  % group 6/6, setting specific contributions from organic carbon
    ResetDegasDICContribution    = 1;
    DegasDICContributionMin      = -25;   % fraction relative to the DIC in solution, not the normalization
    DegasDICContributionMax      = +0.00; % fraction relative to the DIC in solution, not the normalization
    MinFractionalContribution0   = MinFractionalContribution;
    MaxFractionalContribution0   = MaxFractionalContribution; 



% scenario 4: AK_scenario4_carbonate_slctindi_degas_2p5_highfrac
elseif isequal(ScenarioList{ScenarioListIndex},'AK_scenario4_carbonate_slctindi_degas_2p5_highfrac')     
  % group 1/6, river observations - 12 variables
    Riverdatasource              = 'PCKAlaskaData';                                           % name for the source of the river observations
    AdjustRiverObs               = 1;                                                         % whether river measurements should be adjusted to reflect analytical uncertainty. options are: 0 or 1   
    ObsList                      = {'Ca', 'Mg', 'Na', 'K',  'Cl', 'SO4','DIC','d34S','d13C'}; % list of elements, ions, and isotopic ratios to include in the inversion
    CostFunType                  = {'rel','rel','rel','rel','rel','rel','rel','rel', 'rel'};  % whether each variable should be evaluated as proportional cost or absolute cost. options are 'rel' and 'abs'.
    WeightingList                = [1      1     1     1     1     1     1      1      1];    % weighting term for the cost function                     
    ErrorCutMinMB                = [95     95    95    95    95    95    95    -1     -1];    % if IterateOver is "Samples", these are the bounds of mass balance that will be accepted for an inversion result  
    ErrorCutMaxMB                = [105    105   105   105   105   105   105   +1     +1];    % if IterateOver is "Samples", these are the bounds of mass balance that will be accepted for an inversion result             
    nCFList                      = {};                                                        % list of dissolved variables (in ObsList) to not include in the cost function evaluation
    ConvertDelta2RList           = {'d13C','d34S'};                                           % list of isotopic variables (in delta notation) to be converted to isotopic ratios for the analysis (results are converted back to delta)
    AllIonsExplicitlyResolved    = 0;                                                         % Indicate if all cations and anions (all sources of charge) are resolved. options are: 0 or 1    
    ObsInNormalization           = {'Ca','Mg','Na','K','SO4','DIC'};                          % list of variables included in the normalization
    ImposeNormalizationCheck     = 1;                                                         % whether the normalization should be reproduced to the minimum and maximum values for non-isotopic variables 
  % group 2/6, general settings - 8 variables
    Solver                       = 'mldivide_optimize';                                       % type of inversion solution, options are: 'mldivide', 'lsqnonneg', 'mldivide_optimize', 'lsqnonneg_optimize', 'optimize'. enter as string, not a cell.
    IterateOver                  = 'Samples';                                                 % whether to apply the same end-members to all sample ('End-members'), or different end-members to each sample ('Sample'). options are: 'End-members', 'Samples'
    maxiterations                = 5e5;                                                       % if IterateOver equals "Samples", this is the maximum number of inversion attempts     
    maxsuccess                   = 200;                                                       % if IterateOver is "Samples", this is the desired number of successes      
    maxzerohits                  = 1e4;                                                       % if IterateOver is "Samples", this is the maximum number of iterations to try without success before moving on
    numberiterations             = NaN;                                                       % if IterateOver equals "End-members", this is the number of inversions to perform
    MisfitCuts                   = NaN;                                                       % if IterateOver equals "End-members", this percentage of simulations with lowest misfit between model result and data will be kept    
    CullOn                       = 'EachSample';                                              % if IterateOver is "End-members", this must be either 'EachSample' or 'AllSample' and describes HOW to cull the data
    saveuncutdata                = 0;                                                         % if 1, saves all the simulation data (results in very large files)    
  % group 3/6, end-members - 14 variables
    EMdatasource                 = 'PCK23_Alaska_SumCatSO4DIC_uniform_highfrac';              % the name of the end-member group (the entry on column A of the spreadsheet containing end-member information). 
    EMList0                      = {'prec','carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri','evap','corg','degas'};  % list of the end-members to use in the inversion (row 2 of spreadsheet containing end-member information).     
    MinFractionalContribution    = [0       0       0          0         0        0        0       0      0      -9.99]; % the minumum fractional contribution of each end-member to the normalization ion (typically 0, unless secondary mineral formation is in inversion)
    MaxFractionalContribution    = [1       inf     1          1         1        1        1       1      inf    -9.99]; % the maximum fractional contribution of each end-member to the normalization ion (typically 1, unless secondary mineral formation is in inversion) 
    ListNormClosure              = {'Ca','Ca','Ca','Mg','Na','K','SO4','Ca','DIC','DIC'};     % the variable for each end-member that will be calculated by mass balance to ensure interal consistency
    ListChargeClosure            = {};                                                        % the variable for each end-member that will be calculated by charge balance to ensure internal consistency 
    EMUnits                      = 'equi';                                                    % whether the end-member chemistry is given as concentration or charge-equivalent. options are: 'conc', 'equi'    
    EMsources                    = {'prec','carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri','evap','corg'};  % end-members that count as sources of dissolved constituents (typically rocks and other solute sources)
    EMsinks                      = {};                                                        % end-members that count as sinks of dissolved constituents (clays and other secondary phases)
    EndMembersWithNegativeRatios = {};                                                        % end-members with negative chemical ratios    
    CoupleFeS2SO4intoEM          = {};                                                        % list of end-members where SO4 ratios represent pyrite
    CoupleFeS2d34SintoEM         = {};                                                        % list of end-members where SO4 d34S values represent pyrite
    RecordFullFeS2Distribution   = 0;                                                         % whether to record the full distribution of calculated FeS2 values, when d34S is not included in the inversion
    BalanceEvaporite             = 1;                                                         % whether evaporite SO4 = Ca+Mg+Sr and Cl = Na+K.    
  % group 4/6, Cl correction - 2 variables
    PrecProcessing               = 'EndMember';                                               % whether Cl should be treated as a common ion ('EndMember') or if there are Cl critical values ('ClCrit'). options are: 'EndMember', 'ClCrit'                            
    ClCriticalValuesGiven        = 1;                                                         % whether Cl critical values are given. if PrecProcessing is 'ClCrit' and ClCriticalValuesGiven is 1, those ClCritical values are used. if PrecProcessing is 'ClCrit' and ClCriticalValuesGiven is 0, 100% of river Cl is used.        
  % group 5/6, calculation of RZWYC - 8 variables
    CalculateRZCWY               = 1;                                                         % whether model should attempt to calculate R, Z, C, W, and Y. options are: 0 or 1
    R_Numerator_EMList           = {'carb'};                                                  % end-members that should contribute to the numerator of R
    R_Numerator_IonList          = {'Na','Ca','Mg','K'};                                      % ions that should contribute to the numerator of R        
    Z_NumeratorType              = {'ZfromEM'};                                               % how the numerator of Z (SO4 from FeS2) should be calculated. options are: 'ZfromEM', 'ZfromSO4excess', 'ZfromriverSO4', 'Znotcalculated'
    Z_Numerator_EMList           = {'pyri'};                                                  % end-members that should contribute to the numerator of Z (only used if Z_NumeratorType = {'ZfromEM'}).       
    C_Numerator_EMList           = {'corg'};                                                  % End-member that should contribute to the numerator of C
    RZC_Denominator_EMList       = {'carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri'};    % end-members that should contribute to the denominator of R, Z, and C (weathering end-members)
    RZC_Denominator_IonList      = {'Na','Ca','Mg','K'};                                      % ions that should contribute to the denominator of R, Z, and C. we recommend including Na and K.                              
  % group 6/6, setting specific contributions from organic carbon
    ResetDegasDICContribution    = 1;
    DegasDICContributionMin      = -2.5;  % fraction relative to the DIC in solution, not the normalization
    DegasDICContributionMax      = +0.00; % fraction relative to the DIC in solution, not the normalization
    MinFractionalContribution0   = MinFractionalContribution;
    MaxFractionalContribution0   = MaxFractionalContribution; 

    
% scenario 5: AK_scenario5_carbonate_slctindi_degas_2p5_lowfrac
elseif isequal(ScenarioList{ScenarioListIndex},'AK_scenario5_carbonate_slctindi_degas_2p5_lowfrac')     
  % group 1/6, river observations - 12 variables
    Riverdatasource              = 'PCKAlaskaData';                                           % name for the source of the river observations
    AdjustRiverObs               = 1;                                                         % whether river measurements should be adjusted to reflect analytical uncertainty. options are: 0 or 1   
    ObsList                      = {'Ca', 'Mg', 'Na', 'K',  'Cl', 'SO4','DIC','d34S','d13C'}; % list of elements, ions, and isotopic ratios to include in the inversion
    CostFunType                  = {'rel','rel','rel','rel','rel','rel','rel','rel', 'rel'};  % whether each variable should be evaluated as proportional cost or absolute cost. options are 'rel' and 'abs'.
    WeightingList                = [1      1     1     1     1     1     1      1      1];    % weighting term for the cost function                     
    ErrorCutMinMB                = [95     95    95    95    95    95    95    -1     -1];    % if IterateOver is "Samples", these are the bounds of mass balance that will be accepted for an inversion result  
    ErrorCutMaxMB                = [105    105   105   105   105   105   105   +1     +1];    % if IterateOver is "Samples", these are the bounds of mass balance that will be accepted for an inversion result             
    nCFList                      = {};                                                        % list of dissolved variables (in ObsList) to not include in the cost function evaluation
    ConvertDelta2RList           = {'d13C','d34S'};                                           % list of isotopic variables (in delta notation) to be converted to isotopic ratios for the analysis (results are converted back to delta)
    AllIonsExplicitlyResolved    = 0;                                                         % Indicate if all cations and anions (all sources of charge) are resolved. options are: 0 or 1    
    ObsInNormalization           = {'Ca','Mg','Na','K','SO4','DIC'};                          % list of variables included in the normalization
    ImposeNormalizationCheck     = 1;                                                         % whether the normalization should be reproduced to the minimum and maximum values for non-isotopic variables 
  % group 2/6, general settings - 8 variables
    Solver                       = 'mldivide_optimize';                                       % type of inversion solution, options are: 'mldivide', 'lsqnonneg', 'mldivide_optimize', 'lsqnonneg_optimize', 'optimize'. enter as string, not a cell.
    IterateOver                  = 'Samples';                                                 % whether to apply the same end-members to all sample ('End-members'), or different end-members to each sample ('Sample'). options are: 'End-members', 'Samples'
    maxiterations                = 5e5;                                                       % if IterateOver equals "Samples", this is the maximum number of inversion attempts     
    maxsuccess                   = 200;                                                       % if IterateOver is "Samples", this is the desired number of successes      
    maxzerohits                  = 1e4;                                                       % if IterateOver is "Samples", this is the maximum number of iterations to try without success before moving on
    numberiterations             = NaN;                                                       % if IterateOver equals "End-members", this is the number of inversions to perform
    MisfitCuts                   = NaN;                                                       % if IterateOver equals "End-members", this percentage of simulations with lowest misfit between model result and data will be kept    
    CullOn                       = 'EachSample';                                              % if IterateOver is "End-members", this must be either 'EachSample' or 'AllSample' and describes HOW to cull the data
    saveuncutdata                = 0;                                                         % if 1, saves all the simulation data (results in very large files)    
  % group 3/6, end-members - 14 variables
    EMdatasource                 = 'PCK23_Alaska_SumCatSO4DIC_uniform_lowfrac';              % the name of the end-member group (the entry on column A of the spreadsheet containing end-member information). 
    EMList0                      = {'prec','carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri','evap','corg','degas'};  % list of the end-members to use in the inversion (row 2 of spreadsheet containing end-member information).     
    MinFractionalContribution    = [0       0       0          0         0        0        0       0      0      -9.99]; % the minumum fractional contribution of each end-member to the normalization ion (typically 0, unless secondary mineral formation is in inversion)
    MaxFractionalContribution    = [1       inf     1          1         1        1        1       1      inf    -9.99]; % the maximum fractional contribution of each end-member to the normalization ion (typically 1, unless secondary mineral formation is in inversion) 
    ListNormClosure              = {'Ca','Ca','Ca','Mg','Na','K','SO4','Ca','DIC','DIC'};     % the variable for each end-member that will be calculated by mass balance to ensure interal consistency
    ListChargeClosure            = {};                                                        % the variable for each end-member that will be calculated by charge balance to ensure internal consistency 
    EMUnits                      = 'equi';                                                    % whether the end-member chemistry is given as concentration or charge-equivalent. options are: 'conc', 'equi'    
    EMsources                    = {'prec','carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri','evap','corg'};  % end-members that count as sources of dissolved constituents (typically rocks and other solute sources)
    EMsinks                      = {};                                                        % end-members that count as sinks of dissolved constituents (clays and other secondary phases)
    EndMembersWithNegativeRatios = {};                                                        % end-members with negative chemical ratios    
    CoupleFeS2SO4intoEM          = {};                                                        % list of end-members where SO4 ratios represent pyrite
    CoupleFeS2d34SintoEM         = {};                                                        % list of end-members where SO4 d34S values represent pyrite
    RecordFullFeS2Distribution   = 0;                                                         % whether to record the full distribution of calculated FeS2 values, when d34S is not included in the inversion
    BalanceEvaporite             = 1;                                                         % whether evaporite SO4 = Ca+Mg+Sr and Cl = Na+K.    
  % group 4/6, Cl correction - 2 variables
    PrecProcessing               = 'EndMember';                                               % whether Cl should be treated as a common ion ('EndMember') or if there are Cl critical values ('ClCrit'). options are: 'EndMember', 'ClCrit'                            
    ClCriticalValuesGiven        = 1;                                                         % whether Cl critical values are given. if PrecProcessing is 'ClCrit' and ClCriticalValuesGiven is 1, those ClCritical values are used. if PrecProcessing is 'ClCrit' and ClCriticalValuesGiven is 0, 100% of river Cl is used.        
  % group 5/6, calculation of RZWYC - 8 variables
    CalculateRZCWY               = 1;                                                         % whether model should attempt to calculate R, Z, C, W, and Y. options are: 0 or 1
    R_Numerator_EMList           = {'carb'};                                                  % end-members that should contribute to the numerator of R
    R_Numerator_IonList          = {'Na','Ca','Mg','K'};                                      % ions that should contribute to the numerator of R        
    Z_NumeratorType              = {'ZfromEM'};                                               % how the numerator of Z (SO4 from FeS2) should be calculated. options are: 'ZfromEM', 'ZfromSO4excess', 'ZfromriverSO4', 'Znotcalculated'
    Z_Numerator_EMList           = {'pyri'};                                                  % end-members that should contribute to the numerator of Z (only used if Z_NumeratorType = {'ZfromEM'}).       
    C_Numerator_EMList           = {'corg'};                                                  % End-member that should contribute to the numerator of C
    RZC_Denominator_EMList       = {'carb','slct_Ca','slct_Mg','slct_Na','slct_K','pyri'};    % end-members that should contribute to the denominator of R, Z, and C (weathering end-members)
    RZC_Denominator_IonList      = {'Na','Ca','Mg','K'};                                      % ions that should contribute to the denominator of R, Z, and C. we recommend including Na and K.                              
  % group 6/6, setting specific contributions from organic carbon
    ResetDegasDICContribution    = 1;
    DegasDICContributionMin      = -2.5;  % fraction relative to the DIC in solution, not the normalization
    DegasDICContributionMax      = +0.00; % fraction relative to the DIC in solution, not the normalization
    MinFractionalContribution0   = MinFractionalContribution;
    MaxFractionalContribution0   = MaxFractionalContribution; 

else               
    disp(sprintf('catastrophic failure! no parameter values associated with simulation "%s"! check for typos in "ScenarioList" and check that the inversion name matches a set of parameter choices in "MEANDIR_FindIterationParameters"',ScenarioList{ScenarioListIndex}));
    FailCase = 1;
end


%% (2) define the sum of ions in normalization
% read in information on the conversion factors of concentration to equivalents 
[~, ~, all] = xlsread('MEANDIR_UserEntries','MEANDIR_conc2equi'); 
header      = all(1,:);                               for i=1:length(header);    if isnan(header{i});    header{i}    = 'NaN'; end; end
variable    = all(2:end,ismember(header,'Variable')); for i=1:length(variable);  if isnan(variable{i});  variable{i}  = 'NaN'; end; end
charge      = all(2:end,ismember(header,'Charge'));   for i=1:length(charge);    if isnan(charge{i});    charge{i}    = 'NaN'; end; end
% using the information from the sheet "MEANDIR_conc2equi", make the vector IonCharges 
% that contains information on the charge of each variable in the inversion
IonCharges  = {};
for j=1:length(ObsList)      
    acharge = charge{ismember(variable,ObsList{j})};    
    if     isequal(acharge,'positive'); IonCharges{1,j} = '+'; 
    elseif isequal(acharge,'negative'); IonCharges{1,j} = '-'; 
    elseif isequal(acharge,'neutral');  IonCharges{1,j} = '0'; 
    else                                IonCharges{1,j} = 'NaN'; 
    end   
end

% break up the normalization observations into cations, anions, and neutral species
allpos = variable(ismember(charge,'positive'))';                      % find all positive variables in charge
allneg = variable(ismember(charge,'negative'))';                      % find all negative variables in charge
allneu = variable(ismember(charge,'neutral'))';                       % find all neutral variables in charge
CationsInNormalization = allpos(ismember(allpos,ObsInNormalization)); % find all positive variables in the normalization
AnionsInNormalization  = allneg(ismember(allneg,ObsInNormalization)); % find all negative variables in the normalization
NeutralInNormalization = allneu(ismember(allneu,ObsInNormalization)); % find all neutral variables in the normalization
ObsInNormalization     = ObsInNormalization';

% define the normalization type
if  length(ObsInNormalization) == 1
    NormalizationType = ObsInNormalization{1};
else
    NormalizationType = 'SumObs';
    % the value "SumObs" for "NormalizationType" tells MEANDIR 
    % that multiple dissolved variables contribute to the normalization
end


%% (3) the following section runs 44 tests to identify problems with inversion criteria

% 1. check that the solver is one of the 5 acceptable options: mldivide, lsqnonneg, mldivide_optimize, lsqnonneg_optimize, or optimize
if ~(isequal(Solver,'mldivide') | isequal(Solver,'lsqnonneg') | isequal(Solver,'mldivide_optimize') | isequal(Solver,'lsqnonneg_optimize') | isequal(Solver,'optimize'))
   s0 = sprintf('catastrophic failure on test 1 for simulation "%s"! the "Solver" variable is currently "%s" but must be one of the following five options: \n', ScenarioList{ScenarioListIndex},Solver);
   s1 = '     (1) "mldivide"                       : least-squares inversion with a cost function of absolute misfit between model results and observations \n';
   s2 = '     (2) "mldivide_optimize"              : cost function is constrained misfit between model results and observations, with inital condition from least-squares inversion \n';
   s3 = '     (3) "lsqnonneg"                      : non-negative least-squares inversion with a cost function of absolute misfit between model results and observations \n';
   s4 = '     (4) "lsqnonneg_optimize"             : cost function is constrained misfit between model results and observations, with inital condition from non-negative least-squares inversion \n';      
   s5 = '     (5) "optimize"                       : cost function is constrained misfit between model results and observations, with inital condition of 1/n contribution from each of the n end-members \n';
   disp(sprintf([s0 s1 s2 s3 s4 s5]));
   FailCase = 1; 
end

% 2. test that IterateOver is either "Samples" or "End-members"
if ~(isequal(IterateOver,'Samples') | isequal(IterateOver,'End-members'))
    if isnumeric(IterateOver)
        IterateOver = 'NaN';
    end
    disp(sprintf('catastrophic failure for test 2 in simulation "%s"! the parameter "IterateOver" must be equal to either "Samples" or "End-members", but the current value is "%s". this variable controls how end-member compositions are presented to samples.',ScenarioList{ScenarioListIndex},IterateOver));
    FailCase = 1;
end

% 3. test that maxiterations is greater or equal to maxsuccess
if isequal(IterateOver,'Samples') & maxsuccess>maxiterations
   disp(sprintf('catastrophic failure for test 3 in simulation "%s"! the parameter "maxiterations" (current value: %i) must be greater or equation to the parameter "maxsuccess" (current value: %i).',ScenarioList{ScenarioListIndex},maxiterations,maxsuccess)); 
   FailCase = 1; 
end

% 4. check the values of MisfitCuts, ErrorCutMinMB, and ErrorCutMaxMB
if isequal(IterateOver,'End-members') & (MisfitCuts<0 | MisfitCuts>100 | isnan(MisfitCuts)) 
   disp(sprintf('catastrophic failure for test 4a in simulation "%s"! when IterateOver is "End-members", the parameter "MisfitCuts" must be >= 0 and <= 100, but the current value is %s. this parameter controls the fractions of simulations kept on the basis of misfit between model results and observations.',ScenarioList{ScenarioListIndex},num2str(round(MisfitCuts,0))));
   FailCase = 1;  
end   
isotopeposition = (ismember(ObsList,'d34S')   | ismember(ObsList,'d18O')  | ismember(ObsList,'Sr8786') | ismember(ObsList,'d7Li')   | ismember(ObsList,'d26Mg') | ismember(ObsList,'d30Si')  | ...
                   ismember(ObsList,'d42Ca')  | ismember(ObsList,'d44Ca') | ismember(ObsList,'d56Fe')  | ismember(ObsList,'d98Mo')  | ismember(ObsList,'Os8788') | ismember(ObsList,'d13C')  | ismember(ObsList,'Fmod'));
if isequal(IterateOver,'Samples') & (sum(ErrorCutMinMB(~isotopeposition)>100)>0 | sum(ErrorCutMaxMB(~isotopeposition)<100)>0)
    disp(sprintf('catastrophic failure for test 4b in simulation "%s"! when IterateOver is "Samples", the values in "ErrorCutMinMB" must be <= 100 and the values in "ErrorCutMaxMB" must be >= 100. these parameters define the range of mass balances acceptable in simulation results.',ScenarioList{ScenarioListIndex})); 
    FailCase = 1; 
end    
if isequal(IterateOver,'Samples') & (sum(isnan(ErrorCutMinMB))>0 | sum(isnan(ErrorCutMaxMB))>0)
    disp(sprintf('catastrophic failure for test 4c in simulation "%s"! when IterateOver is "Samples", the parameter "ErrorCutMinMB" must be <= 100 and the parameter "ErrorCutMaxMB" must be >= 100. currently either "ErrorCutMinMB" or "ErrorCutMaxMB" appears to contain a NaN value. these parameters define the range of mass balances acceptable in simulation results.',ScenarioList{ScenarioListIndex})); 
    FailCase = 1; 
end   

% 5. check the value of AdjustRiverObs 
if ~(AdjustRiverObs ==0 | AdjustRiverObs  == 1)
    if isnumeric(AdjustRiverObs)
        AdjustRiverObs = num2str(AdjustRiverObs);
    end
    disp(sprintf('catastrophic failure for test 5 in simulation "%s"! "AdjustRiverObs" must have a value of either 0 or 1 (numbers, not strings). the variable "AdjustRiverObs" currently has the value: %s.',ScenarioList{ScenarioListIndex},AdjustRiverObs)); 
    FailCase = 1;
end

% 6. test that entries in CationsInNormalization, AnionsInNormalization, and NeutralInNormalization are also found in ObsList
if isequal(NormalizationType,'SumObs')
    if ~isequal(CationsInNormalization{1},'NaN')
        if sum(~ismember(CationsInNormalization,ObsList))>0
            i  = find(~ismember(CationsInNormalization,ObsList)); 
            disp(sprintf('catastrophic failure for test 6a in simulation "%s"! parameter "%s" is in "CationsInNormalization" but not "ObsList".',ScenarioList{ScenarioListIndex},CationsInNormalization{i}));
            FailCase = 1; 
        end
    end    
    if ~isempty(AnionsInNormalization)
        if sum(~ismember(AnionsInNormalization,ObsList))>0
            i = find(~ismember(AnionsInNormalization,ObsList)); 
            disp(sprintf('catastrophic failure for test 6b in simulation "%s"! parameter "%s" is in "AnionsInNormalization" but not "ObsList".',ScenarioList{ScenarioListIndex},AnionsInNormalization{i}));
            FailCase = 1; 
        end    
    end    
    if ~isempty(NeutralInNormalization)
        if sum(~ismember(NeutralInNormalization,ObsList))>0
            i = find(~ismember(NeutralInNormalization,ObsList)); 
            disp(sprintf('catastrophic failure for test 6c in simulation "%s"! parameter "%s" is in "NeutralInNormalization" but not "ObsList".',ScenarioList{ScenarioListIndex},NeutralInNormalization{i}));
            FailCase = 1; 
        end    
    end    
end

% 7. test that entries of CationsInNormalization, AnionsInNormalization, and NeutralInNormalization are not duplicated
if length(unique(CationsInNormalization))~=length(CationsInNormalization);  disp(sprintf('catastrophic failure for test 7a! the number of unique entries in "CationsInNormalization" does not match its length. check for duplicate entries in "CationsInNormalization" in simulation "%s"', ScenarioList{ScenarioListIndex})); FailCase = 1; end
if length(unique(AnionsInNormalization)) ~=length(AnionsInNormalization);   disp(sprintf('catastrophic failure for test 7b! the number of unique entries in "AnionsInNormalization" does not match its length. check for duplicate entries in "AnionsInNormalization" in simulation "%s"',   ScenarioList{ScenarioListIndex})); FailCase = 1; end
if length(unique(NeutralInNormalization))~=length(NeutralInNormalization);  disp(sprintf('catastrophic failure for test 7c! the number of unique entries in "NeutralInNormalization" does not match its length. check for duplicate entries in "NeutralInNormalization" in simulation "%s"', ScenarioList{ScenarioListIndex})); FailCase = 1; end
if length(unique(ObsList))~=length(ObsList);                                disp(sprintf('catastrophic failure for test 7d! the number of unique entries in "ObsList" does not match its length. problem with "ObsList" in simulation "%s", check for duplicate entries',                    ScenarioList{ScenarioListIndex})); FailCase = 1; end

% 8. test if the requested end-member ratios can be found
EM = MEANDIR_ReadEndMembersIntoMatlab;    
for i=1:length(ObsList)    
for j=1:length(EMList0)         
    t1 = sprintf('EM.%s.%s.NOR_Men%s%s;',EMdatasource,EMList0{j},ObsList{i},NormalizationType);
    t2 = sprintf('EM.%s.%s.UNI_Min%s%s;',EMdatasource,EMList0{j},ObsList{i},NormalizationType);
    t3 = sprintf('EM.%s.%s.LGU_Min%s%s;',EMdatasource,EMList0{j},ObsList{i},NormalizationType);    
    try eval(t1);
    catch
        try eval(t2);
        catch
            try eval(t3);
            catch
                disp(sprintf('catastrophic failure for test 8 in simulation "%s"! there is no variable %s normalized by %s for end-member %s in end-member set %s. check the excel sheet of end-member compositions and the names of variables in "ObsList".',ScenarioList{ScenarioListIndex},ObsList{i},NormalizationType,EMList0{j},EMdatasource))
                FailCase = 1; 
            end
        end
    end
    
end
end

% 9. test the selected error style for the end-members
if ~(isequal(EMUnits,'conc') | isequal(EMUnits,'equi'))
    disp(sprintf('catastrophic failure for test 9 in simulation %s! "EMUnits" is currently "%s" but must be either "conc" or "equi"', ScenarioList{ScenarioListIndex},EMUnits)); 
    FailCase = 1;
end

% 10. check for duplicate entries to EMList0
if length(unique(EMList0))~=length(EMList0); 
    disp(sprintf('catastrophic failure for test 10! problem with "EMList0" in simulation "%s", check for duplicate entries',     ScenarioList{ScenarioListIndex}));  
    FailCase = 1;
end

% 11. check if the PrecProcessing entry is acceptable
if ~(isequal(PrecProcessing,'ClCrit') | isequal(PrecProcessing,'EndMember'))
    disp(sprintf('catastrophic failure for test 11 in simulation %s! "PrecProcessing" is currently "%s" but must be either "ClCrit" or "EndMember"', ScenarioList{ScenarioListIndex},PrecProcessing));
    FailCase = 1;
end

% 12. check the value of ClCriticalValuesGiven
if ~(ClCriticalValuesGiven==0 | ClCriticalValuesGiven == 1)
    if isnumeric(ClCriticalValuesGiven)
        ClCriticalValuesGiven = num2str(ClCriticalValuesGiven);
    end
    disp(sprintf('catastrophic failure for test 12 in simulation "%s"! "ClCriticalValuesGiven" must have a value of either 0 or 1 (numbers, not strings). the current value is: %s',ScenarioList{ScenarioListIndex},ClCriticalValuesGiven)); 
    FailCase = 1;
end

%  13. check if CalculateRZCWY is 0 or 1
if ~(CalculateRZCWY ==0 | CalculateRZCWY == 1)
    if isnumeric(CalculateRZCWY)
        CalculateRZCWY = num2str(CalculateRZCWY);
    end        
    disp(sprintf('catastrophic failure for test 13 in simulation "%s"! "CalculateRZ" must have a value of either 0 or 1 (numbers, not strings). the current value is: %s',ScenarioList{ScenarioListIndex},CalculateRZCWY)); 
    FailCase       = 1;
    CalculateRZCWY = NaN; 
end

if CalculateRZCWY== 1
    % 14. test that all cations in sublists are not duplicates and occur in ObsList
        if length(unique(R_Numerator_IonList))~=length(R_Numerator_IonList);                                                     disp(sprintf('catastrophic failure for test 14a! problem with "R_Numerator_IonList" in simulation "%s", check for duplicate entries',       ScenarioList{ScenarioListIndex}));                                 FailCase = 1; end
        if sum(~ismember(R_Numerator_IonList,ObsList))>0;                 i = find(~ismember(R_Numerator_IonList,ObsList));      disp(sprintf('catastrophic failure for test 14b in simulation "%s"! Variable "%s" is in "R_Numerator_IonList" but not "ObsList"',           ScenarioList{ScenarioListIndex},R_Numerator_IonList{min(i)}));     FailCase = 1; end
        if length(unique(RZC_Denominator_IonList))~=length(RZC_Denominator_IonList);                                             disp(sprintf('catastrophic failure for test 14c! problem with "RZ_Denominator_IonList" in simulation "%s", check for duplicate entries \n', ScenarioList{ScenarioListIndex}));                                 FailCase = 1; end
        if sum(~ismember(RZC_Denominator_IonList,ObsList))>0;              i = find(~ismember(RZC_Denominator_IonList,ObsList)); disp(sprintf('catastrophic failure for test 14d in simulation "%s"! Variable "%s" is in "RZ_Denominator_IonList" but not "ObsList" \n',     ScenarioList{ScenarioListIndex},RZC_Denominator_IonList{min(i)})); FailCase = 1; end    
    
    % 15. test that all end-members in sublists are also found in EMList0 and are not duplicates
        if sum(~ismember(R_Numerator_EMList,EMList0))    >0; i = find(~ismember(R_Numerator_EMList,EMList0));     disp(sprintf('catastrophic failure for test 15a in simulation "%s"! end-member "%s" is in "R_Numerator_EMList" but not "EMList0"',   ScenarioList{ScenarioListIndex},R_Numerator_EMList{i}));     FailCase = 1; end
        if sum(~ismember(RZC_Denominator_EMList,EMList0))>0; i = find(~ismember(RZC_Denominator_EMList,EMList0)); disp(sprintf('catastrophic failure for test 15b in simulation "%s"! end-member "%s" is in "RZ_Denominator_EMList" but not "EMList0"',ScenarioList{ScenarioListIndex},RZC_Denominator_EMList{i})); FailCase = 1; end
        if sum(~ismember(Z_Numerator_EMList,EMList0))    >0; i = find(~ismember(Z_Numerator_EMList,EMList0));     disp(sprintf('catastrophic failure for test 15c in simulation "%s"! end-member "%s" is in "Z_Numerator_EMList" but not "EMList0"',   ScenarioList{ScenarioListIndex},Z_Numerator_EMList{i}));     FailCase = 1; end
        if sum(~ismember(C_Numerator_EMList,EMList0))    >0; i = find(~ismember(C_Numerator_EMList,EMList0));     disp(sprintf('catastrophic failure for test 15d in simulation "%s"! end-member "%s" is in "C_Numerator_EMList" but not "EMList0"',   ScenarioList{ScenarioListIndex},C_Numerator_EMList{i}));     FailCase = 1; end
        if length(unique(R_Numerator_EMList))    ~=length(R_Numerator_EMList);                                    disp(sprintf('catastrophic failure for test 15e! "R_Numerator_EMList" in simulation "%s", may contain duplicate values' ,            ScenarioList{ScenarioListIndex}));                           FailCase = 1; end
        if length(unique(RZC_Denominator_EMList))~=length(RZC_Denominator_EMList);                                disp(sprintf('catastrophic failure for test 15f! "RZ_Denominator_EMList" in simulation "%s", may contain duplicate values',          ScenarioList{ScenarioListIndex}));                           FailCase = 1; end
        if length(unique(Z_Numerator_EMList))    ~=length(Z_Numerator_EMList);                                    disp(sprintf('catastrophic failure for test 15g! "ZnumEM" in simulation "%s", may contain duplicate values',                         ScenarioList{ScenarioListIndex}));                           FailCase = 1; end
        if length(unique(C_Numerator_EMList))    ~=length(C_Numerator_EMList);                                    disp(sprintf('catastrophic failure for test 15h! "CnumEM" in simulation "%s", may contain duplicate values',                         ScenarioList{ScenarioListIndex}));                           FailCase = 1; end        
    
    % 16. check Z numerator   
        ZfromSO4excess = 0; ZfromriverSO4  = 0; ZfromEM  = 0; Znotcalculated = 0;     
        if length(Z_NumeratorType)~=1; disp(sprintf('catastrophic failure for test 16a in simulation "%s"! Z_NumeratorType must only have 1 value from the following options: "ZfromEM", "ZfromSO4excess", "ZfromriverSO4", "Znotcalculated".',ScenarioList{ScenarioListIndex})); end
        if     isequal(Z_NumeratorType{1},'ZfromEM');        ZfromEM        = 1;
        elseif isequal(Z_NumeratorType{1},'ZfromSO4excess'); ZfromSO4excess = 1;
        elseif isequal(Z_NumeratorType{1},'ZfromriverSO4');  ZfromriverSO4  = 1;
        elseif isequal(Z_NumeratorType{1},'Znotcalculated'); Znotcalculated = 1;
        else   disp(sprintf('catastrophic failure for test 16b in "%s"! "Z_NumeratorType" is "%s", but must be equal to either "ZfromEM", "ZfromSO4excess", "ZfromriverSO4", or "Znotcalculated". Check the parameter settings in "MEANDIR_FindIteration"',ScenarioList{ScenarioListIndex},Z_NumeratorType{1}))
        end
        
    % 17. check that all elements of R are also in the denominators
        if sum(~ismember(R_Numerator_IonList,RZC_Denominator_IonList))>0; i = find(~ismember(R_Numerator_IonList,RZC_Denominator_IonList));  disp(sprintf('catastrophic failure for test 17a in simulation "%s"! parameter "%s" is in "R_Numerator_IonList" but not "RZ_Denominator_IonList"', ScenarioList{ScenarioListIndex},R_Numerator_IonList{min(i)})); FailCase = 1; end
        if sum(~ismember(R_Numerator_EMList,RZC_Denominator_EMList))>0;   i = find(~ismember(R_Numerator_EMList,RZC_Denominator_EMList));    disp(sprintf('catastrophic failure for test 17b in simulation "%s"! parameter "%s" is in "R_Numerator_EMList" but not "RZ_Denominator_EMList"',   ScenarioList{ScenarioListIndex},R_Numerator_EMList{min(i)}));  FailCase = 1; end
else    
    R_Numerator_EMList      = {'NaN'};                    
    R_Numerator_IonList     = {'NaN'};
    RZC_Denominator_EMList  = {'NaN'};
    RZC_Denominator_IonList = {'NaN'};
    Z_NumeratorType         = {'NaN'};
    Z_Numerator_EMList      = {'NaN'};
    ZfromEM                 = 0; 
    ZfromSO4excess          = 0; 
    ZfromriverSO4           = 0;
    Znotcalculated          = 0; 
end

% 18. test CoupleFeS2SO4intoEMSO4 and CoupleFeS2d34SintoEMSO4
if length(unique(CoupleFeS2SO4intoEM)) ~=length(CoupleFeS2SO4intoEM);  disp(sprintf('catastrophic failure for test 18a! problem with "CoupleFeS2SO4intoEM" in simulation "%s". the number of unique entries does not match the length of "CoupleFeS2SO4intoEM", check for duplicate entries',   ScenarioList{ScenarioListIndex})); FailCase = 1; end
if length(unique(CoupleFeS2d34SintoEM))~=length(CoupleFeS2d34SintoEM); disp(sprintf('catastrophic failure for test 18b! problem with "CoupleFeS2d34SintoEM" in simulation "%s". the number of unique entries does not match the length of "CoupleFeS2d34SintoEM", check for duplicate entries', ScenarioList{ScenarioListIndex})); FailCase = 1; end
if sum(~ismember(CoupleFeS2SO4intoEM,EMList0))>0;  i = find(~ismember(CoupleFeS2SO4intoEM,EMList0)); disp(sprintf('catastrophic failure for test 18c in simulation "%s"! End-member "%s" is in "CoupleFeS2SO4intoEM" but not "EMList0"',   ScenarioList{ScenarioListIndex},CoupleFeS2SO4intoEM{min(i)}));  FailCase = 1; end
if sum(~ismember(CoupleFeS2d34SintoEM,EMList0))>0; i = find(~ismember(CoupleFeS2d34SintoEM,EMList0)); disp(sprintf('catastrophic failure for test 18d in simulation "%s"! End-member "%s" is in "CoupleFeS2d34SintoEM" but not "EMList0"', ScenarioList{ScenarioListIndex},CoupleFeS2d34SintoEM{min(i)})); FailCase = 1; end

% 19. check if HCO3 and DIC are being requested
if sum(ismember(ObsList,'DIC'))>0 & sum(ismember(ObsList,'HCO3'))>0
FailCase = 1; disp(sprintf('catastrophic failure for test 19 in simulation "%s"! "ObsList" cannot contain both DIC and HCO3. choose 1 way of representing carbon in rivers.',ScenarioList{ScenarioListIndex}));      
carbonisotopematch = 'NaN';
else
   if   sum(ismember(ObsList,'DIC'))>0; carbonisotopematch = 'DIC';
   else;                                carbonisotopematch = 'HCO3';
   end
end

% 20. test if isotopes are matched by elements
if sum(ismember(ObsList,'d34S'))  >0 & sum(ismember(ObsList,'SO4'))==0;                                       FailCase = 1; disp(sprintf('catastrophic failure for test 20a in simulation "%s"! isotope d34S is included in inversion without SO4, but the isotope mass balance equation requires both',ScenarioList{ScenarioListIndex}));   end
if sum(ismember(ObsList,'d7Li'))  >0 & sum(ismember(ObsList,'Li'))==0;                                        FailCase = 1; disp(sprintf('catastrophic failure for test 20b in simulation "%s"! isotope d7Li is included in inversion without Li, but the isotope mass balance equation requires both',ScenarioList{ScenarioListIndex}));    end
if sum(ismember(ObsList,'d26Mg')) >0 & sum(ismember(ObsList,'Mg'))==0;                                        FailCase = 1; disp(sprintf('catastrophic failure for test 20c in simulation "%s"! isotope d26Mg is included in inversion without Mg, but the isotope mass balance equation requires both',ScenarioList{ScenarioListIndex}));   end
if sum(ismember(ObsList,'d42Ca')) >0 & sum(ismember(ObsList,'Ca'))==0;                                        FailCase = 1; disp(sprintf('catastrophic failure for test 20d in simulation "%s"! isotope d42Ca is included in inversion without Ca, but the isotope mass balance equation requires both',ScenarioList{ScenarioListIndex}));   end
if sum(ismember(ObsList,'d44Ca')) >0 & sum(ismember(ObsList,'Ca'))==0;                                        FailCase = 1; disp(sprintf('catastrophic failure for test 20e in simulation "%s"! isotope d44Ca is included in inversion without Ca, but the isotope mass balance equation requires both',ScenarioList{ScenarioListIndex}));   end
if sum(ismember(ObsList,'d56Fe')) >0 & sum(ismember(ObsList,'Fe'))==0;                                        FailCase = 1; disp(sprintf('catastrophic failure for test 20f in simulation "%s"! isotope d56Fe is included in inversion without Fe, but the isotope mass balance equation requires both',ScenarioList{ScenarioListIndex}));   end
if sum(ismember(ObsList,'Sr8786'))>0 & sum(ismember(ObsList,'Sr'))==0;                                        FailCase = 1; disp(sprintf('catastrophic failure for test 20g in simulation "%s"! isotope Sr8786 is included in inversion without Sr, but the isotope mass balance equation requires both',ScenarioList{ScenarioListIndex}));  end
if sum(ismember(ObsList,'d98Mo')) >0 & sum(ismember(ObsList,'Mo'))==0;                                        FailCase = 1; disp(sprintf('catastrophic failure for test 20h in simulation "%s"! isotope d98Mo is included in inversion without Mo, but the isotope mass balance equation requires both',ScenarioList{ScenarioListIndex}));   end
if sum(ismember(ObsList,'Os8788'))>0 & sum(ismember(ObsList,'Os'))==0;                                        FailCase = 1; disp(sprintf('catastrophic failure for test 20i in simulation "%s"! isotope Os8788 is included in inversion without Os, but the isotope mass balance equation requires both',ScenarioList{ScenarioListIndex}));  end
if sum(ismember(ObsList,'d13C'))  >0 & isequal(carbonisotopematch,'HCO3') & sum(ismember(ObsList,'HCO3'))==0; FailCase = 1; disp(sprintf('catastrophic failure for test 20j in simulation "%s"! isotope d13C is included in inversion with HCO3 or DIC',ScenarioList{ScenarioListIndex}));                                                   end
if sum(ismember(ObsList,'Fmod'))  >0 & isequal(carbonisotopematch,'HCO3') & sum(ismember(ObsList,'HCO3'))==0; FailCase = 1; disp(sprintf('catastrophic failure for test 20k in simulation "%s"! isotope Fmod is included in inversion without HCO3 or DIC',ScenarioList{ScenarioListIndex}));                                                end

% 21. check length of ListNormClosure and ListChargeClosure
if length(ListNormClosure)  ~=length(EMList0) & isequal(NormalizationType,'SumObs');  FailCase = 1; disp(sprintf('catastrophic failure for test 21a in simulation "%s"! length of ListNormClosure does not match length of EMList0',ScenarioList{ScenarioListIndex}));                                       end
if length(ListChargeClosure)~=length(EMList0) & AllIonsExplicitlyResolved==1;         FailCase = 1; disp(sprintf('catastrophic failure for test 21b in simulation "%s"! length of ListChargeClosure does not match length of EMList0 while "AllIonsExplicitlyResolved"=1',ScenarioList{ScenarioListIndex})); end

% 22. check that ListNormClosure and ListChargeClosure are in ObsList
testListNormClosure   = ListNormClosure(~ismember(ListNormClosure,'NaN'));
testListChargeClosure = ListChargeClosure(~ismember(ListChargeClosure,'NaN'));
if sum(ismember(testListNormClosure,ObsList))~=length(testListNormClosure);     FailCase = 1; disp(sprintf('catastrophic failure for test 22a in simulation "%s"! entries of "ListNormClosure" must all be entries of "ObsList"',ScenarioList{ScenarioListIndex}));   end
if sum(ismember(testListChargeClosure,ObsList))~=length(testListChargeClosure); FailCase = 1; disp(sprintf('catastrophic failure for test 22b in simulation "%s"! entries of "ListChargeClosure" must all be entries of "ObsList"',ScenarioList{ScenarioListIndex})); end

% 23. check that elements of nCFList are in ObsList
if sum(ismember(nCFList,ObsList)==0)>0
    val = nCFList(~ismember(nCFList,ObsList));
    val = val{1};
    FailCase = 1; disp(sprintf('catastrophic failure for test 23 in simulation "%s"! element "%s" in nCFList is not in ObsList. ',ScenarioList{ScenarioListIndex},val));
end

% 24. if IterateOver is EndMembers, check that CullOn is acceptable
if  isequal(IterateOver,'End-members')
       if ~(isequal(CullOn,'EachSample') | isequal(CullOn,'AllSample'))
          FailCase = 1; 
          disp(sprintf('catastrophic failure for test 24 in simulation "%s"! when IterateOver is End-members, "CullOn" must either "EachSample" or "AllSample"',ScenarioList{ScenarioListIndex}));        
       end
elseif isequal(IterateOver,'Samples')
       CullOn = NaN; 
end

% 25. if dealing with normalization and charge, check that the charge terms are not in the normalization
if isequal(NormalizationType,'SumObs') & AllIonsExplicitlyResolved == 1
    if ~isempty(intersect(ObsInNormalization,ListChargeClosure))       
        l = intersect(ObsInNormalization,ListChargeClosure);
        lp = ''; 
        for j=1:length(l)
            lp = [lp l{j} ', ']; 
        end
        lp = lp(1:end-2); % This removes the final ', ' from the list to make the display nicer
        FailCase = 1; disp(sprintf('catastrophic failure for test 25 in simulation "%s"! when "AllIonsExplicitlyResolved"=1, "ListChargeClosure" must only contain variables not included in the normalization. replace ion(s): %s in "ListChargeClosure".',ScenarioList{ScenarioListIndex},lp));
    end
end

% 26. check that elements of ListNormClosure are in ObsInNormalization
if sum(ismember(ListNormClosure,ObsInNormalization))~=length(ListNormClosure) & isequal(NormalizationType,'SumObs')
   FailCase = 1; disp(sprintf('catastrophic failure for test 26 in simulation "%s"! all elements of "ListNormClosure" must appear in "ObsInNormalization".',ScenarioList{ScenarioListIndex}));
end

% 27. check that ion charges are +, -, or 0, and balanced
if BalanceEvaporite == 1 | AllIonsExplicitlyResolved == 1
    if  length(IonCharges)~=length(ObsList)
        FailCase = 1; disp(sprintf('catastrophic failure for test 27a in simulation "%s"! the length of "IonCharges" must equal the length of "ObsList" when either "BalanceEvaporite"=1 or "AllIonsExplicitlyResolved"=1',ScenarioList{ScenarioListIndex}));
    end
    if  sum(ismember(IonCharges,'+') | ismember(IonCharges,'-') | ismember(IonCharges,'0'))~=length(IonCharges)
        FailCase = 1; disp(sprintf('catastrophic failure for test 27b in simulation "%s"! "IonCharges" may only contain "+", "-", and "0"',ScenarioList{ScenarioListIndex}));   
    end
else
    IonCharges = {};
end

% 28. check the length of ErrorCutMinMB and ErrorCutMaxMB
if isequal(IterateOver,'Samples')
    if length(ErrorCutMinMB)~=length(ObsList)       
       disp(sprintf('catastrophic failure for test 28a in simulation "%s"! the length of "ErrorCutMinMB" must equal the length of "ObsList"',ScenarioList{ScenarioListIndex}));
       FailCase = 1; 
    end
    if length(ErrorCutMaxMB)~=length(ObsList)       
       disp(sprintf('catastrophic failure for test 28b in simulation "%s"! the length of "ErrorCutMaxMB" must equal the length of "ObsList"',ScenarioList{ScenarioListIndex}));
       FailCase = 1; 
    end
end

% 29. check the length of WeightingList and length and content of CostFunType
if (isequal(Solver,'mldivide_optimize') | isequal(Solver,'lsqnonneg_optimize') | isequal(Solver,'optimize'))
    if length(WeightingList)~=length(ObsList)
       FailCase = 1; disp(sprintf('catastrophic failure for test 29a in simulation "%s"! the length of "WeightingList" must equal the length of "ObsList"',ScenarioList{ScenarioListIndex}));
    end
    if length(CostFunType)~=length(ObsList)
       FailCase = 1; disp(sprintf('catastrophic failure for test 29b in simulation "%s"! the length of "CostFunType" must equal the length of "ObsList"',ScenarioList{ScenarioListIndex}));
    end
    if sum(~ismember(CostFunType,'abs') & ~ismember(CostFunType,'rel'))>0
       badindx = find(~ismember(CostFunType,'abs') & ~ismember(CostFunType,'rel'));
       FailCase = 1; disp(sprintf('catastrophic failure for test 29c in simulation "%s"! "CostFunType" at position %i has value "%s", but must be either "rel" or "abs"',ScenarioList{ScenarioListIndex},badindx(1),CostFunType{badindx}));        
    end
end
if ~(isequal(Solver,'mldivide_optimize') | isequal(Solver,'lsqnonneg_optimize') | isequal(Solver,'optimize'))
   if ~isempty(WeightingList); disp(sprintf('in simulation "%s", the variable "WeightingList" has values but will not be used because the solver is not "mldivide_optimize", "lsqnonneg_optimize", or "optimize"',ScenarioList{ScenarioListIndex})); end
   if ~isempty(CostFunType);   disp(sprintf('in simulation "%s", the variable "CostFunType" has values but will not be used because the solver is not "mldivide_optimize", "lsqnonneg_optimize", or "optimize"',ScenarioList{ScenarioListIndex}));   end
end
    
% 30. check mass balance will be maintained even with evaporite balance
% this test is difficult to understand, but it is important to pass, so do not remove!
if BalanceEvaporite == 1 & sum(ismember(EMList0,'evap'))>0 & isequal(NormalizationType,'SumObs')
   evappos = ismember(EMList0,'evap');      
   if isequal(ListNormClosure{evappos},'SO4') | isequal(ListNormClosure{evappos},'Cl')
      FailCase = 1; disp(sprintf('catastrophic failure for test 30 in simulation "%s"! when "BalanceEvaporite"=1, the closure variable for the evaporite end-member cannot be SO4 or Cl',ScenarioList{ScenarioListIndex}));
   end   
   activeion = ListNormClosure{evappos};   
end

% 31. check is EMsources is empty
if isempty(EMsources) & CalculateRZCWY==1   
    EMsources = EMList0;
    disp('EMsources was empty! setting EMsources equal to EMList0, meaning the inversion will consider all end-members as sources.');
end

% 32. check entries of EMsources and EMsinks can be found in EMList0
if ~isempty(EMsources)
    if sum(~ismember(EMsources,EMList0))>0
           problemEM = EMsources{~ismember(EMsources,EMList0)};
           disp(sprintf('catastrophic failure for test 32a in simulation "%s"! "EMsources" is attempting to use end-member "%s", which is not in "EMList0"',ScenarioList{ScenarioListIndex},problemEM));    
           FailCase = 1;
    end
end
if ~isempty(EMsinks)
    if sum(~ismember(EMsinks,EMList0))>0
           problemEM = EMsinks{~ismember(EMsinks,EMList0)};
           disp(sprintf('catastrophic failure for test 32b in simulation "%s"! "EMsinks" is attempting to use end-member "%s", which is not in "EMList0"',ScenarioList{ScenarioListIndex},problemEM));    
           FailCase = 1;
    end
end
 
% 33. check if ConvertDelta2RList exists
if ~exist('ConvertDelta2RList','var')
   disp(sprintf('catastrophic failure for test 33 in simulation "%s"! "ConvertDelta2RList" does not exist as a variable!',ScenarioList{ScenarioListIndex}));    
   FailCase = 1;
end

% 34. check entries of ConvertDelta2RList occur in ObsList
ConvertDelta2RList = ConvertDelta2RList(~ismember(ConvertDelta2RList,' '));
ConvertDelta2RList = ConvertDelta2RList(~ismember(ConvertDelta2RList,''));
if ~isempty(ConvertDelta2RList)
    problist = ConvertDelta2RList(~ismember(ConvertDelta2RList,ObsList));
    if ~isempty(problist)
        for j=1:length(problist)
            disp(sprintf('catastrophic failure for test 34a in simulation "%s"! "ConvertDelta2RList" is attempting to convert isotopic information "%s", which is not in "ObsList"',ScenarioList{ScenarioListIndex},problist{j}));
        end
       FailCase = 1; 
    end
end
if sum(ismember(ConvertDelta2RList,'d7Li')) >0 & sum(ismember(ObsList,'Li')) ==0; disp(sprintf('catastrophic failure for test 34b in simulation "%s"! "ConvertDelta2RList" is attempting to convert isotopic information "d7Li" but "Li" is not in "ObsList"',ScenarioList{ScenarioListIndex}));  FailCase = 1; end
if sum(ismember(ConvertDelta2RList,'d18O')) >0 & sum(ismember(ObsList,'SO4'))==0; disp(sprintf('catastrophic failure for test 34c in simulation "%s"! "ConvertDelta2RList" is attempting to convert isotopic information "d18O" but "SO4" is not in "ObsList"',ScenarioList{ScenarioListIndex})); FailCase = 1; end
if sum(ismember(ConvertDelta2RList,'d26Mg'))>0 & sum(ismember(ObsList,'Mg')) ==0; disp(sprintf('catastrophic failure for test 34d in simulation "%s"! "ConvertDelta2RList" is attempting to convert isotopic information "d26Mg" but "Mg" is not in "ObsList"',ScenarioList{ScenarioListIndex})); FailCase = 1; end
if sum(ismember(ConvertDelta2RList,'d30Si'))>0 & sum(ismember(ObsList,'Si')) ==0; disp(sprintf('catastrophic failure for test 34e in simulation "%s"! "ConvertDelta2RList" is attempting to convert isotopic information "d30Si" but "Si" is not in "ObsList"',ScenarioList{ScenarioListIndex})); FailCase = 1; end
if sum(ismember(ConvertDelta2RList,'d34S')) >0 & sum(ismember(ObsList,'SO4'))==0; disp(sprintf('catastrophic failure for test 34f in simulation "%s"! "ConvertDelta2RList" is attempting to convert isotopic information "d34S" but "SO4" is not in "ObsList"',ScenarioList{ScenarioListIndex})); FailCase = 1; end
if sum(ismember(ConvertDelta2RList,'d42Ca'))>0 & sum(ismember(ObsList,'Ca')) ==0; disp(sprintf('catastrophic failure for test 34g in simulation "%s"! "ConvertDelta2RList" is attempting to convert isotopic information "d42Ca" but "Ca" is not in "ObsList"',ScenarioList{ScenarioListIndex})); FailCase = 1; end
if sum(ismember(ConvertDelta2RList,'d44Ca'))>0 & sum(ismember(ObsList,'Ca')) ==0; disp(sprintf('catastrophic failure for test 34h in simulation "%s"! "ConvertDelta2RList" is attempting to convert isotopic information "d44Ca" but "Ca" is not in "ObsList"',ScenarioList{ScenarioListIndex})); FailCase = 1; end
if sum(ismember(ConvertDelta2RList,'d56Fe'))>0 & sum(ismember(ObsList,'Fe')) ==0; disp(sprintf('catastrophic failure for test 34i in simulation "%s"! "ConvertDelta2RList" is attempting to convert isotopic information "d56Fe" but "Fe" is not in "ObsList"',ScenarioList{ScenarioListIndex})); FailCase = 1; end
if sum(ismember(ConvertDelta2RList,'d98Mo'))>0 & sum(ismember(ObsList,'Mo')) ==0; disp(sprintf('catastrophic failure for test 34j in simulation "%s"! "ConvertDelta2RList" is attempting to convert isotopic information "d98Mo" but "Mo" is not in "ObsList"',ScenarioList{ScenarioListIndex})); FailCase = 1; end

% 35. reshape the MinFractionalContribution and MaxFractionalContribution vectors
MinFractionalContribution  = reshape(MinFractionalContribution,length(EMList0),1);
MaxFractionalContribution  = reshape(MaxFractionalContribution,length(EMList0),1); 
MinFractionalContribution0 = reshape(MinFractionalContribution0,length(EMList0),1);
MaxFractionalContribution0 = reshape(MaxFractionalContribution0,length(EMList0),1); 

    
% 36. find the isotope data
isotopeposition = (ismember(ObsList,'d34S')  | ismember(ObsList,'d18O')  | ismember(ObsList,'Sr8786') | ismember(ObsList,'d7Li') | ismember(ObsList,'d26Mg') | ismember(ObsList,'d30Si')  |...
                   ismember(ObsList,'d42Ca') | ismember(ObsList,'d44Ca') | ismember(ObsList,'d56Fe')  | ismember(ObsList,'d98Mo') | ismember(ObsList,'Os8788') | ismember(ObsList,'d13C') | ismember(ObsList,'Fmod'));
indxisopos      = find(isotopeposition);    % find the positions of isotope data within ObsList

% 37. when IterateOver equals "End-members" define arbitrarily wide criteria for
% acceptance. when IterOver equals "Samples" set MisfitCuts to NaN.
if     isequal(IterateOver,'End-members')        
       ErrorCutMinMB(~isotopeposition) = -1e10;  ErrorCutMinMB(isotopeposition) = -1e10;
       ErrorCutMaxMB(~isotopeposition) = 1e10;   ErrorCutMaxMB(isotopeposition) =  1e10;
elseif isequal(IterateOver,'Samples')
       MisfitCuts = NaN;    
end

% 38. set empty charge-specific normalization vectors to NaN
if isempty(AnionsInNormalization);  AnionsInNormalization  = 'NaN'; end
if isempty(CationsInNormalization); CationsInNormalization = 'NaN'; end
if isempty(NeutralInNormalization); NeutralInNormalization = 'NaN'; end

% 39. check if SO4 analysis is possible
EMHaveSO4 = 1; 
for j=1:length(EMList0)
 t1 = sprintf('EM.%s.%s.NOR_Men%s%s',EMdatasource,EMList0{j},'SO4',NormalizationType);
 t2 = sprintf('EM.%s.%s.UNI_Min%s%s',EMdatasource,EMList0{j},'SO4',NormalizationType);
 t3 = sprintf('EM.%s.%s.LGU_Min%s%s',EMdatasource,EMList0{j},'SO4',NormalizationType);
 try eval([t1 ';']);
 catch 
     try eval([t2 ';']);
     catch
         try eval([t3 ';']);
         catch
             EMHaveSO4 = 0; 
         end
     end
 end
end

% 40. check if d34S analysis is possible
EMHaved34S = 1; 
for j=1:length(EMList0)     
 t1 = sprintf('EM.%s.%s.NOR_Men%s%s',EMdatasource,EMList0{j},'d34S',NormalizationType);
 t2 = sprintf('EM.%s.%s.UNI_Min%s%s',EMdatasource,EMList0{j},'d34S',NormalizationType);
 t3 = sprintf('EM.%s.%s.LGU_Min%s%s',EMdatasource,EMList0{j},'d34S',NormalizationType);
 try eval([t1 ';']);
 catch 
     try eval([t2 ';']);
     catch
         try eval([t3 ';']);
         catch
             EMHaved34S = 0; 
         end
     end
 end
end

% 42. only allow EMHaved34S to equal 1 if EMHaveSO4 is not equal to 0
if EMHaveSO4==0
   EMHaved34S = 0;
end

% 42. check if an appropriate pyrite end-member exists
if ~isempty(CoupleFeS2SO4intoEM)
    if   length(fieldnames(eval(sprintf('EM.%s.pyri',EMdatasource))))>0
    else disp(['catastrophic failure! CoupleFeS2SO4intoEM contains entries, but no ranges for normalized FeS2 SO4 are given in EM group ' EMdatasource '. check end-members spreadsheet']); 
         FailCase = 1;
    end
end

% 43. check if the choices of Z make sense
if CalculateRZCWY== 1
    if ZfromEM == 1        && sum(ismember(ObsList,'SO4'))<1;  disp(sprintf('catastrophic failure in simulation "%s"! Z is from End-members, but SO4 is not included in the river inversion.',                           ScenarioList{ScenarioListIndex})); FailCase = 1; end
    if ZfromSO4excess == 1 && sum(ismember(ObsList,'SO4'))>0;  disp(sprintf('catastrophic failure in simulation "%s"! Z is from SO4-excess, but SO4 is included in the river inversion. there will not be a SO4 excess', ScenarioList{ScenarioListIndex})); FailCase = 1; end
    if ZfromSO4excess == 1 && EMHaveSO4 == 0;                  disp(sprintf('catastrophic failure in simulation "%s"! calculation of SO4 excess is impossible because the end-members lack SO4 ratios',                  ScenarioList{ScenarioListIndex})); FailCase = 1; end
end

% 44. check for ResetDegasDICContribution, DegasDICContributionMin, and DegasDICContributionMax
if ~exist('ResetDegasDICContribution','var'); ResetDegasDICContribution = 0;  end
if ~exist('DegasDICContributionMin','var');   DegasDICContributionMin   = 0;  end
if ~exist('DegasDICContributionMax','var');   DegasDICContributionMax   = 0;  end


%% (4) display the run information
% format a string displaying the simulation information
if  FailCase == 0 
    str_ObsList = []; for ionlistindx=1:length(ObsList); str_ObsList = [str_ObsList  ObsList{ionlistindx} ', ']; end
    str_EMlist  = []; for ionlistindx=1:length(EMList0); str_EMlist  = [str_EMlist   EMList0{ionlistindx} ', ']; end
    str_ObsList  = str_ObsList(1:end-2); str_ObsList = [' ' str_ObsList];
    str_EMlist   = str_EMlist(1:end-2); 
    
    if    isequal(PrecProcessing,'ClCrit')
          if     ClCriticalValuesGiven==1; PrecProcessingText = 'with given constant ClCrit';    
          elseif ClCriticalValuesGiven==0  PrecProcessingText = 'with ClCrit values equal to sample Cl';
          end
    else                                   PrecProcessingText = 'as a standard end-member';
    end    
    if     AdjustRiverObs ==1; erroradjusttext = ' accounting for analytical error';
    elseif AdjustRiverObs ==0  erroradjusttext = ' not accounting for analytical error'; end

    if     isequal(IterateOver,'Samples');     txt1 = sprintf('running %se%s max simulations and %s max sucesses with datasource: %s',num2str(floor(maxiterations/(10^floor(log10(maxiterations))))),num2str(floor(log10(maxiterations))),num2str(maxsuccess),Riverdatasource);                     
    elseif isequal(IterateOver,'End-members'); txt1 = sprintf('running %s simulations with datasource: %s',num2str(round(numberiterations,1)),Riverdatasource); 
    end        
    txt2 = ['simulation variables are' str_ObsList ',' erroradjusttext ', ' sprintf('normalized to %s',NormalizationType)]; % keep all four % signs to actually print 1
    txt3 = sprintf('simulation end-members are %s from end-member set %s in units of %s',str_EMlist,EMdatasource,EMUnits);
    txt4 = sprintf('treating precipitation %s, and solving simulation with %s',PrecProcessingText,Solver);        
    
    normobs = [];
    for lll=1:length(ObsInNormalization)
        normobs = [normobs ObsInNormalization{lll} ', '];
    end
    normobs = normobs(1:end-2);
    txt5 = sprintf('observations included in the normalization are %s',normobs);
    
    printstring = [txt1  '\n'  txt2 '\n' txt5 '\n' txt3 '\n' txt4 '\n'];
    
    % if printing the string was requested, display the relevant scenario information
    if printtxt==1
       fprintf(printstring);
    end
end

% remove clutter from the workspace
clear EM; 
clear printstring; 
clear str_EMlist;
clear str_ObsList;
clear test;
clear testListChargeClosure;
clear testListNormClosure;
clear txt1;
clear txt2;
clear txt3;
clear txt4;
clear txt5;
