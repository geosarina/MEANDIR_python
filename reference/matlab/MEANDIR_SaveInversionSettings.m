    % this script saves some of the information associated with the inversion
    
    % (1) save scenario information
    river.info.successfullsimulations                 = successfullsimulations;              clear successfullsimulations;
    river.info.badmbcount                             = badmbcount;                          clear badmbcount;
    river.info.iterationcounts                        = iterationcounts;                     clear iterationcounts;
    river.info.s                                      = s;                                   clear s;
    
    % (2) save the river observations and end-member information
    river.RiverMatrix                                 = RiverMatrix0;                        clear RiverMatrix0;
    river.EndMembers.EM                               = eval(sprintf('EM.%s',EMdatasource)); clear EM;
    river.EndMembers.distributions.truncated          = TrEMdist;                            clear TrEMdist;
    river.EndMembers.distributions.original           = OrEMdist;                            clear OrEMdist;
        
    % (3) save inversion settings
    river.settings.Solver                             = Solver;                              clear Solver;
    river.settings.IterateOver                        = IterateOver;                         clear IterateOver;
    river.settings.maxiterations                      = maxiterations;                       clear maxiterations;
    river.settings.maxsuccess                         = maxsuccess;                          clear maxsuccess;
    river.settings.maxzerohits                        = maxzerohits;                         clear maxzerohits;
    river.settings.CostFunType                        = CostFunType;                         clear CostFunType;
    river.settings.WeightingList                      = WeightingList;                       clear WeightingList;
    river.settings.ErrorCutMinMB                      = ErrorCutMinMB;                       clear ErrorCutMinMB;
    river.settings.ErrorCutMaxMB                      = ErrorCutMaxMB;                       clear ErrorCutMaxMB;
    river.settings.numberiterations                   = numberiterations;                    clear numberiterations;
    river.settings.MisfitCuts                         = MisfitCuts;                          clear MisfitCuts;
    river.settings.CullOn                             = CullOn;                              clear CullOn;
    river.settings.saveuncutdata                      = saveuncutdata;                       clear saveuncutdata;
    river.settings.Riverdatasource                    = Riverdatasource;                     clear Riverdatasource;
    river.settings.AdjustRiverObs                     = AdjustRiverObs;                      clear AdjustRiverObs;
    river.settings.ObsList                            = ObsList;                             clear ObsList;
    river.settings.nCFList                            = nCFList;                             clear nCFList;
    river.settings.NormalizationType                  = NormalizationType;                   clear NormalizationType;
    river.settings.CationsInNormalization             = CationsInNormalization;              clear CationsInNormalization;
    river.settings.AnionsInNormalization              = AnionsInNormalization;               clear AnionsInNormalization;
    river.settings.NeutralInNormalization             = NeutralInNormalization;              clear NeutralInNormalization;
    river.settings.ObsInNormalization                 = ObsInNormalization;                  clear ObsInNormalization;
    river.settings.AllIonsExplicitlyResolved          = AllIonsExplicitlyResolved;           clear AllIonsExplicitlyResolved;
    river.settings.IonCharges                         = IonCharges;                          clear IonCharges;
    river.settings.EMdatasource                       = EMdatasource;                        clear EMdatasource;
    river.settings.EMUnits                            = EMUnits;                             clear EMUnits;
    river.settings.EMList0                            = EMList0;                             clear EMList0;
    river.settings.ListNormClosure                    = ListNormClosure;                     clear ListNormClosure;
    river.settings.ListChargeClosure                  = ListChargeClosure;                   clear ListChargeClosure;
    river.settings.MinFractionalContribution          = MinFractionalContribution;           clear MinFractionalContribution;
    river.settings.MaxFractionalContribution          = MaxFractionalContribution;           clear MaxFractionalContribution;
    river.settings.CoupleFeS2SO4intoEM                = CoupleFeS2SO4intoEM;                 clear CoupleFeS2SO4intoEM;
    river.settings.CoupleFeS2d34SintoEM               = CoupleFeS2d34SintoEM;                clear CoupleFeS2d34SintoEM;
    river.settings.RecordFullFeS2Distribution         = RecordFullFeS2Distribution;          clear RecordFullFeS2Distribution;
    river.settings.BalanceEvaporite                   = BalanceEvaporite;                    clear BalanceEvaporite;
    river.settings.EndMembersWithNegativeRatios       = EndMembersWithNegativeRatios;        clear EndMembersWithNegativeRatios;
    river.settings.PrecProcessing                     = PrecProcessing;                      clear PrecProcessing;
    river.settings.ClCriticalValuesGiven              = ClCriticalValuesGiven;               clear ClCriticalValuesGiven;
    river.settings.CalculateRZCWY                     = CalculateRZCWY;                      clear CalculateRZCWY;
    river.settings.R_Numerator_EMList                 = R_Numerator_EMList;                  clear R_Numerator_EMList;
    river.settings.R_Numerator_IonList                = R_Numerator_IonList;                 clear R_Numerator_IonList;
    river.settings.RZC_Denominator_EMList             = RZC_Denominator_EMList;              clear RZC_Denominator_EMList;
    river.settings.RZC_Denominator_IonList            = RZC_Denominator_IonList;             clear RZC_Denominator_IonList;
    river.settings.Z_NumeratorType                    = Z_NumeratorType;                     clear Z_NumeratorType;
    river.settings.Z_Numerator_EMList                 = Z_Numerator_EMList;                  clear Z_Numerator_EMList;
    river.settings.C_Numerator_EMList                 = C_Numerator_EMList;                  clear C_Numerator_EMList;
    river.settings.EMsources                          = EMsources;                           clear EMsources;
    river.settings.EMsinks                            = EMsinks;                             clear EMsinks;
    river.settings.indxisopos                         = indxisopos;                          clear indxisopos;
    river.settings.isotopeposition                    = isotopeposition;                     clear isotopeposition;
    river.settings.ZfromSO4excess                     = ZfromSO4excess;                      clear ZfromSO4excess;
    river.settings.ZfromriverSO4                      = ZfromriverSO4;                       clear ZfromriverSO4;
    river.settings.ZfromEM                            = ZfromEM;                             clear ZfromEM;
    river.settings.Znotcalculated                     = Znotcalculated;                      clear Znotcalculated;
    river.settings.EMHaved34S                         = EMHaved34S;                          clear EMHaved34S;
    river.settings.remaininginstances                 = remaininginstances;                  clear remaininginstances;
    river.settings.ScenarioList                       = ScenarioList;
    river.settings.ScenarioListIndex                  = ScenarioListIndex;    
    river.settings.carbonisotopematch                 = carbonisotopematch;                  clear carbonisotopematch; 
    