function [successcounter, smpres, EMSucessIter0Matrix, EMSucessIter0Matrix_frac, RiverColumnSucessMatrix, smpresmisfit, smpfunmisfit, badmbcounter, InvRivDatsave] = MEANDIR_EvaluateInversionInstance(i, X,nEM,nOL,EMIterInst0, EMIterInst0_updated, RiverColumn0,ObsList,ErrorCutMinMB,ErrorCutMaxMB,successcounter,InvRivDat,smpres,NormalizationType, badmbcounter, EMSucessIter0Matrix, EMSucessIter0Matrix_frac, RiverColumnSucessMatrix, smpresmisfit, smpfunmisfit, river, EM, EMList0, EMdatasource, InvRivDatsave, EMHaveSO4, MinFractionalContribution, MaxFractionalContribution, Solver, isotopeposition, carbonisotopematch, SolveCFList, abspos, relpos, WeightingList, ConvertDelta2RList, ConvertDelta2R, ImposeNormalizationCheck, EndMembersWithNegativeRatios, go, functioncost)

         % this function evaluates the results of a single inversion instance. first, the function evaluates whether the
         % fractional contributions to the normalization variables are within the user-defined ranges, as defined in 
         % MinFractionalContribution and MaxFractionalContribution. if so, the function then calculates the fractional
         % contribution of each end-member to each observation and ensures that the sum of these contributions is
         % within the user-defined ranges set in ErrorCutMinMB and ErrorCutMaxMB. in this case the instance is saved 
         % as a successful simulation result.
         
         % confirmation the the fractional contributions are within the specified ranges for the simulation 
         if     sum(X>=MinFractionalContribution & X<=MaxFractionalContribution)==nEM         
             
                % unpack the current sample results and temporary results data
                                   smpresultsd34S   = smpres.smpresultsd34S;   uniresultsd34S   = NaN(nEM,1);     % d34S
                                   smpresultsSO4    = smpres.smpresultsSO4;    uniresultsSO4    = NaN(nEM,1);     % SO4
                                   smpresultsnorm   = smpres.smpresultsnorm;   uniresultsnorm   = NaN(nEM,1);     % norm
                if go.ALK    == 1; smpresultsALK    = smpres.smpresultsALK;    uniresultsALK    = NaN(nEM,1); end % ALK
                if go.DIC    == 1; smpresultsDIC    = smpres.smpresultsDIC;    uniresultsDIC    = NaN(nEM,1); end % DIC
                if go.Ca     == 1; smpresultsCa     = smpres.smpresultsCa;     uniresultsCa     = NaN(nEM,1); end % Ca
                if go.Mg     == 1; smpresultsMg     = smpres.smpresultsMg;     uniresultsMg     = NaN(nEM,1); end % Mg
                if go.Na     == 1; smpresultsNa     = smpres.smpresultsNa;     uniresultsNa     = NaN(nEM,1); end % Na
                if go.K      == 1; smpresultsK      = smpres.smpresultsK;      uniresultsK      = NaN(nEM,1); end % K
                if go.Sr     == 1; smpresultsSr     = smpres.smpresultsSr;     uniresultsSr     = NaN(nEM,1); end % Sr
                if go.Fe     == 1; smpresultsFe     = smpres.smpresultsFe;     uniresultsFe     = NaN(nEM,1); end % Fe
                if go.Cl     == 1; smpresultsCl     = smpres.smpresultsCl;     uniresultsCl     = NaN(nEM,1); end % Cl
                if go.NO3    == 1; smpresultsNO3    = smpres.smpresultsNO3;    uniresultsNO3    = NaN(nEM,1); end % NO3
                if go.PO4    == 1; smpresultsPO4    = smpres.smpresultsPO4;    uniresultsPO4    = NaN(nEM,1); end % PO4
                if go.Si     == 1; smpresultsSi     = smpres.smpresultsSi;     uniresultsSi     = NaN(nEM,1); end % Si
                if go.Ge     == 1; smpresultsGe     = smpres.smpresultsGe;     uniresultsGe     = NaN(nEM,1); end % Ge
                if go.Li     == 1; smpresultsLi     = smpres.smpresultsLi;     uniresultsLi     = NaN(nEM,1); end % Li
                if go.F      == 1; smpresultsF      = smpres.smpresultsF;      uniresultsF      = NaN(nEM,1); end % F
                if go.B      == 1; smpresultsB      = smpres.smpresultsB;      uniresultsB      = NaN(nEM,1); end % B
                if go.Re     == 1; smpresultsRe     = smpres.smpresultsRe;     uniresultsRe     = NaN(nEM,1); end % Re
                if go.Mo     == 1; smpresultsMo     = smpres.smpresultsMo;     uniresultsMo     = NaN(nEM,1); end % Mo
                if go.Os     == 1; smpresultsOs     = smpres.smpresultsOs;     uniresultsOs     = NaN(nEM,1); end % Os
                if go.HCO3   == 1; smpresultsHCO3   = smpres.smpresultsHCO3;   uniresultsHCO3   = NaN(nEM,1); end % HCO3
                if go.d7Li   == 1; smpresultsd7Li   = smpres.smpresultsd7Li;   uniresultsd7Li   = NaN(nEM,1); end % d7Li
                if go.d13C   == 1; smpresultsd13C   = smpres.smpresultsd13C;   uniresultsd13C   = NaN(nEM,1); end % d13C
                if go.d18O   == 1; smpresultsd18O   = smpres.smpresultsd18O;   uniresultsd18O   = NaN(nEM,1); end % d18O
                if go.d26Mg  == 1; smpresultsd26Mg  = smpres.smpresultsd26Mg;  uniresultsd26Mg  = NaN(nEM,1); end % d26Mg
                if go.d30Si  == 1; smpresultsd30Si  = smpres.smpresultsd30Si;  uniresultsd30Si  = NaN(nEM,1); end % d30Si
                if go.d42Ca  == 1; smpresultsd42Ca  = smpres.smpresultsd42Ca;  uniresultsd42Ca  = NaN(nEM,1); end % d42Ca
                if go.d44Ca  == 1; smpresultsd44Ca  = smpres.smpresultsd44Ca;  uniresultsd44Ca  = NaN(nEM,1); end % d44Ca
                if go.d56Fe  == 1; smpresultsd56Fe  = smpres.smpresultsd56Fe;  uniresultsd56Fe  = NaN(nEM,1); end % d56Fe
                if go.Sr8786 == 1; smpresultsSr8786 = smpres.smpresultsSr8786; uniresultsSr8786 = NaN(nEM,1); end % Sr8786
                if go.d98Mo  == 1; smpresultsd98Mo  = smpres.smpresultsd98Mo;  uniresultsd98Mo  = NaN(nEM,1); end % d98Mo
                if go.Os8788 == 1; smpresultsOs8788 = smpres.smpresultsOs8788; uniresultsOs8788 = NaN(nEM,1); end % Os8788
                if go.Fmod   == 1; smpresultsFmod   = smpres.smpresultsFmod;   uniresultsFmod   = NaN(nEM,1); end % Fmod

                % before iterating over the relevant ion ratios, assign the normalization ion
                uniresultsnorm = X;

                mb_vector = [];  mb_count = 0;  % define an empty matrix for this sample
                for k = 1:nOL                   % k is the k'th ion ratio
                    workingresult = NaN(nEM,1); % empty the working results vector   
                    for j = 1:nEM               % j is the j'th end-member   
                        workingresult(j,1) = X(j).*EMIterInst0_updated(k,j)./RiverColumn0(k); 
                    end

                    % define the uniresults 
                    if     isequal(ObsList{k},'ALK');    uniresultsALK    = workingresult; % ALK
                    elseif isequal(ObsList{k},'DIC');    uniresultsDIC    = workingresult; % DIC
                    elseif isequal(ObsList{k},'Ca');     uniresultsCa     = workingresult; % Ca
                    elseif isequal(ObsList{k},'Mg');     uniresultsMg     = workingresult; % Mg
                    elseif isequal(ObsList{k},'Na');     uniresultsNa     = workingresult; % Na
                    elseif isequal(ObsList{k},'K');      uniresultsK      = workingresult; % K
                    elseif isequal(ObsList{k},'Sr');     uniresultsSr     = workingresult; % Sr
                    elseif isequal(ObsList{k},'Fe');     uniresultsFe     = workingresult; % Fe
                    elseif isequal(ObsList{k},'Cl');     uniresultsCl     = workingresult; % Cl
                    elseif isequal(ObsList{k},'SO4');    uniresultsSO4    = workingresult; % SO4
                    elseif isequal(ObsList{k},'NO3');    uniresultsNO3    = workingresult; % NO3
                    elseif isequal(ObsList{k},'PO4');    uniresultsPO4    = workingresult; % PO4
                    elseif isequal(ObsList{k},'Si');     uniresultsSi     = workingresult; % Si
                    elseif isequal(ObsList{k},'Ge');     uniresultsGe     = workingresult; % Ge
                    elseif isequal(ObsList{k},'Li');     uniresultsLi     = workingresult; % Li
                    elseif isequal(ObsList{k},'F');      uniresultsF      = workingresult; % F
                    elseif isequal(ObsList{k},'B');      uniresultsB      = workingresult; % B
                    elseif isequal(ObsList{k},'Re');     uniresultsRe     = workingresult; % Re
                    elseif isequal(ObsList{k},'Mo');     uniresultsMo     = workingresult; % Mo
                    elseif isequal(ObsList{k},'Os');     uniresultsOs     = workingresult; % Os
                    elseif isequal(ObsList{k},'HCO3');   uniresultsHCO3   = workingresult; % HCO3
                    elseif isequal(ObsList{k},'d7Li');   uniresultsd7Li   = workingresult; % d7Li
                    elseif isequal(ObsList{k},'d13C');   uniresultsd13C   = workingresult; % d13C
                    elseif isequal(ObsList{k},'d18O');   uniresultsd18O   = workingresult; % d18O
                    elseif isequal(ObsList{k},'d26Mg');  uniresultsd26Mg  = workingresult; % d26Mg
                    elseif isequal(ObsList{k},'d30Si');  uniresultsd30Si  = workingresult; % d30Si
                    elseif isequal(ObsList{k},'d34S');   uniresultsd34S   = workingresult; % d34S
                    elseif isequal(ObsList{k},'d42Ca');  uniresultsd42Ca  = workingresult; % d42Ca
                    elseif isequal(ObsList{k},'d44Ca');  uniresultsd44Ca  = workingresult; % d44Ca
                    elseif isequal(ObsList{k},'d56Fe');  uniresultsd56Fe  = workingresult; % d56Fe
                    elseif isequal(ObsList{k},'Sr8786'); uniresultsSr8786 = workingresult; % Sr8786
                    elseif isequal(ObsList{k},'d98Mo');  uniresultsd98Mo  = workingresult; % d98Mo
                    elseif isequal(ObsList{k},'Os8788'); uniresultsOs8788 = workingresult; % Os8788
                    elseif isequal(ObsList{k},'Fmod');   uniresultsFmod   = workingresult; % Fmod
                    end   
                    mb_count              = mb_count+1;
                    mb_vector(mb_count,:) = sum(workingresult); % store the data 
                end

                 % make copies of isotopeposition, ErrorCutMinMB, andErrorCutMaxMB
                 ErrorCutMinMB_mbcheck   = ErrorCutMinMB;
                 ErrorCutMaxMB_mbcheck   = ErrorCutMaxMB;
                 isotopeposition_mbcheck = isotopeposition;

                 % determine if you need to save the normalization
                 if  isequal(NormalizationType,'SumObs') & ImposeNormalizationCheck==1
                     uniresultsSumObs                    = X;                                    % define the uniresultsSumObs
                     mb_count                            = mb_count+1;                           % add 1 to the mb count
                     mb_vector(mb_count,:)               = sum(uniresultsSumObs);                % store the sum
                     ErrorCutMinMB_mbcheck(1,mb_count)   = min(ErrorCutMinMB(~isotopeposition)); % get the minimum at a non-isotope position
                     ErrorCutMaxMB_mbcheck(1,mb_count)   = max(ErrorCutMaxMB(~isotopeposition)); % get the maximum at a non-isotope position
                     isotopeposition_mbcheck(mb_count,1) = 0;
                 end                     

                 % only keep data within specified error window, and do not yet consider whether isotopes are faithfully reconstructed
                 if sum(mb_vector(~isotopeposition_mbcheck) >= ErrorCutMinMB_mbcheck(~isotopeposition_mbcheck)'/100  & mb_vector(~isotopeposition_mbcheck) <= ErrorCutMaxMB_mbcheck(~isotopeposition_mbcheck)'/100) == mb_count-sum(isotopeposition_mbcheck)

                     % check whether the isotopic ratios are within the specified bounds for the inversion to count as successful
                     isotopepasscriteria = 1; 
                     if sum(isotopeposition)>0  
                         isotopenames = ObsList(isotopeposition);
                         isotopeindx  = find(isotopeposition);
                         reconstructed_river = EMIterInst0_updated*uniresultsnorm;
                         for ll = 1:length(isotopeindx)
                             activeisoname = isotopenames{ll};
                             if isequal(activeisoname,'d7Li')  ;                                      reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'Li'));   river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'Li'));   end % d7Li
                             if isequal(activeisoname,'d18O')  ;                                      reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'SO4'));  river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'SO4'));  end % d18O
                             if isequal(activeisoname,'d26Mg') ;                                      reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'Mg'));   river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'Mg'));   end % d26Mg
                             if isequal(activeisoname,'d30Si') ;                                      reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'Si'));   river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'Si'));   end % d30Si
                             if isequal(activeisoname,'d34S')  ;                                      reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'SO4'));  river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'SO4'));  end % d34S
                             if isequal(activeisoname,'d42Ca') ;                                      reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'Ca'));   river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'Ca'));   end % d42Ca 
                             if isequal(activeisoname,'d44Ca') ;                                      reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'Ca'));   river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'Ca'));   end % d44Ca
                             if isequal(activeisoname,'d56Fe') ;                                      reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'Fe'));   river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'Fe'));   end % d56Fe
                             if isequal(activeisoname,'Sr8786');                                      reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'Sr'));   river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'Sr'));   end % Sr8786
                             if isequal(activeisoname,'d98Mo') ;                                      reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'Mo'));   river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'Mo'));   end % d98Mo
                             if isequal(activeisoname,'Os8788') ;                                     reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'Os'));   river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'Os'));   end % Os8788
                             if isequal(activeisoname,'d13C')  & isequal(carbonisotopematch,'HCO3');  reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'HCO3')); river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'HCO3')); end % d13C for HCO3
                             if isequal(activeisoname,'Fmod')  & isequal(carbonisotopematch,'HCO3');  reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'HCO3')); river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'HCO3')); end % Fmod for HCO3
                             if isequal(activeisoname,'d13C')  & isequal(carbonisotopematch,'DIC');   reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'DIC'));  river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'DIC'));  end % d13C for DIC
                             if isequal(activeisoname,'Fmod')  & isequal(carbonisotopematch,'DIC');   reconstructed_iso = reconstructed_river(isotopeindx(ll))./reconstructed_river(ismember(ObsList,'DIC'));  river_iso = RiverColumn0(isotopeindx(ll))./RiverColumn0(ismember(ObsList,'DIC'));  end % Fmod for DIC

                             Minisotopetargetdifference = ErrorCutMinMB(isotopeindx(ll)); % this value should be negative
                             Maxisotopetargetdifference = ErrorCutMaxMB(isotopeindx(ll)); % this value should be positive

                             if   sum(ismember(ConvertDelta2RList,ObsList{isotopeindx(ll)}))>0
                                  conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{isotopeindx(ll)}));
                                  Minisotopetargetdifference = Minisotopetargetdifference/1000*conv;
                                  Maxisotopetargetdifference = Maxisotopetargetdifference/1000*conv;
                             end                                         
                             if   ~((reconstructed_iso >= river_iso+Minisotopetargetdifference) & (reconstructed_iso <= river_iso+Maxisotopetargetdifference))
                                  isotopepasscriteria = 0;
                             end
                         end
                     end

                     if isotopepasscriteria==1

                        % if the code reaches this point, all of the mass balance criteria for dissolved abundances and all of the precision
                        % criteria for re-creating isotopic information have been satisfied. MEANDIR will now calculate the fractional
                        % contribution of each end-member to each variable and save the results of this simulation

                        % to check the result: [EMIterInst0_updated*X RiverColumn0]
                        successcounter = successcounter+1;

                        % make a SO4 results variable if SO4 was not in the simulation and the end-members have SO4 values.
                        % note that this means each sucessful inversion is associated with one pull from the normalized SO4
                        % distributions (as if FeS2/other were an independent end-member). however, it is possible to imagine an
                        % alternative calculation in which each simulation should be associated with a distribution of values, each of
                        % which could then be assoicated with its own distribution of d34S values.
                        if sum(ismember(ObsList,'SO4'))==0  & ~isequal(NormalizationType,'SO4') & EMHaveSO4 == 1
                            RiverSO4Data = river.model_variable.SO4./InvRivDat.norm; 
                            for j = 1:nEM        
                                EMvalue          = MEANDIR_FindEndMemberSO4Ratio(EM,EMList0,EMdatasource,NormalizationType,j, EndMembersWithNegativeRatios);
                                uniresultsSO4(j) = uniresultsnorm(j).*EMvalue./RiverSO4Data(i); 
                            end
                        end

                        % store the results   
                                           smpresultsnorm(:,successcounter)   = uniresultsnorm;   smpres.smpresultsnorm   = smpresultsnorm;   InvRivDatsave.norm(successcounter,1)   = InvRivDat.norm;       % normalization
                                           smpresultsSO4(:,successcounter)    = uniresultsSO4;    smpres.smpresultsSO4    = smpresultsSO4;    InvRivDatsave.SO4(successcounter,1)    = InvRivDat.SO4;        % SO4
                                           smpresultsd34S(:,successcounter)   = uniresultsd34S;   smpres.smpresultsd34S   = smpresultsd34S;   InvRivDatsave.d34S(successcounter,1)   = InvRivDat.d34S;       % d34S
                        if go.ALK    == 1; smpresultsALK(:,successcounter)    = uniresultsALK;    smpres.smpresultsALK    = smpresultsALK;    InvRivDatsave.ALK(successcounter,1)    = InvRivDat.ALK;    end % ALK
                        if go.DIC    == 1; smpresultsDIC(:,successcounter)    = uniresultsDIC;    smpres.smpresultsDIC    = smpresultsDIC;    InvRivDatsave.DIC(successcounter,1)    = InvRivDat.DIC;    end % DIC
                        if go.Ca     == 1; smpresultsCa(:,successcounter)     = uniresultsCa;     smpres.smpresultsCa     = smpresultsCa;     InvRivDatsave.Ca(successcounter,1)     = InvRivDat.Ca;     end % Ca
                        if go.Mg     == 1; smpresultsMg(:,successcounter)     = uniresultsMg;     smpres.smpresultsMg     = smpresultsMg;     InvRivDatsave.Mg(successcounter,1)     = InvRivDat.Mg;     end % Mg
                        if go.Na     == 1; smpresultsNa(:,successcounter)     = uniresultsNa;     smpres.smpresultsNa     = smpresultsNa;     InvRivDatsave.Na(successcounter,1)     = InvRivDat.Na;     end % Na
                        if go.K      == 1; smpresultsK(:,successcounter)      = uniresultsK;      smpres.smpresultsK      = smpresultsK;      InvRivDatsave.K(successcounter,1)      = InvRivDat.K;      end % K
                        if go.Sr     == 1; smpresultsSr(:,successcounter)     = uniresultsSr;     smpres.smpresultsSr     = smpresultsSr;     InvRivDatsave.Sr(successcounter,1)     = InvRivDat.Sr;     end % Sr
                        if go.Fe     == 1; smpresultsFe(:,successcounter)     = uniresultsFe;     smpres.smpresultsFe     = smpresultsFe;     InvRivDatsave.Fe(successcounter,1)     = InvRivDat.Fe;     end % Fe
                        if go.Cl     == 1; smpresultsCl(:,successcounter)     = uniresultsCl;     smpres.smpresultsCl     = smpresultsCl;     InvRivDatsave.Cl(successcounter,1)     = InvRivDat.Cl;     end % Cl
                        if go.NO3    == 1; smpresultsNO3(:,successcounter)    = uniresultsNO3;    smpres.smpresultsNO3    = smpresultsNO3;    InvRivDatsave.NO3(successcounter,1)    = InvRivDat.NO3;    end % NO3
                        if go.PO4    == 1; smpresultsPO4(:,successcounter)    = uniresultsPO4;    smpres.smpresultsPO4    = smpresultsPO4;    InvRivDatsave.PO4(successcounter,1)    = InvRivDat.PO4;    end % PO4
                        if go.Si     == 1; smpresultsSi(:,successcounter)     = uniresultsSi;     smpres.smpresultsSi     = smpresultsSi;     InvRivDatsave.Si(successcounter,1)     = InvRivDat.Si;     end % Si
                        if go.Ge     == 1; smpresultsGe(:,successcounter)     = uniresultsGe;     smpres.smpresultsGe     = smpresultsGe;     InvRivDatsave.Ge(successcounter,1)     = InvRivDat.Ge;     end % Ge
                        if go.Li     == 1; smpresultsLi(:,successcounter)     = uniresultsLi;     smpres.smpresultsLi     = smpresultsLi;     InvRivDatsave.Li(successcounter,1)     = InvRivDat.Li;     end % Li
                        if go.F      == 1; smpresultsF(:,successcounter)      = uniresultsF;      smpres.smpresultsF      = smpresultsF;      InvRivDatsave.F(successcounter,1)      = InvRivDat.F;      end % F
                        if go.B      == 1; smpresultsB(:,successcounter)      = uniresultsB;      smpres.smpresultsB      = smpresultsB;      InvRivDatsave.B(successcounter,1)      = InvRivDat.B;      end % B
                        if go.Re     == 1; smpresultsRe(:,successcounter)     = uniresultsRe;     smpres.smpresultsRe     = smpresultsRe;     InvRivDatsave.Re(successcounter,1)     = InvRivDat.Re;     end % Re
                        if go.Mo     == 1; smpresultsMo(:,successcounter)     = uniresultsMo;     smpres.smpresultsMo     = smpresultsMo;     InvRivDatsave.Mo(successcounter,1)     = InvRivDat.Mo;     end % Mo
                        if go.Os     == 1; smpresultsOs(:,successcounter)     = uniresultsOs;     smpres.smpresultsOs     = smpresultsOs;     InvRivDatsave.Os(successcounter,1)     = InvRivDat.Os;     end % Os
                        if go.HCO3   == 1; smpresultsHCO3(:,successcounter)   = uniresultsHCO3;   smpres.smpresultsHCO3   = smpresultsHCO3;   InvRivDatsave.HCO3(successcounter,1)   = InvRivDat.HCO3;   end % HCO3
                        if go.d7Li   == 1; smpresultsd7Li(:,successcounter)   = uniresultsd7Li;   smpres.smpresultsd7Li   = smpresultsd7Li;   InvRivDatsave.d7Li(successcounter,1)   = InvRivDat.d7Li;   end % d7Li
                        if go.d13C   == 1; smpresultsd13C(:,successcounter)   = uniresultsd13C;   smpres.smpresultsd13C   = smpresultsd13C;   InvRivDatsave.d13C(successcounter,1)   = InvRivDat.d13C;   end % d13C
                        if go.d18O   == 1; smpresultsd18O(:,successcounter)   = uniresultsd18O;   smpres.smpresultsd18O   = smpresultsd18O;   InvRivDatsave.d18O(successcounter,1)   = InvRivDat.d18O;   end % d18O
                        if go.d26Mg  == 1; smpresultsd26Mg(:,successcounter)  = uniresultsd26Mg;  smpres.smpresultsd26Mg  = smpresultsd26Mg;  InvRivDatsave.d26Mg(successcounter,1)  = InvRivDat.d26Mg;  end % d26Mg
                        if go.d30Si  == 1; smpresultsd30Si(:,successcounter)  = uniresultsd30Si;  smpres.smpresultsd30Si  = smpresultsd30Si;  InvRivDatsave.d30Si(successcounter,1)  = InvRivDat.d30Si;  end % d30Si
                        if go.d42Ca  == 1; smpresultsd42Ca(:,successcounter)  = uniresultsd42Ca;  smpres.smpresultsd42Ca  = smpresultsd42Ca;  InvRivDatsave.d42Ca(successcounter,1)  = InvRivDat.d42Ca;  end % d42Ca
                        if go.d44Ca  == 1; smpresultsd44Ca(:,successcounter)  = uniresultsd44Ca;  smpres.smpresultsd44Ca  = smpresultsd44Ca;  InvRivDatsave.d44Ca(successcounter,1)  = InvRivDat.d44Ca;  end % d44Ca
                        if go.d56Fe  == 1; smpresultsd56Fe(:,successcounter)  = uniresultsd56Fe;  smpres.smpresultsd56Fe  = smpresultsd56Fe;  InvRivDatsave.d56Fe(successcounter,1)  = InvRivDat.d56Fe;  end % d56Fe
                        if go.Sr8786 == 1; smpresultsSr8786(:,successcounter) = uniresultsSr8786; smpres.smpresultsSr8786 = smpresultsSr8786; InvRivDatsave.Sr8786(successcounter,1) = InvRivDat.Sr8786; end % Sr8786
                        if go.d98Mo  == 1; smpresultsd98Mo(:,successcounter)  = uniresultsd98Mo;  smpres.smpresultsd98Mo  = smpresultsd98Mo;  InvRivDatsave.d98Mo(successcounter,1)  = InvRivDat.d98Mo;  end % d98Mo
                        if go.Os8788 == 1; smpresultsOs8788(:,successcounter) = uniresultsOs8788; smpres.smpresultsOs8788 = smpresultsOs8788; InvRivDatsave.Os8788(successcounter,1) = InvRivDat.Os8788; end % Os8788
                        if go.Fmod   == 1; smpresultsFmod(:,successcounter)   = uniresultsFmod;   smpres.smpresultsFmod   = smpresultsFmod;   InvRivDatsave.Fmod(successcounter,1)   = InvRivDat.Fmod;   end % Fmod

                        % this is an important step: the matrix with fractionation factors (EMIterInst0) is being saved into
                        % "EMSucessIter0Matrix_frac", while the matrix with substitued values for the isotopic composition of
                        % fractionated end-members (EMIterInst0_updated) is saved into "EMSucessIter0Matrix". for an inversion
                        % without fractionation, these two matricies should be identical
                        EMSucessIter0Matrix(:,:,successcounter)      = EMIterInst0_updated; % record the end-member matrix, updated for fractionation effects
                        EMSucessIter0Matrix_frac(:,:,successcounter) = EMIterInst0;         % record the end-member matrix, still with fractionation factors
                        RiverColumnSucessMatrix(successcounter,:)    = RiverColumn0;        % save the chemistry of the river water sample

                        % calculate the misfit vector for model results and observations
                        if     isequal(Solver,'mldivide')  | isequal(Solver,'lsqnonneg')
                               misfitvec = (RiverColumn0 - EMIterInst0_updated*X).^2; 
                               smpresmisfit(1,successcounter) = sqrt(sum(misfitvec(SolveCFList)));
                        elseif isequal(Solver,'optimize') | isequal(Solver,'mldivide_optimize') | isequal(Solver,'lsqnonneg_optimize')
                               misfitvec_rel = ((RiverColumn0 - EMIterInst0_updated*X)./RiverColumn0).^2;
                               misfitvec_abs = (RiverColumn0  - EMIterInst0_updated*X).^2;
                               smpresmisfit(1,successcounter) = sqrt(sum([WeightingList(SolveCFList & relpos).*misfitvec_rel(SolveCFList & relpos); WeightingList(SolveCFList & abspos).*misfitvec_abs(SolveCFList & abspos)]));
                        end

                        % save the result of the cost function evaluation
                        smpfunmisfit(1,successcounter) = functioncost; 
                        % without ClCritical and where are observations are given to the cost 
                        % function, smpfunmisfit and misfitvec should have the same values.

                    else % if isotopepasscritera is not equal to 1
                         badmbcounter = badmbcounter + 1;
                    end
             else  % if the mass balance check did not work
                  badmbcounter = badmbcounter + 1;
             end  % end mb check
        end       % end fractional contribution check
end               % end function