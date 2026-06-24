function [RiverColumn, EMIterInst, EMIterInst0, Xdirect, EMList, distEMdist, iterationcounter] = MEANDIR_ClCriticalCorrection(river, InvRivDat, EMIterInst0, Xdirect, IterateOver, ClCriticalValuesGiven, nOL, EMList0, EM, EMdatasource, EndMembersWithNegativeRatios, ObsList, IonCharges, isotopeposition,indxisopos, i, iterationcounter, maxiterations, RiverColumn0, NormalizationType, ObsInNormalization, CationsInNormalization, AnionsInNormalization, NeutralInNormalization, ListNormClosure, ListChargeClosure,  CoupleFeS2d34SintoEM, AllIonsExplicitlyResolved, BalanceEvaporite, carbonisotopematch, TrEMdist, distEMdist0, distEMdist0min, distEMdist0max, smpadddist_mean, smpadddist_std, smpadddist_min, smpadddist_max, ConvertDelta2RList, ConvertDelta2R, lgu_sign, MinFractionalContribution)

        % this function applies ClCritical values, given by the user, to river observations. the function outputs: 
        % RiverColumn      : a set of river ratios adjusted for ClCritical values
        % EMIterInst       : a set of end-member ratios adjusted for ClCritical values (no precipitation end-member)
        % EMIterInst0      : the full end-member matrix (can be reset during this script)
        % Xdirect          : a vector that will ultimately hold fractional contributions from end-members, with a non-NaN value for precipitation
        % EMList           : an updated vector of end-members without precipitation
        % distEMdist       : a matrix of end-member distributions (no precipitation end-member)
        % iterationcounter : if ClCritical is NaN, iterationcounter if re-set to maxiterations so that the model moves to the next sample

        RiverColumn          = RiverColumn0;                                     % define RiverColumn, which will have information on updated chemistry of this sample
        ClOverNormRiverRatio = RiverColumn0(ismember(ObsList,'Cl'),1);           % recall that the data here are already updated for analytical error and normalized
        ClRiverConcEqui      = InvRivDat.Cl;                                     % error-adjusted Cl concentration in the river water sample
        precEM_indx          = ismember(EMList0,'prec');                         % the index position of precipitation within the list of end-members
        Cl_indx              = ismember(ObsList,'Cl');                           % the index position of Cl within the list of variables
        maxredefinition      = 100;                                              % the maximum number of times to re-define the end-member matrix   

        ClCrit = NaN;                                                            % define an initial NaN value for ClCritical, to then be overwritten
        if    ClCriticalValuesGiven == 1; ClCrit = river.model_variable.ClCr(i); % if ClCritical values are given in the simulation parameters, use them
        else;                             ClCrit = ClRiverConcEqui;              % if the precipitation should be processed with ClCrit values but none are given (ClCriticalValuesGiven=0), assume 100% of Cl- comes from precipitation
        end

        % (1) check if removal of the specific amount of Cl will result in negative ratios for the other dissolved constituents. 
        % if so, re-pull the precipitation chemistry
        if    ~isnan(ClCrit)
               % first, test if the current end-member will return negative values for other dissolved variables when Cl is removed later in the script.
               % if so, re-define the end-member matrix up to the number of times specified by maxreiterations.
               precratios = MEANDIR_QuantifyPrecRatios(ObsList, EMIterInst0, InvRivDat, ClCrit, precEM_indx, Cl_indx);
               if sum(precratios>1)>0 & isequal(IterateOver,'Samples') & sum(MinFractionalContribution<0)==0 % if IterateOver is EndMembers, keep the current set of end-member values
                  stillimpossible        = 1;
                  stillimpossiblecounter = 0;
                  while stillimpossible == 1
                      % if the current matrix of end-member chemistry, for the given ClCritical value, would produce a 
                      % negative concentration for one of the non-isotope variables in ObsList, make maxredefinition
                      % attempts to find an EM matrix that produces physicallly realizable results 
                      stillimpossiblecounter  = stillimpossiblecounter + 1;
                      [EMIterInst0, FailCase] = MEANDIR_PullEndMemberRatios(EMList0, EM, EMdatasource, EndMembersWithNegativeRatios, ObsList, IonCharges, isotopeposition,indxisopos, i, iterationcounter, maxiterations, ...
                                                RiverColumn0, NormalizationType, ObsInNormalization, CationsInNormalization, AnionsInNormalization, NeutralInNormalization, ...
                                                ListNormClosure, ListChargeClosure,  CoupleFeS2d34SintoEM, AllIonsExplicitlyResolved, BalanceEvaporite, ...
                                                carbonisotopematch, TrEMdist, distEMdist0, distEMdist0min, distEMdist0max, smpadddist_mean, smpadddist_std, smpadddist_min, smpadddist_max, ConvertDelta2RList, ConvertDelta2R, lgu_sign);
                      if FailCase == 1
                         disp('warning! problem finding end-members in MEANDIR_ClCriticalCorrection. program ended'); 
                         return;
                      end
                      precratios = MEANDIR_QuantifyPrecRatios(ObsList, EMIterInst0, InvRivDat, ClCrit, precEM_indx, Cl_indx);
                      if sum(precratios>1)==0
                         stillimpossible  = 0;
                      end
                      if stillimpossiblecounter>maxredefinition
                         RiverColumn      = NaN(size(RiverColumn)); 
                         stillimpossible  = 0;
                      end
                  end
               end  
        elseif isnan(ClCrit) 
               % if the ClCrit value is NaN, just short-circuit the looping so we don't have to wait for maxiterations number of failed solutions
               iterationcounter =  maxiterations;
        end

         % (2) now pull the Cl-normalized value of the acceptable precipitation term
         ClOverNormPrecRatio = EMIterInst0(Cl_indx,precEM_indx);

         % calculate fractional contribution of precipitation to normalization ion
         if     ClRiverConcEqui <= ClCrit
                % when there is less Cl in the river than the ClCrit value 100% of the Cl comes from precipitation
                f_norm_prec = ClOverNormRiverRatio/ClOverNormPrecRatio; 

         elseif ClRiverConcEqui > ClCrit
                % when there is more Cl in the river than the ClCrit value, 100% of the ClCrit comes
                % from precipitation, and the remainder has to be taken up by something else. 
                f_norm_prec = (ClCrit/ClRiverConcEqui)*ClOverNormRiverRatio/ClOverNormPrecRatio;

         elseif isnan(ClRiverConcEqui) || isnan(ClCrit)
                % if the ClCrit is NaN, then call the fractional contribution also NaN
                f_norm_prec = NaN;
         end

         % (3) update the river data for the removal of ClCrit
         for j=1:nOL                             
               RiverColumn(j,1) = [RiverColumn0(j,1) - EMIterInst0(j,precEM_indx)*f_norm_prec];
         end

         % the following code sets river entries arbitrarily close to zero to equal zero.
         if sum(RiverColumn<5e-16 & RiverColumn>-5e-16) > 0
               RiverColumn(RiverColumn<5e-16 & RiverColumn>-5e-16) = 0;        
         end

         % remove the precipitation variable from the EMIterInst0
         EMIterInst  = [];                          % set EMIterInst to be empty
         EMIterInst  = EMIterInst0(:,~precEM_indx); % remove precipitation from the EMIter matrix  
         distEMdist  = distEMdist0(:,~precEM_indx); % remove precipitation from the distEMdist matrix 

         % (4) now update Xdirect and EMlist
         Xdirect(ismember(EMList0,'prec')) = f_norm_prec;   
         EMList = EMList0(~ismember(EMList0,'prec'));     

end      % end of function