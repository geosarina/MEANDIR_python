function [EMIterInst, FailCase] = MEANDIR_PullEndMemberRatios(EMList0, EM, EMdatasource, EndMembersWithNegativeRatios, ObsList, IonCharges, isotopeposition,indxisopos, i, iterationcounter, maxiterations, RiverColumn0, NormalizationType, ObsInNormalization, CationsInNormalization, AnionsInNormalization, NeutralInNormalization, ListNormClosure, ListChargeClosure, CoupleFeS2d34SintoEM, AllIonsExplicitlyResolved, BalanceEvaporite, carbonisotopematch, TrEMdist, distEMdist0, distEMdist0_min, distEMdist0_max, smpadddist_mean, smpadddist_std, smpadddist_min, smpadddist_max, ConvertDelta2RList, ConvertDelta2R, lgu_sign)

    % this function produces an internally consistent set of end-members
    
    [nIR, nEM]           = size(TrEMdist); % use the cell containing truncated distributions to find the number of end-members and observations
    EMIterInst           = NaN(nIR,nEM);   % "EMIterInst" is one possible set of end-member chemical compositions
    EMsetisnotacceptable = 1;              % value to indicate if a set of acceptable end-member compositions has been found
    FailCase             = 0;              % value to determine if end-members were found successfully
    trycount             = 0;              % counter for the number of attempts to generate an internally consistent end-member
    maxtrycount          = 500;            % the maximum number of attempts to generate an internally consistent end-member   
            
    while EMsetisnotacceptable == 1        % run the while loop until a set of acceptable end-member chemistry has been found
          trycount   = trycount+1;         % add 1 to the counter tracking the number of attempts to find a reasonable set of end-member chemistries
          listallNaN = 1;                  % if there is a problem with the end-members, list the NaN positions, unless all values are NaN
            
          % for each observation and end-member, pull a random value from the specified distribution (saved in TrEMdist)
          for ii=1:nIR % iterate over observations
          for jj=1:nEM % iterate over end-members
            
              % (1) first consider a user-defined uniform or normal distribution with absolute ranges. this is the 
              % code called in the vast majority of cases and it simply pulls a value from the pre-defined distributions
              if isequal(distEMdist0{ii,jj},'UNI') | isequal(distEMdist0{ii,jj},'NOR') | isequal(distEMdist0{ii,jj},'EPS-UNI') | isequal(distEMdist0{ii,jj},'EPS-NOR')
                  EMIterInst(ii,jj) = random(TrEMdist{ii,jj},1,1);
              end            
              
              % for a user-defined log-uniform distribution with absolute ranges
              % recall that the corresponding entry of TrEMdist is a uniform distribution 
              % bounded between the natural log of the end-member ranges, so the actual
              % value is e raised to the value pulled from the distribution.             
              if isequal(distEMdist0{ii,jj},'LGU') | isequal(distEMdist0{ii,jj},'EPS-LGU')
                  EMIterInst(ii,jj) = exp(random(TrEMdist{ii,jj},1,1))*lgu_sign(ii,jj);
              end          
            
              % next, consider the case where the end-member chemistries are defined relative to the samples. 
              % this could include a normal distribution offset from the value of each sample or where a 
              % uniform or log-uniform distribution has one or both of its end-member values relative to the sample.
              if   isequal(distEMdist0{ii,jj},'SMP-NOR') | isequal(distEMdist0{ii,jj},'SMP-UNI') | isequal(distEMdist0{ii,jj},'SMP-LGU')
             
                   % because these distributions are defined relative to the sample, we have to isolate the relevant chemical parameter for 
                   % the river data. in the case of isotopes, convert the product of abundance and isotopic ratios to the isotopic data. 
                   if     ~ismember(ii,indxisopos);                                         riverentry = RiverColumn0(ii);
                   elseif isequal(ObsList{ii},'d7Li');                                      riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'Li'));   % d7Li
                   elseif isequal(ObsList{ii},'d18O');                                      riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'SO4'));  % d18O
                   elseif isequal(ObsList{ii},'d26Mg');                                     riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'Mg'));   % d26Mg
                   elseif isequal(ObsList{ii},'d30Si');                                     riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'Si'));   % d30Si
                   elseif isequal(ObsList{ii},'d34S');                                      riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'SO4'));  % d34S
                   elseif isequal(ObsList{ii},'d42Ca');                                     riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'Ca'));   % d42Ca
                   elseif isequal(ObsList{ii},'d44Ca');                                     riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'Ca'));   % d44Ca
                   elseif isequal(ObsList{ii},'d56Fe');                                     riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'Fe'));   % d56Fe
                   elseif isequal(ObsList{ii},'Sr8786');                                    riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'Sr'));   % Sr8786
                   elseif isequal(ObsList{ii},'Os8788');                                    riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'Os'));   % Os8788
                   elseif isequal(ObsList{ii},'d98Mo');                                     riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'Mo'));   % d98Mo
                   elseif isequal(ObsList{ii},'d13C') & isequal(carbonisotopematch,'HCO3'); riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'HCO3')); % d13C for HCO3
                   elseif isequal(ObsList{ii},'Fmod') & isequal(carbonisotopematch,'HCO3'); riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'HCO3')); % Fmod for HCO3
                   elseif isequal(ObsList{ii},'d13C') & isequal(carbonisotopematch,'DIC');  riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'DIC'));  % d13C for DIC
                   elseif isequal(ObsList{ii},'Fmod') & isequal(carbonisotopematch,'DIC');  riverentry = RiverColumn0(ii)./RiverColumn0(ismember(ObsList,'DIC'));  % Fmod for DIC
                   end                                        
              
                   % for a normal distribution with mean offset from the value of the sample
                   if     isequal(distEMdist0{ii,jj},'SMP-NOR')
                          EMIterInst(ii,jj) = riverentry + smpadddist_mean(ii,jj) + randn*smpadddist_std(ii,jj);
                               
                   % if the distribution is uniform, with one or both end-members defined relative to the samples
                   elseif isequal(distEMdist0{ii,jj},'SMP-UNI')
                          % find the new minimum and maximum of the ranges
                          if     isequal(distEMdist0_min{ii,jj},'SMP'); minval = riverentry + smpadddist_min(ii,jj);
                          elseif isequal(distEMdist0_min{ii,jj},'UNI'); minval = smpadddist_min(ii,jj);
                          end
                          if     isequal(distEMdist0_max{ii,jj},'SMP'); maxval = riverentry + smpadddist_max(ii,jj);
                          elseif isequal(distEMdist0_max{ii,jj},'UNI'); maxval = smpadddist_max(ii,jj);
                          end
                          % if the minimum is still less than the maximum, pull the value
                          if      minval <= maxval
                                  EMIterInst(ii,jj) = minval + rand*(maxval-minval);
                          else
                                  disp(sprintf('warning! error in MEANDIR_PullEndMemberRatios! the requested minimum value for variable %s in end-member %s is larger than the maximum value. MEANDIR will call this entry NaN.',ObsList{ii},EMList0{jj}));
                                  EMIterInst(ii,jj) = NaN;
                          end  
               
                   % if the distribution is log-uniform, with one or both end-members defined relative to the samples
                   elseif isequal(distEMdist0{ii,jj},'SMP-LGU')
                          % find the new minimum and maximum of the ranges
                          if     isequal(distEMdist0_min{ii,jj},'SMP'); minval = riverentry + smpadddist_min(ii,jj);
                          elseif isequal(distEMdist0_min{ii,jj},'LGU'); minval = smpadddist_min(ii,jj);
                          end
                          if     isequal(distEMdist0_max{ii,jj},'SMP'); maxval = riverentry + smpadddist_max(ii,jj);
                          elseif isequal(distEMdist0_max{ii,jj},'LGU'); maxval = smpadddist_max(ii,jj);
                          end
                          if     minval <= maxval
                                 if minval<=0; minval = 1e-10; disp(sprintf('warning! check %s end-member range for variable "%s". minval is negative and log would be imaginary, so setting value equal to 1e-10. ',EMList0{jj},ObsList{ii})); end
                                 if maxval<=0; maxval = 1e-10; disp(sprintf('warning! check %s end-member range for variable "%s". maxval is negative and log would be imaginary, so setting value equal to 1e-10. ',EMList0{jj},ObsList{ii})); end
                                 v1 = log(minval);
                                 v2 = log(maxval);
                                 if v1<v2; lowval = v1; highval = v2; end
                                 if v2<v1; lowval = v2; highval = v1; end
                                 EMIterInst(ii,jj) = exp((lowval + rand*(highval-lowval)));
                          else
                                 disp(sprintf('error in MEANDIR_PullEndMemberRatios! the requested minimum value for variable %s in end-member %s is larger than the maximum value. MEANDIR will call this entry NaN.',ObsList{ii},EMList0{jj}));
                                 EMIterInst(ii,jj) = NaN;
                          end                    
                   end % end of check which type distribution defined relative to each sample this is
              end      % end of checking if this is a distribution defined relative to each sample
          end          % end of iteration on end-members (jj)
          end          % end of iteration on observations (ii)
       
          % (2) if requested by the user (if they have set BalanceEvaporite==1) and when the inversion contains an end-member called "evap", the following
          % lines of code balance the evaporite end-member to be a stoichiometric NaCl-(Ca,Mg,Sr)SO4 mineral. the SO4 ratio is set equal
          % to the sum of Ca, Mg, and Sr, while the Cl ratio is set equal to the sum of the other cation ratios.
          if BalanceEvaporite==1 & sum(ismember(EMList0,'evap'))>0
             % first, match the SO4 ratio to the sum of Ca+Mg+Sr
             if sum(ismember(ObsList,'SO4'))>0 & (sum(ismember(ObsList,'Ca'))>0 | sum(ismember(ObsList,'Mg'))>0 | sum(ismember(ObsList,'Sr'))>0)
                newSO4values = EMIterInst(ismember(ObsList,'Ca') | ismember(ObsList,'Mg') | ismember(ObsList,'Sr'),ismember(EMList0,'evap'));
                EMIterInst(ismember(ObsList,'SO4'),ismember(EMList0,'evap')) = sum(newSO4values);
             end            
             % second, match the Cl ratio to the sum of every other cation. 
             if sum(ismember(ObsList,'Cl'))>0                   
                newClvalue = EMIterInst(ismember(IonCharges','+') & ~(ismember(ObsList,'Ca') | ismember(ObsList,'Mg') | ismember(ObsList,'Sr')),ismember(EMList0,'evap'));
                EMIterInst(ismember(ObsList,'Cl'),ismember(EMList0,'evap')) = sum(newClvalue);
             end
          end
          
          % (3) if AllIonsExplicitlyResolved = 1 and the normalization is to a single value, scale the values
          if AllIonsExplicitlyResolved==1
             for i=1:length(EMList0)
                   charbalanceratio = sum(EMIterInst(ismember(ObsList,ObsInNormalization),i)); % this should equal 1 if the end-member is consistent for its normalization. however, by setting the chargebalanceratio equal to the actual user-defined
                                                                                               % end-member values, we allow for the end-member to be charge-balanced even if it is not internally consistent for its normalization. 
                   offiopos = find(ismember(ObsList,ListChargeClosure{i}));                            
                   % test that the variable used to ensure internal consistency of charge is not included in the normalization.
                   if sum(ismember(ObsInNormalization,ListChargeClosure{i}))>0
                      disp(sprintf('warning! the value of "ListChargeClosure" for end-member "%s" is "%s", which is included in the normalization.',EMList0{i},ListChargeClosure{i})) 
                      disp('please select an observation not included in the normalization to ensure charge consistency of the end-members.')
                   end                     
                   for j=1:length(ObsList)
                       % the subsequent code calculates the normalized ratio for a variable, which is not included in the normalization, to ensure charge balance of the
                       % end-member. to understand the origin of the calculation, consider an example in which the inversion includes Si, DIC, Ca, Mg, Na, K, Cl, SO4, and HCO3,
                       % normalized to either Na or Si, and where Cl is used to ensure overall charge balance. first, write a statement of charge balance:
                       % Ca + Mg + Na + K = Cl + SO4 + HCO3 (note that the neutral species DIC and Si do not appear at this point)
                       % second, add the normalization for these scenarios (Na or Si)
                       % ex 1: Ca/Na + Mg/Na + Na/Na + K/Na = Cl/Na + SO4/Na + HCO3/Na
                       % ex 2: Ca/Si + Mg/Si + Na/Si + K/Si = Cl/Si + SO4/Si + HCO3/Si
                       % third, add species to the left-hand side and remove to the right-hand side in order to substitute '1' on the left 
                       % ex 1: Na/Na = 1 = Cl/Na + SO4/Na + HCO3/Na - Ca/Na -  Mg/Na - K/Na 
                       % ex 2: Si/Si = 1 = Cl/Si + SO4/Si + HCO3/Si - Ca/Si - Mg/Si - Na/Si - K/Si + Si/Si
                       % finally, rearrange the equation to solve for the term used to ensure charge consistency, which here is the normalized Cl ratio:
                       % Cl/Na = 1 - SO4/Na - HCO3/Na + Ca/Na + Mg/Na + K/Na
                       % Cl/Si = 1 - SO4/Si - HCO3/Si + Ca/Si + Mg/Si + + K/Si + Na/Si - Si/Si
                       if j~=offiopos                               
                          if     ismember(IonCharges{j},'+') & sum(ismember(CationsInNormalization,ObsList{j}))==1; charbalanceratio = charbalanceratio;                     % (1) for a cation in the normalization, do not change the value             (Na in the first example)
                          elseif ismember(IonCharges{j},'+') & sum(ismember(CationsInNormalization,ObsList{j}))==0; charbalanceratio = charbalanceratio + EMIterInst(j,i);   % (2) for a cation not in the normalization, add the normalized ratio        (Ca in the first example, Na in the second)
                          elseif ismember(IonCharges{j},'-') & sum(ismember(AnionsInNormalization, ObsList{j}))==1; charbalanceratio = charbalanceratio - 2*EMIterInst(j,i); % (3) for an anion in the normalization, subtract twice the normalized value (not in the examples above)
                          elseif ismember(IonCharges{j},'-') & sum(ismember(AnionsInNormalization, ObsList{j}))==0; charbalanceratio = charbalanceratio - EMIterInst(j,i);   % (4) for an anion not in the normalization, submit the normalized value     (SO4, HCO3 in both above example)
                          elseif ismember(IonCharges{j},'0') & sum(ismember(NeutralInNormalization,ObsList{j}))==1; charbalanceratio = charbalanceratio - EMIterInst(j,i);   % (5) for a neutral in the normalization, subtract the normalized value      (Si in the second example example)
                          elseif ismember(IonCharges{j},'0') & sum(ismember(NeutralInNormalization,ObsList{j}))==0; charbalanceratio = charbalanceratio;                     % (6) for a neutral not in the normalization, do not change the value        (Si in the first example example)
                          end
                       end
                   end
                   if isequal(IonCharges{offiopos},'-'); sign = 1;  end % if the ion is an ion, the calculation is as presented
                   if isequal(IonCharges{offiopos},'+'); sign = -1; end % if the ion is a cation, multiply by -1. In the example above, this would be the case of using Mg or K instead of Cl
                   EMIterInst(offiopos,i) = charbalanceratio*sign;      % save the off-charge ion entry
             end 
          end       
                  
          % (4) if the normalization is to a single observation, such as Na, and not to the sum of multiple observations,
          % the matrix at this point should be interally consistent for all end-members. in this case, set
          % 'EMsetisnotacceptable' to equal 0 to end the while loop and continue with the inversion. 
          if ~isequal(NormalizationType,'SumObs')
              EMsetisnotacceptable = 0;
          end   
        
          % (5) if the normalization is to the the sum of multiple observations, 
          % the internal consistency of the end-members must now be ensured.
          % there are multiple steps to this process:
          %   (5.1) calculating a single ratio by mass balance for each end-member
          %   (5.2) re-setting the evaporite end-member to the desired chemistry
          %   (5.3) setting a single ratio to ensure charge balance
          %   (5.4) set values arbitrarily close to zero to equal zero
          %   (5.5) confirm that the end-member matrix meets expectations
          %   (5.6) check if the end-member chemistry has unexpected negative values
          if isequal(NormalizationType,'SumObs')     
              
              % (5.1) % calculate one position by mass balance to ensure internal conistency of each end-member for the normalization 
             for i=1:length(EMList0)
                     offion = ismember(ObsList,ListNormClosure{i});           % find the value that should be re-defined
                     onions = ismember(ObsList,ObsInNormalization) & ~offion; % find the values that should be preserved
                     EMIterInst(offion,i) = 1 - sum(EMIterInst(onions,i),1);  % calculate the new value to ensure internal consistency
             end                        
                       
             % (5.2) now adjust the evaporite chemistry to account for closure. in the case that the anions SO4 and Cl are NOT included
             % in the normalization variable, you simply update their values. in the case that they are, you need to simultaneously
             % solve the problem of keeping the evaporite end-member interally consistent and account for the full mineral stoichiometry
             if  BalanceEvaporite==1 & sum(ismember(EMList0,'evap'))>0 & sum(ismember(ObsList,'SO4'))>0 & sum(ismember(ObsList,'Cl'))>0
                 closureion = ListNormClosure{ismember(EMList0,'evap')}; % find the variable used for internal consistent of the evaporite end-member
                                                                         % whether this is the Ca-Mg-Sr-SO4 or Na-Cl system will determine how to solve the problem
                 if  (isequal(closureion,'Ca') | isequal(closureion,'Mg') | isequal(closureion,'Sr'))==1
                     % access this code if the closure ion is in the Ca-Mg-Sr system. this means that, in the step above,
                     % one of these three systems are varied in order to ensure mass balance of the evaporite end-member.
                     
                     if     sum(ismember(ObsInNormalization,'SO4'))==0
                            % if SO4 is not in the normalization, then re-set the SO4 ratio to equal the sum
                            % of the Ca, Mg, and Sr ratios following the above correction for internal consistency
                            newSO4values = EMIterInst(ismember(ObsList,'Ca') | ismember(ObsList,'Mg') | ismember(ObsList,'Sr'),ismember(EMList0,'evap'));
                            EMIterInst(ismember(ObsList,'SO4'),ismember(EMList0,'evap')) = sum(newSO4values);
                  
                     elseif sum(ismember(ObsInNormalization,'SO4'))==1
                            % this is the harder case in which SO4 is included in the normalization
                            % first, calculate the fractional amount of the normalization variable that still has to be explained (massremaining)
                            % this is the sum of the closure variable and SO4, massremaining = closure + SO4
                            massremaining = 1-sum(EMIterInst(ismember(ObsList,ObsInNormalization(~(ismember(ObsInNormalization,'SO4') | ismember(ObsInNormalization,closureion)))),ismember(EMList0,'evap')));
                            
                            % second, calculate the amount of SO4-paired charge already accounted for (chargeaccounted). 
                            % this quantity plus the closure variable must equal the SO4 ratio: chargeaccounted + closure = SO4
                            chargeaccounted = 0; 
                            if ~isequal(closureion,'Ca') & sum(ismember(ObsList,'Ca'))>0; chargeaccounted = chargeaccounted+EMIterInst(ismember(ObsList,'Ca'),ismember(EMList0,'evap')); end
                            if ~isequal(closureion,'Mg') & sum(ismember(ObsList,'Mg'))>0; chargeaccounted = chargeaccounted+EMIterInst(ismember(ObsList,'Mg'),ismember(EMList0,'evap')); end
                            if ~isequal(closureion,'Sr') & sum(ismember(ObsList,'Sr'))>0; chargeaccounted = chargeaccounted+EMIterInst(ismember(ObsList,'Sr'),ismember(EMList0,'evap')); end
                            
                            % combining the equations, we find that closure = (massremaining - chargeaccounted)/2
                            closure = 1/2*(massremaining-chargeaccounted);                         
                            EMIterInst(ismember(ObsList,closureion),ismember(EMList0,'evap')) = closure;     
                            EMIterInst(ismember(ObsList,'SO4'),ismember(EMList0,'evap'))      = chargeaccounted+closure;
                     end                     
                 else 
                     % access this code if the closure ion is in the Na-Cl-other system. this means that, in the step above,
                     % one of these systems are varied in order to ensure mass balance of the evaporite end-member
                     
                     if     sum(ismember(ObsInNormalization,'Cl'))==0
                            % if Cl is not in the normalization, then re-set the Cl ratio to equal the sum
                            % of the Na and other (such as K) ratios following correction for internal consistency
                            newClvalue = EMIterInst(ismember(IonCharges','+') & ~(ismember(ObsList,'Ca') | ismember(ObsList,'Mg') | ismember(ObsList,'Sr')),ismember(EMList0,'evap'));
                            EMIterInst(ismember(ObsList,'Cl'),ismember(EMList0,'evap')) = sum(newClvalue);
                            
                     elseif sum(ismember(ObsInNormalization,'Cl'))==1
                            % this is the harder case in which Cl is included in the normalization.
                            % first, calculate the fractional amount of the normalization variable that still has to be
                            % explained (massremaining). This is the sum of the closure variable and Cl, massremaining = closure + SO4.
                            massremaining = 1-sum(EMIterInst(ismember(ObsList,ObsInNormalization(~(ismember(ObsInNormalization,'Cl') | ismember(ObsInNormalization,closureion)))),ismember(EMList0,'evap')));
                            
                            % second, calculate the amount of Cl-paired charge already accounted for (chargeaccounted). 
                            % this quantity plus the closure variable must equal the Cl ratio: chargeaccounted + closure = SO4
                            chargeaccounted = sum(EMIterInst(ismember(IonCharges','+') & ~(ismember(ObsList,'Ca') | ismember(ObsList,'Mg') | ismember(ObsList,'Sr') | ismember(ObsList,closureion)),ismember(EMList0,'evap')));
                            
                            % combining the equations, we find that closure = (massremaining - chargeaccounted)/2
                            closure = 1/2*(massremaining-chargeaccounted);
                            EMIterInst(ismember(ObsList,closureion),ismember(EMList0,'evap')) = closure;
                            EMIterInst(ismember(ObsList,'Cl'),ismember(EMList0,'evap'))       = chargeaccounted+closure;
                     end  % end loop on whether Cl is in the normalization
                 end      % end loop on the closure variable
             end          % end loop on balance evaporite
    
             
             % (5.3) if all ions are explicitly resolved, account for overall charge balance by solving for 1 ratio using the other
             % ratios. the variable solved by charge balance should not be included in the normalization variable. 
             if AllIonsExplicitlyResolved==1    
                for i=1:length(EMList0)
                     charbalanceratio = sum(EMIterInst(ismember(ObsList,ObsInNormalization),i)); % this should equal 1 if the end-member is consistent for its normalization. however, by setting the chargebalanceratio equal to the actual user-defined
                                                                                                 % end-member values, we allow for the end-member to be charge-balanced even if it is not internally consistent for its normalization.
                     offiopos = find(ismember(ObsList,ListChargeClosure{i}));                            
                     % test that the variable used to ensure internal consistency of charge is not included in the normalization.
                     if sum(ismember(ObsInNormalization,ListChargeClosure{i}))>0
                        disp(sprintf('warning! the value of "ListChargeClosure" for end-member "%s" is "%s", which is included in the normalization.',EMList0{i},ListChargeClosure{i})) 
                        disp('please select an observation not included in the normalization to ensure charge consistency of the end-members.')
                     end                     
                     for j=1:length(ObsList)
                         % the subsequent code calculates the normalized ratio for a variable, which is not included in the normalization, to ensure charge balance of the
                         % end-member. to understand the origin of the calculation, consider in which the inversion includes Si, DIC, Ca, Mg, Na, K, Cl, SO4, and HCO3,
                         % normalized to the sum of Ca, Na, DIC, and SO4, and where Cl is used to ensure overall charge balance.  in this example, note that Cl is not 
                         % included in the normalization but that the normalization does contain cations, anions, and neutral species. first, write a statement of charge balance:
                         % Ca + Mg + Na + K = Cl + SO4 + HCO3 (note that the neutral species DIC and Si do not appear at this point)
                         % second, add the normalization for this scenario (Ca+Na+DIC+SO4)
                         % Ca/(Ca+Na+DIC+SO4) + Mg/(Ca+Na+DIC+SO4) + Na/(Ca+Na+DIC+SO4) + K/(Ca+Na+DIC+SO4) = Cl/(Ca+Na+DIC+SO4) + SO4/(Ca+Na+DIC+SO4) + HCO3/(Ca+Na+DIC+SO4)
                         % third, add species to the left-hand side and remove to the right-hand side in order to substitute '1' on the left                     % 
                         % Ca/(Ca+Na+DIC+SO4) + Na/(Ca+Na+DIC+SO4) +  SO4/(Ca+Na+DIC+SO4) + DIC/(Ca+Na+DIC+SO4) = Cl/(Ca+Na+DIC+SO4) + 2*SO4/(Ca+Na+DIC+SO4) + HCO3/(Ca+Na+DIC+SO4) - Mg/(Ca+Na+DIC+SO4) -  K/(Ca+Na+DIC+SO4) + DIC/(Ca+Na+DIC+SO4)
                         % 1 = Cl/(Ca+Na+DIC+SO4) + 2*SO4/(Ca+Na+DIC+SO4) + HCO3/(Ca+Na+DIC+SO4) - Mg/(Ca+Na+DIC+SO4) -  K/(Ca+Na+DIC+SO4) + DIC/(Ca+Na+DIC+SO4)
                         % finally, rearrange the equation to solve for the term used to ensure charge consistency, which here is the normalized Cl ratio:
                         % Cl/(Ca+Na+DIC+SO4) = 1 - 2*SO4/(Ca+Na+DIC+SO4) - HCO3/(Ca+Na+DIC+SO4) + Mg/(Ca+Na+DIC+SO4) + K/(Ca+Na+DIC+SO4) - DIC/(Ca+Na+DIC+SO4)
                         % to generalize the above example, start with 1 and then add or subtract values depending on whether the ion is a cation, anion, or neutral, and
                         % depending on whether it is included or excluded from the normalization.                      
                         if j~=offiopos                                                                                        
                            if     ismember(IonCharges{j},'+') & sum(ismember(CationsInNormalization,ObsList{j}))==1; charbalanceratio = charbalanceratio;                     % (1) for a cation in the normalization, do not change the value             (Ca in the above example)
                            elseif ismember(IonCharges{j},'+') & sum(ismember(CationsInNormalization,ObsList{j}))==0; charbalanceratio = charbalanceratio + EMIterInst(j,i);   % (2) for a cation not in the normalization, add the normalized ratio        (Mg in the above example)
                            elseif ismember(IonCharges{j},'-') & sum(ismember(AnionsInNormalization, ObsList{j}))==1; charbalanceratio = charbalanceratio - 2*EMIterInst(j,i); % (3) for an anion in the normalization, subtract twice the normalized value (SO4 in the above example)
                            elseif ismember(IonCharges{j},'-') & sum(ismember(AnionsInNormalization, ObsList{j}))==0; charbalanceratio = charbalanceratio - EMIterInst(j,i);   % (4) for an anion not in the normalization, submit the normalized value     (HCO3 in the above example)
                            elseif ismember(IonCharges{j},'0') & sum(ismember(NeutralInNormalization,ObsList{j}))==1; charbalanceratio = charbalanceratio - EMIterInst(j,i);   % (5) for a neutral in the normalization, subtract the normalized value      (DIC in the above example)
                            elseif ismember(IonCharges{j},'0') & sum(ismember(NeutralInNormalization,ObsList{j}))==0; charbalanceratio = charbalanceratio;                     % (6) for a neutral not in the normalization, do not change the value        (Si in the above example)
                            end
                         end
                      end
                      if isequal(IonCharges{offiopos},'-'); sign = 1;  end % if the ion is an ion, the calculation is as presented
                      if isequal(IonCharges{offiopos},'+'); sign = -1; end % if the ion is a cation, multiply by -1. In the example above, this would be the case of using Mg or K instead of Cl
                      EMIterInst(offiopos,i) = charbalanceratio*sign;      % save the off-charge ion entry
                end 
             end                                      
                                          
             % (5.4) set negative things arbitrarily close to zero to equal zero
             % this is required because sometimes values like -1e-10, which is effectively zero, will be read as negative.
             EMIterInst(EMIterInst>-1e-12 & EMIterInst<0)=0; 
                        
             % (5.5) check that the end-member matrix fits expectations:     
             % (5.5.1) check for mass balance
             if ~sum(sum(EMIterInst(ismember(ObsList,ObsInNormalization),:))==1)==nEM
                disp('warning! end-members are internally inconsistent for mass balance!'); 
                FailCase = 1;
             end                 
              % (5.5.2) check for evaporite balance
             if BalanceEvaporite==1 & sum(ismember(EMList0,'evap'))>0
                if sum(ismember(ObsList,'SO4'))>0 & (sum(ismember(ObsList,'Ca'))>0 | sum(ismember(ObsList,'Mg'))>0 | sum(ismember(ObsList,'Sr'))>0)
                   CatSO4values = EMIterInst(ismember(ObsList,'Ca') | ismember(ObsList,'Mg') | ismember(ObsList,'Sr'),ismember(EMList0,'evap'));
                   SO4values    = EMIterInst(ismember(ObsList,'SO4'),ismember(EMList0,'evap'));
                   if sum(CatSO4values)~=SO4values
                       FailCase = 1;
                       disp('warning! normalized evaporite Ca+Mg+Sr is not equal to normalized SO4! this error is identified in "MEANDIR_PullEndMemberRatios".'); FailCase = 1; 
                   end
                end
                if sum(ismember(ObsList,'Cl'))>0                   
                   CatClvalues = EMIterInst(ismember(IonCharges','+') & ~(ismember(ObsList,'Ca') | ismember(ObsList,'Mg') | ismember(ObsList,'Sr')),ismember(EMList0,'evap'));
                   Clvalues    = EMIterInst(ismember(ObsList,'Cl'),ismember(EMList0,'evap'));
                   if  round(sum(CatClvalues),10)~=round(Clvalues,10)
                        FailCase = 1;
                       disp('warning! normalized evaporite Na is not equal to normalized non-Ca+Mg+Sr! this error is identified in "MEANDIR_PullEndMemberRatios".'); FailCase = 1;
                   end
                end
             end    
             % (5.5.3) check for charge balance
             if AllIonsExplicitlyResolved==1 
                 catpos = ismember(IonCharges','+');
                 anipos = ismember(IonCharges','-');
                 cb = sum(EMIterInst(catpos,:)) - sum(EMIterInst(anipos,:)); 
                 if ~sum(cb==0)==nEM
                     disp('warning! end-members are internally inconsistent for charge balance! this error is identified in "MEANDIR_PullEndMemberRatios".');  
                     FailCase = 1;
                 end
             end            
             
             % (5.6) check that nothing in EMIterInst that should be positive is negative.
              % if the matrix is set for inversion, set 'stillbad' to 0 to continue with the inversion 
             if  sum(sum(EMIterInst(~isotopeposition,~ismember(EMList0,EndMembersWithNegativeRatios))<0))==0
                 EMsetisnotacceptable = 0; 
             end
         end % end corrections for when normalizing by SumCations
               
         % (6) if MEANDIR is struggling to generate reasonable end-members given the inputs, diplay some information to the user
         if  trycount==maxtrycount
             EMsetisnotacceptable = 0;
             disp(sprintf('\non attempt %i/%i, sample %i, MEANDIR was not able to generate a viable set of end-members given the number of tries specified by "maxtrycount" (currently set to: %i)',iterationcounter,maxiterations,i,maxtrycount))
             disp('if internally consistent end-members are expected to be difficult to find, such as when many elements are included in the normalization, raise the value of "maxtrycount" within the program "MEANDIR_PullEndMemberRatios"')
             disp('the last attempted matrix of end-member ratios is given below. try to identify any errors:')
             EMIterInst
                                       
             negEMlist = EMList0(sum(EMIterInst<0,1)>0);
             if ~isempty(negEMlist)
                endmemberOneglist = []; for lll=1:length(negEMlist); endmemberOneglist = [endmemberOneglist negEMlist{lll} ', ']; end; endmemberOneglist = endmemberOneglist(1:end-2);
                endmemberKneglist = []; for lll=1:length(EndMembersWithNegativeRatios); endmemberKneglist = [endmemberKneglist EndMembersWithNegativeRatios{lll} ', ']; end; endmemberKneglist = endmemberKneglist(1:end-2);
                if sum(ismember(negEMlist,EndMembersWithNegativeRatios))~=length(negEMlist)
                    disp(sprintf('problem! MEANDIR has found that the following end-members have negative ratios: %s',endmemberOneglist))
                    disp(sprintf('currently only the following end-members are allowed to have negative ratios: %s',endmemberKneglist))
                    disp('note that MEANDIR does not accept negative end-member ratios except where defined in "EndMembersWithNegativeRatios".')
                end
             else
                disp('if this error occurs for all inversion attempts, there may be a problem with the selection of end-members that is preventing physically-realizable compositions')
                disp('if normalizing to a sum, make sure the ratio solved by mass balance can be positive')
                disp('if all cations and anions are included in the inversion, ensure that the ratio solved by charge balance can be positive')
             end                         
             
             EMIterInst = nan(size(EMIterInst));          
             listallNaN = 0;
             
         end  % end of displaying information when trycount>maxtrycount
     end      % end of while loop  
           
     % (7) if H2SO4 should be combined with rock weathering and those end-members given the same d34S, 
     % define the d34S of the FeS2 using the "pyri" end-member of the "MEANDIR_Endmembers" spreadsheet
     % the reason this calculation is done here, as opposed to when initially defining the end-members, is
     % that it is desirable for all end-member where SO4 represents FeS2 oxidation to have the same d34S.     
      if sum(ismember(ObsList,'SO4'))== 1 & sum(ismember(ObsList,'d34S'))== 1 & ~isempty(CoupleFeS2d34SintoEM)
        % match pyrite d34S values across end-members
        pos_obs_d34S = ismember(ObsList,'d34S'); 
        pos_ems_d34S = ismember(EMList0,CoupleFeS2d34SintoEM);
        if  sum(pos_obs_d34S) > 0
            disttype = eval(sprintf('EM.%s.disttype.FeS2d34S',EMdatasource));
            % (7.1) if FeS2 d34S has a normal distribution
            if     isequal(disttype,'NOR')
                   pyr_d34S_mean = eval(sprintf('EM.%s.pyri.%s_MenFeS2d34S%s',EMdatasource,disttype,NormalizationType));
                   pyr_d34S_std  = eval(sprintf('EM.%s.pyri.%s_SigFeS2d34S%s',EMdatasource,disttype,NormalizationType));
                   if   sum(ismember(ConvertDelta2RList,'d34S'))>0
                        conv        = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d34S'));
                        record_mean = (pyr_d34S_mean/1000+1)*conv;
                        record_std  = (pyr_d34S_std/1000)*conv;
                   else
                        record_mean  = pyr_d34S_mean;
                        record_std   = pyr_d34S_std;
                   end
                   EMIterInst(pos_obs_d34S,pos_ems_d34S) = record_mean + record_std*randn;
            % (7.2) if FeS2 d34S has a uniform distribution
            elseif isequal(disttype,'UNI')
                   pyr_d34S_min = eval(sprintf('EM.%s.pyri.%s_MinFeS2d34S%s',EMdatasource,disttype,NormalizationType));
                   pyr_d34S_max = eval(sprintf('EM.%s.pyri.%s_MaxFeS2d34S%s',EMdatasource,disttype,NormalizationType));
                   if   sum(ismember(ConvertDelta2RList,'d34S'))>0
                        conv        = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d34S'));
                        record_min  = (pyr_d34S_min/1000+1)*conv;
                        record_max  = (pyr_d34S_max/1000+1)*conv;
                   else
                        record_min  = pyr_d34S_min;
                        record_max  = pyr_d34S_max;
                   end
                   EMIterInst(pos_obs_d34S,pos_ems_d34S) = record_min + rand*(record_max-record_min);
            end
        end
      end

     % (8) adjust the EM matrix for any isotopes in the system, so that its entries are now the product of
     % the isotope ratio and the concentration information (this product is what mixes approximately linearly)
     for j=1:sum(isotopeposition) % iterate over the number of isotope ratio in the system   
          activeisopos = indxisopos(j);
         activeisoname = ObsList{activeisopos};
         if isequal(activeisoname,'d7Li');                                       activeionpos = find(ismember(ObsList,'Li'));   end % if d7Li,   find Li
         if isequal(activeisoname,'d18O');                                       activeionpos = find(ismember(ObsList,'SO4'));  end % if d18O,   find SO4
         if isequal(activeisoname,'d26Mg');                                      activeionpos = find(ismember(ObsList,'Mg'));   end % if d26Mg,  find Mg
         if isequal(activeisoname,'d30Si');                                      activeionpos = find(ismember(ObsList,'Si'));   end % if d30Si,  find Si
         if isequal(activeisoname,'d34S');                                       activeionpos = find(ismember(ObsList,'SO4'));  end % if d34S,   find SO4
         if isequal(activeisoname,'d42Ca');                                      activeionpos = find(ismember(ObsList,'Ca'));   end % if d42Ca,  find Ca
         if isequal(activeisoname,'d44Ca');                                      activeionpos = find(ismember(ObsList,'Ca'));   end % if d44Ca,  find Ca
         if isequal(activeisoname,'d56Fe');                                      activeionpos = find(ismember(ObsList,'Fe'));   end % if d56Fe,  find Fe
         if isequal(activeisoname,'Sr8786');                                     activeionpos = find(ismember(ObsList,'Sr'));   end % if Sr8786, find Sr
         if isequal(activeisoname,'Os8788');                                     activeionpos = find(ismember(ObsList,'Os'));   end % if Os8788, find Os
         if isequal(activeisoname,'d98Mo');                                      activeionpos = find(ismember(ObsList,'Mo'));   end % if d98Mo,  find Mo
         if isequal(activeisoname,'d13C') & isequal(carbonisotopematch,'HCO3');  activeionpos = find(ismember(ObsList,'HCO3')); end % if d13C and carbon type is HCO3
         if isequal(activeisoname,'Fmod') & isequal(carbonisotopematch,'HCO3');  activeionpos = find(ismember(ObsList,'HCO3')); end % if Fmod and carbon type is HCO3
         if isequal(activeisoname,'d13C') & isequal(carbonisotopematch,'DIC');   activeionpos = find(ismember(ObsList,'DIC'));  end % if d13C and carbon type is DIC
         if isequal(activeisoname,'Fmod') & isequal(carbonisotopematch,'DIC');   activeionpos = find(ismember(ObsList,'DIC'));  end % if Fmod and carbon type is DIC
         EMIterInst(activeisopos,:) = EMIterInst(activeisopos,:).*EMIterInst(activeionpos,:);
     end  
    
     % (9) in case NaN values made it this far, tell the user
     if sum(isnan(EMIterInst(:)))>0 & listallNaN==1
        [a, b] = find(isnan(EMIterInst));
        if     length(a)==1
               disp(sprintf('warning! %i entry of the end-member matrix (EMIterInst) is NaN.',length(a)));   
        elseif length(a)>1  
               disp(sprintf('warning! %i entries of the end-member matrix (EMIterInst) are NaN.',length(a))); 
        end
        for aa=1:length(a)
            disp(sprintf('         observation "%s" in end-member "%s" is NaN',ObsList{a(aa)},EMList0{b(aa)}));
        end   
     end    
end