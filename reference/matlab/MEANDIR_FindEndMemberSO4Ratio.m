function EMvalue = MEANDIR_FindEndMemberSO4Ratio(EM,EMList0,EMdatasource,NormalizationType,j, EndMembersWithNegativeRatios)

        % this script finds SO4 end-member ratios when SO4 is not included
        % in the inversion (such as when calculating SO4 excess)
        % this function is called from MEANDIR_EvaluateInversionInstance

        % first, identify the type of SO4 distribution associated with this end-member
        disttype = eval(sprintf('EM.%s.disttype.SO4',EMdatasource));

        % if the end-member SO4 distribution is normal
        if     isequal(disttype,'NOR')
               EMvaluemean   = eval(sprintf('EM.%s.%s.NOR_MenSO4%s',EMdatasource,EMList0{j},NormalizationType)); % mean of normal distribution
               EMvaluesig    = eval(sprintf('EM.%s.%s.NOR_SigSO4%s',EMdatasource,EMList0{j},NormalizationType)); % standard deviation of normal distribution
               if    EMvaluesig == 0
                     EMvalue    = EMvaluemean; 
               else
                     pd = makedist('Normal','mu',EMvaluemean,'sigma',EMvaluesig);
                     if  ismember(EMList0{j},EndMembersWithNegativeRatios)==0
                         % if this end-member cannot have negative chemical
                         % ratios, truncate the distribution
                         tr = truncate(pd,0,inf);
                         EMvalue = random(tr,1);
                     else
                         EMvalue = random(pd,1);
                     end
               end

        % if the end-member SO4 distribution is uniform
        elseif isequal(disttype,'UNI')
               EMvaluemin = eval(sprintf('EM.%s.%s.UNI_MinSO4%s',EMdatasource,EMList0{j},NormalizationType)); % pull the minimum value
               EMvaluemax = eval(sprintf('EM.%s.%s.UNI_MaxSO4%s',EMdatasource,EMList0{j},NormalizationType)); % pull the maximum value
               EMvalue    = EMvaluemin + (EMvaluemax-EMvaluemin)*rand(1);                                     % get a value pulled from that uniform distribution

        % if the end-member SO4 distribution is log-uniform
        elseif isequal(disttype,'LGU')        
               EMvaluemin = eval(sprintf('EM.%s.%s.LGU_MinSO4%s',EMdatasource,EMList0{j},NormalizationType)); % pull the minimum value
               EMvaluemax = eval(sprintf('EM.%s.%s.LGU_MaxSO4%s',EMdatasource,EMList0{j},NormalizationType)); % pull the maximum value

               if   EMvaluemin==EMvaluemax
                    EMvalue = (EMvaluemin+EMvaluemax)/2; % if the values are equal, there is only a single value.
               else                                                                            
                    if   EMvaluemin==0 | EMvaluemax==0   % if only one value is zero, a log-uniform distribution is not well-defined
                         disp('  '); 
                         disp('error! cannot define a log-normal distribution when only one of the values is zero.'); 
                         disp(sprintf('check variable "SO4" in end-member "%s"',EMList0{j}));
                         disp('  ');         
                         EMvalue = NaN; 
                    else    
                         % proceed in the case that the minimum is less than the maximum
                         if   EMvaluemin<EMvaluemax
                              if   EMvaluemin<0 | EMvaluemax < 0; isneg = 1; mult_ifneg = -1; % note if this chemical range is negative
                              else;                               isneg = 0; mult_ifneg = +1; % note if this chemical range is positive
                              end
                              if   isneg == 1 
                                   v1 = log(-EMvaluemin); 
                                   v2 = log(-EMvaluemax);
                                   if v1<v2; lowval = v1; highval = v2; end
                                   if v2<v1; lowval = v2; highval = v1; end
                              else
                                   v1 = log(EMvaluemin); 
                                   v2 = log(EMvaluemax);
                                   if v1<v2; lowval = v1; highval = v2; end
                                   if v2<v1; lowval = v2; highval = v1; end
                              end
                              EMvalue = exp(lowval+(highval-lowval)*rand(1))*mult_ifneg;
                         else                    
                              disp(sprintf('\nproblem with end-members! the lower normalized ratio for SO4 in %s is higher than its higher normalized ratio.',EMList0{j}));
                         end
                    end
               end  % end of "if EMvaluemin==EMvaluemax"

        end % end if statement on the type of distribution
end         % end function