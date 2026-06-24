function cost = MEANDIR_CostFunction(X0,EMIterInst0, EMIterInst_r, RiverColumn0, RiverColumn_r, SolveCFList_r, abspos_r, relpos_r, WeightingList_r, fractionation, Xdirect, sources)

         % the cost function solved during optimization. this function defines the actual inversion
         % step, in which matlab optimizes a set of X values to minimize the function evaluation.
         % without fractionation, this script is just a cost function
         % with fractionation, the isotopic composition of certain end-members
         % has to be updated to reflect the inputs of all other end-members

         for i = 1:fractionation.red.n_r
             % access this code in the case that there is fractionation in the model. iterate over each pair of isotope and end-member
             % start with the following equations, where Fclay is negative and a negative Delta indicates a normal isotope effect:
             %   deta_obs   = delta_source + Fclay/Sum(Fsources)*(Delta) 
             %   delta_clay = (delta_source + Delta) + Fclay/Sum(Fsources)*(Delta)

             Xfull                  = Xdirect;                                                                                                                                       % set X as the known solution vector for the normalization ion
             Xfull(isnan(Xdirect))  = X0;                                                                                                                                            % fill in gaps in X with the current set of fractional contributions
             Xactive                = EMIterInst0(fractionation.red.ionpos0(i),:).*Xfull'./RiverColumn0(fractionation.red.ionpos0(i));                                               % convert the X values to fraction of the active ion
             endmemberisotopevalues = EMIterInst0(fractionation.red.isopos0(i),:)./EMIterInst0(fractionation.red.ionpos0(i),:);                                                      % calculate the isotopic information of the end-members
             sourceindx             = sources & Xactive~=0;                                                                                                                          % find the position of sources with non-zero contributions
             sourceisotopevalue     = sum(Xactive(sourceindx).*endmemberisotopevalues(sourceindx))./sum(Xactive(sourceindx));                                                        % calculate the isotopic information value of the inputs
             activefractionation    = EMIterInst0(fractionation.red.isopos0(i),fractionation.red.empos0(i))./EMIterInst0(fractionation.red.ionpos0(i),fractionation.red.empos0(i));  % calculate the current fractionation
             newisotopevalue        = (sourceisotopevalue + activefractionation) + activefractionation.*Xactive(fractionation.red.empos0(i))/sum(Xactive(sourceindx));               % calculate the isotope information of the secondary phase
             EMIterInst_r(fractionation.red.isopos_r(i),fractionation.red.empos_r(i)) = EMIterInst_r(fractionation.red.ionpos_r(i),fractionation.red.empos_r(i)).*newisotopevalue;   % replace the value of the secondary phase in the inversion matrix
         end

         % now calculate the cost function
         cost_relative = ((RiverColumn_r - EMIterInst_r*X0)./RiverColumn_r).^2; % calculate the cost function for relative change
         cost_absolute = (RiverColumn_r  - EMIterInst_r*X0).^2;                 % calculate the cost function for absolute change

         % sum the cost function values at the user-requested positions
         cost          = sqrt(sum([WeightingList_r(SolveCFList_r & relpos_r).*cost_relative(SolveCFList_r & relpos_r); WeightingList_r(SolveCFList_r & abspos_r).*cost_absolute(SolveCFList_r & abspos_r)])); 

end      % end of function