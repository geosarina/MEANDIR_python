function [RiverColumn_r, EMIterInst_r, EMList_r, ObsList_r, Xdirect, distEMdist_r] = MEANDIR_PrepareUpdatedDataForInversion(Xdirect,RiverColumn,ObsList,EMList,EMList0,relpos,EMIterInst,Solver, distEMdist)

         % this function prepares data post-ClCritical correction for inversion. in the case of absolute error functions, 
         % nothing must be removed from the inversion and you can simply reset the Cl-corrected data to have "_r" in the variable names.
         % in the case of optimization, where there are zero values evaluated with relative error, remove those entries of the river
         % likewise, remove from the matrix of end-members the row corresponding to this value, and set any end-members with
         % positive contributions to that parameter to have a fractional contribution of zero to the normalization ion. 

         if     isequal(Solver,'mldivide') | isequal(Solver,'lsqnonneg')

                RiverColumn_r = RiverColumn;
                EMIterInst_r  = EMIterInst;
                EMList_r      = EMList;
                ObsList_r     = ObsList;
                distEMdist_r  = distEMdist;
             
         elseif isequal(Solver,'mldivide_optimize') | isequal(Solver,'lsqnonneg_optimize') | isequal(Solver,'optimize')

                % (1) if the river column has a zero entry and should be evaluated with relative
                % misfit, remove that entry from the inversion and record which positions were zero
                RiverColumnZeroRemove_index = (RiverColumn==0 & relpos);                                          % get the index of zero values for the river
                RiverColumnZeroRemove_pos   = find(RiverColumnZeroRemove_index);                                  % get the position of those zero values
                if sum(RiverColumnZeroRemove_index)>0; RiverColumn_r = RiverColumn(~RiverColumnZeroRemove_index); % if there are zero values, remove them from the RiverColumn 
                else;                                  RiverColumn_r = RiverColumn;                               % if there are no zero values, re-define the RiverColumn
                end 

                % (2) address the case where the river has a zero value evaluated with 
                % relative error, but there are end-members attempting to source that ratio. 
                EMcol_removelist = [];        
                for jj=1:length(RiverColumnZeroRemove_pos)            % iterate over index where river data is zero
                     activeposition  = RiverColumnZeroRemove_pos(jj); % get the active row
                     EMnonzerocolumn = find(EMIterInst(activeposition,:)~=0);
                     if     isempty(EMnonzerocolumn)
                            % if EMcol is empty, the river and EM matrix are all zero 
                            % in this row, and no fractional contributions can be determined
                     elseif ~isempty(EMnonzerocolumn) 
                            % if EMcol is not empty, end-members (columns) with non-zero values in the
                            % active position where the river is zero cannot be contributing to the river
                            for kk=1:length(EMnonzerocolumn)
                                Xdirect(ismember(EMList0,EMList{EMnonzerocolumn(kk)})) = 0;
                                EMcol_removelist = [EMcol_removelist EMnonzerocolumn(kk)]; % this allows for multiple end-members to be removed
                            end
                     end
                end

                % (3) make "_r" version of the variables
                EMrmve       = ismember(1:length(EMList),EMcol_removelist);
                EMIterInst_r = EMIterInst(~RiverColumnZeroRemove_index,~EMrmve);
                EMList_r     = EMList(~EMrmve);
                ObsList_r    = ObsList(~RiverColumnZeroRemove_index);
                distEMdist_r = distEMdist(~RiverColumnZeroRemove_index,~EMrmve);

         end  % checking the value of the solver

         % quality contol check asking if the reduced end-member matrix has any zero values
         if  sum(sum(EMIterInst_r,2)==0) > 0
             findoff = find(sum(EMIterInst_r,2)==0);
             if  isequal(Solver,'mldivide_optimize')
                 disp(' ')
                 for jj=1:length(findoff)
                     disp(sprintf('warning! the reduced matrix of end-members has only zero values for %s',ObsList_r{findoff(jj)}))
                 end
                 disp('this situation may result from asking the inversion to invert a sample with zero concentration using a relative cost function') 
                 disp('as a result, the matrix for inversion may be singular') 
                 disp('this error arises in "MEANDIR_PrepareUpdatedDataForInversion"') 
                 disp(' ')
             end
         end

end           % end of the function