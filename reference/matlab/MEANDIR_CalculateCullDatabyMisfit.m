% this script culls simulations to isolate the fraction of results with lowest misfit between model and observations
% when IterateOver is "End-members", the script saves the fraction of results specified with the user-defined variable 
% "MisfitCuts", and evaluates misfit at the scale of either all samples or each sample individually as specified by 
% the user-defined variable "CullOn". when IterateOver is "Samples", the script sets remaininginstances = maxsuccess

if isequal(IterateOver,'End-members') % only cut the data by MisfitCuts when IterateOver is "End-members"
   
    if     isequal(CullOn,'AllSample')          
           % calculate the error across all of these samples
           erroracrosssamples         = sqrt(nansum(results.misfit.^2,1)); 
           [sortedlist, sortindx]     = sort(erroracrosssamples);    
           remaininginstances         = floor(numberiterations*MisfitCuts/100);
           kpindx                     = ismember(1:numberiterations,sortindx(1:remaininginstances));

           % update the matrix of end member chemistry
           EMSucessIterMatrixAll      = EMSucessIterMatrixAll(:,:,kpindx,:);
           EMSucessIterMatrixAll_frac = EMSucessIterMatrixAll_frac(:,:,kpindx,:);
  
           % update the river data
           fn = fieldnames(InvRivDataAgg);
           for i=1:length(fn)  
               currentdata      = eval(sprintf('InvRivDataAgg.%s',fn{i}));
               InvRivDataAgg    = setfield(InvRivDataAgg,fn{i},currentdata(:,kpindx));
               if saveuncutdata == 1
                  InvRivDataAgg = setfield(InvRivDataAgg,'all_uncut',fn{i},currentdata);
               end
           end

           % update the results structure
           fn = fieldnames(results);
           for i=1:length(fn)
               currentdata      = eval(sprintf('results.%s',fn{i}));
               results          = setfield(results,fn{i},currentdata(:,kpindx,:));
               if saveuncutdata == 1
                  results       = setfield(results,'all_uncut',fn{i},currentdata);
               end
           end

           successfullsimulations  = sum(successfullsimulations(:,kpindx),2);
           badmbcount              = sum(badmbcount(:,kpindx),2);
           iterationcounts         = sum(iterationcounts(:,kpindx),2);

   elseif isequal(CullOn,'EachSample')

           % start by defining some useful working vectors and matricies
           EMSucessIterMatrixAll_temp      = [];
           EMSucessIterMatrixAll_frac_temp = [];
           successfullsimulations_temp     = [];
           badmbcount_temp                 = [];
           iterationcounts_temp            = [];

           % make a temporary structure to hold InvRivDataAgg
           InvRivDataAgg_temp = struct;
           fn = fieldnames(InvRivDataAgg);
           for i=1:length(fn)
               InvRivDataAgg_temp = setfield(InvRivDataAgg_temp,fn{i},NaN(s,floor(numberiterations*MisfitCuts/100)));
           end
           if saveuncutdata==1
              for i=1:length(fn)
                  InvRivDataAgg_temp = setfield(InvRivDataAgg_temp,'all_uncut',fn{i},NaN(s,numberiterations));
              end
           end

           % make a temporary structure to hold results
           results_temp = struct;
           fn = fieldnames(results);
           fn = fn(~(ismember(fn,'misfit')));
           fn = fn(~(ismember(fn,'misfit_costfunction')));
           for i=1:length(fn)
               results_temp = setfield(results_temp,fn{i},NaN(nEM,floor(numberiterations*MisfitCuts/100),s));
           end
           results_temp.misfit              = NaN(s,floor(numberiterations*MisfitCuts/100));
           results_temp.misfit_costfunction = NaN(s,floor(numberiterations*MisfitCuts/100));
           if saveuncutdata==1
              for i=1:length(fn)
                  results_temp                               = setfield(results_temp,'all_uncut',fn{i},NaN(nEM,numberiterations,s));
                  results_temp.all_uncut.misfit              = NaN(s,numberiterations);
                  results_temp.all_uncut.misfit_costfunction = NaN(s,numberiterations);
              end
           end

           % iterate over each sample
           for k=1:s
                % calculate the error across all of these samples
                erroracrosssample      = results.misfit(k,:);
                [sortedlist, sortindx] = sort(erroracrosssample);
                remaininginstances     = floor(numberiterations*MisfitCuts/100);
                kpindx                 = ismember(1:numberiterations,sortindx(1:remaininginstances));

                % update the end-members
                EMSucessIterMatrixAll_temp(:,:,1:remaininginstances,k)      = EMSucessIterMatrixAll(:,:,kpindx,k);
                EMSucessIterMatrixAll_frac_temp(:,:,1:remaininginstances,k) = EMSucessIterMatrixAll_frac(:,:,kpindx,k);

                % update the river data
                fn = fieldnames(InvRivDataAgg);
                for i=1:length(fn)
                    currentdata           = eval(sprintf('InvRivDataAgg.%s',fn{i}));
                    currentdata           = currentdata(k,:);
                    currentsaveddata      = eval(sprintf('InvRivDataAgg_temp.%s',fn{i}));
                    currentsaveddata(k,:) = currentdata(kpindx);
                    InvRivDataAgg_temp    = setfield(InvRivDataAgg_temp,fn{i},currentsaveddata);
                    if saveuncutdata == 1
                       currentallsaveddata      = eval(sprintf('InvRivDataAgg_temp.all_uncut.%s',fn{i}));
                       currentallsaveddata(k,:) = currentdata;
                       InvRivDataAgg_temp       = setfield(InvRivDataAgg_temp,'all_uncut',fn{i},currentallsaveddata);
                    end
                end

                % update the results structure
                fn = fieldnames(results);
                fn = fn(~(ismember(fn,'misfit')));
                fn = fn(~(ismember(fn,'misfit_costfunction')));
                for i=1:length(fn)
                    currentresults             = eval(sprintf('results.%s',fn{i}));
                    currentresults             = currentresults(:,:,k);
                    currentresultsddata        = eval(sprintf('results_temp.%s',fn{i}));
                    currentresultsddata(:,:,k) = currentresults(:,kpindx);
                    results_temp = setfield(results_temp,fn{i},currentresultsddata);
                    if saveuncutdata==1
                        currentallresddata        = eval(sprintf('results_temp.all_uncut.%s',fn{i}));
                        currentallresddata(:,:,k) = currentresults;
                        results_temp              = setfield(results_temp,'all_uncut',fn{i},currentallresddata);
                    end
                end

                % for results, now do the misfit calculation
                currentresults           = results.misfit;
                currentresults           = currentresults(k,:);
                currentresultsddata      = results_temp.misfit;
                currentresultsddata(k,:) = currentresults(:,kpindx);
                results_temp.misfit      = currentresultsddata;
                if saveuncutdata == 1
                   currentresultsddata           = results_temp.all_uncut.misfit;
                   currentresultsddata(k,:)      = currentresults;
                   results_temp.all_uncut.misfit = currentresultsddata;
                end
            
                % for results, now do the other misfit cost function calculation
                currentresults_cf                = results.misfit_costfunction;
                currentresults_cf                = currentresults_cf(k,:);
                currentresultsddata_cf           = results_temp.misfit_costfunction;
                currentresultsddata_cf(k,:)      = currentresults_cf(:,kpindx);
                results_temp.misfit_costfunction = currentresultsddata_cf;
                if saveuncutdata == 1
                   currentresultsddata_cf                     = results_temp.all_uncut.misfit_costfunction;
                   currentresultsddata_cf(k,:)                = currentresults_cf;
                   results_temp.all_uncut.misfit_costfunction = currentresultsddata_cf;
                end
             
                successfullsimulations_temp(k,1)  = sum(successfullsimulations(k,kpindx),2);
                badmbcount_temp(k,1)              = sum(badmbcount(k,kpindx),2);
                iterationcounts_temp(k,1)         = sum(iterationcounts(k,kpindx),2);
           end    

           % now re-set the matricies that were culled in the loop on k
           EMSucessIterMatrixAll      = EMSucessIterMatrixAll_temp;      clear EMSucessIterMatrixAll_temp;
           EMSucessIterMatrixAll_frac = EMSucessIterMatrixAll_frac_temp; clear EMSucessIterMatrixAll_frac_temp;
           InvRivDataAgg              = InvRivDataAgg_temp;              clear InvRivDataAgg_temp;
           results                    = results_temp;                    clear results_temp;
           successfullsimulations     = successfullsimulations_temp;     clear successfullsimulations_temp;
           badmbcount                 = badmbcount_temp;                 clear badmbcount_temp;
           iterationcounts            = iterationcounts_temp;            clear iterationcounts_temp;
    end

elseif isequal(IterateOver,'Samples')
       remaininginstances = maxsuccess;
       % remaininginstances records how many successful simulations exist
       % by using remaininginstances in the subsequent calculations,
       % MEANDIR can use the same functions when IterateOver equals 'Samples' or 'End-members'

end % end of primary if statement