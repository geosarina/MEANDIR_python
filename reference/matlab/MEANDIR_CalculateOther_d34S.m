% this script calculates an expected value for the d34S of other SO4 
clear notpyrite_list; clear pyrite_list;

fnobs = fieldnames(river.observations);        % get the fieldnames of river observations
if sum(ismember(fnobs,'d34S_conc'))>0          % only run this code if the river data contains a d34S variable
if sum(~isnan(river.observations.d34S_conc))>0 % only run this code if the river data contains a d34S variable
if sum(ismember(ObsList,'d34S'))==0            % only run this code if the inversion does not contain d34S 
if EMHaved34S == 1                             % only run this code is the end-members have d34S values

    % (1) define the list of end-members where SO4 contributions should count as other SO4 or as SO4 from a known end-member (such as precipitation) 
    % if SO4 is included in the inversion, specify using the user-input Z_Numerator_EMList which end-members should count as known
    % When SO4 is not included in the inversion (as is normal when using SO4 excess) call all end-members known
    if     sum(ismember(ObsList,'SO4'))>0
           knownSO4_list = EMList0(~ismember(EMList0,Z_Numerator_EMList));
           otherSO4_list = EMList0(ismember(EMList0,Z_Numerator_EMList));
    elseif sum(ismember(ObsList,'SO4'))==0
           knownSO4_list = EMList0;
           otherSO4_list = {};
    end

    % (2) calculate the SO4 contribution of "other". do this by iterating over the entries of
    % otherSO4_list and adding together their fractional contributions to SO4 in each simulation
    otherSO4_frac = 0;
    for jj=1:length(otherSO4_list)
        if   exist('results','var'); active_otherSO4frac = reshape(results.SO4(ismember(EMList0,otherSO4_list{jj}),:,:),remaininginstances,s);            % makes a matrix of size (remaininginstances x s)
        else;                        active_otherSO4frac = reshape(river.fraction.all.SO4(ismember(EMList0,otherSO4_list{jj}),:,:),remaininginstances,s); % makes a matrix of size (remaininginstances x s)
        end
        otherSO4_frac = otherSO4_frac + active_otherSO4frac;
    end

    % set up an iteration of length n for the number of times to 
    % calculate the d34S of 'other' SO4
    n = 1000; % calculate other d34S

    % (3) isolate the end-members that contribute SO4   
    knownSO4_fracxd34S = NaN(remaininginstances,s,n); % size (remaininginstances x s x n)
    knownSO4_frac      = NaN(remaininginstances,s,n); % size (remaininginstances x s x n)

    % iterate over n, which is the number of times to calculate the d34S of other SO4
    for kk=1:n
        temp_knownSO4_fracxd34S = 0; % reset the temporary product of known SO4 fraction and d34S
        temp_knownSO4_frac      = 0; % reset the temporary product of known SO4 fraction

        % iterate over end-members in knownSO4_list to calculate the known contributions of
        % SO4 to river, as well as that isotopic impact
        for jj=1:length(knownSO4_list)    

            % pull the d34S distribution for this knwon SO4 end-member, then find a value for its d34S by pulling from its distribution
            disttype = eval(sprintf('EM.%s.disttype.d34S',EMdatasource));
            if     isequal(disttype,'NOR')
                   active_mean = eval(sprintf('EM.%s.%s.NOR_Mend34S%s',EMdatasource,knownSO4_list{jj},NormalizationType));
                   active_unc  = eval(sprintf('EM.%s.%s.NOR_Sigd34S%s',EMdatasource,knownSO4_list{jj},NormalizationType));
                   active_d34S = active_mean + randn(remaininginstances,s)*active_unc;
            elseif isequal(disttype,'UNI')
                   active1  = eval(sprintf('EM.%s.%s.UNI_Mind34S%s',EMdatasource,knownSO4_list{jj},NormalizationType));
                   active2  = eval(sprintf('EM.%s.%s.UNI_Maxd34S%s',EMdatasource,knownSO4_list{jj},NormalizationType));
                   active_d34S = active1 + rand(remaininginstances,s)*(active2-active1);
            end      

            % pull the fractional SO4 contribution for this known SO4 end-member
            if   exist('results','var');  active_knownSO4frac = reshape(results.SO4(ismember(EMList0,knownSO4_list{jj}),:,:),remaininginstances,s);            % makes a matrix of size (remaininginstances x s)
            else;                         active_knownSO4frac = reshape(river.fraction.all.SO4(ismember(EMList0,knownSO4_list{jj}),:,:),remaininginstances,s); % makes a matrix of size (remaininginstances x s)
            end  

            % adjust temp_knownSO4_frac andtemp_knownSO4_frac
            temp_knownSO4_fracxd34S = temp_knownSO4_fracxd34S + active_d34S.*active_knownSO4frac; % the evolving product of SO4 fractional contribution and d34S (remaininginstances x s)
            temp_knownSO4_frac      = temp_knownSO4_frac      + active_knownSO4frac;              % the evolving value of  SO4 fractional contribution (remaininginstances x s)
        end     

    knownSO4_fracxd34S(1:remaininginstances,1:s,kk) = temp_knownSO4_fracxd34S; % the saved product of SO4 fractional contribution and d34S (remaininginstances x s x n)
    knownSO4_frac(1:remaininginstances,1:s,kk)      = temp_knownSO4_frac;      % the saved product of SO4 fractional contribution (remaininginstances x s x n)
    end

    % (4) calculate the d34S of the other/FeS2/SO4 excess and remove impossible values
    otherSO4_d34S = NaN(remaininginstances,s,n); % preallocate this variable, size (remaininginstances x s x n)
    if     sum(ismember(ObsList,'SO4'))>0
           % (4.1a) first consider the case where SO4 was included in the inversion
           % this would be typical for putting FeS2 into end-members, or having it be a stand-alone end-member
           for kk=1:n
               otherSO4_d34S(:,:,kk) = ((otherSO4_frac+knownSO4_frac(:,:,kk)).*repmat(river.observations.d34S_conc',remaininginstances,1) - knownSO4_fracxd34S(:,:,kk))./(otherSO4_frac);
               % double-check test for mass balance in the sum of otherSO4_frac and knownSO4_frac
               mnmb = 100*min(min(otherSO4_frac+knownSO4_frac(:,:,kk)));
               mxmb = 100*max(max(otherSO4_frac+knownSO4_frac(:,:,kk)));
               if mnmb<ErrorCutMinMB(ismember(ObsList,'SO4')) | mxmb>ErrorCutMaxMB(ismember(ObsList,'SO4'))
                  disp('Warning! Mass balance error for SO4 in "MEANDR_Calculateother_d34S"'); 
               end           
           end 
    elseif sum(ismember(ObsList,'SO4'))==0           
           % (4.1b) next consider the case where SO4 was not included in the
           % inversion. this would be typical for calculating FeS2 inputs as SO4
           % excess. in this case, we assume that the total amount of fractional
           % SO4 is 1, unlike the case above where we calculate FeS2 using the
           % model's reconstructed amounts of SO4.
           for kk=1:n
               otherSO4_d34S(:,:,kk) = (repmat(river.observations.d34S_conc',remaininginstances,1) - knownSO4_fracxd34S(:,:,kk))./(1 - knownSO4_frac(:,:,kk));
           end
    end
   
    % (4.2) only keep the 'other' data when the other end-member contributes above a critical threshold (here that is taken as 1%)
    % the code also removes simulations where d34S is less than 1000, as this would correspond to negative 34S atoms
    criticalother     = 0.01;          % the critical 'other SO4' contribution that must be met for calculating d34S
    otherSO4_d34S_cut = otherSO4_d34S; % define a new variable to contain the cut data
    for kk=1:n
        activeotherd34Ssheet  = otherSO4_d34S_cut(:,:,kk);       % find the active sheet of otherSO4_d34S_cut
        activeotherd34Ssheet(otherSO4_frac<criticalother) = NaN; % remove simulations where otherSO4_frac was below the threshold
        activeotherd34Ssheet(activeotherd34Ssheet<-1000)  = NaN; % remove simulations where d34S was calculated to be < - 1000
        otherSO4_d34S_cut(:,:,kk) = activeotherd34Ssheet;        % reset the sheet within otherSO4_d34S_cut
    end  

    % (5) reshape and save the data    
    mean_d34S     = NaN(s,1); mean_d34S_cut     = NaN(s,1); % preallocate mean vectors
    median_d34S   = NaN(s,1); median_d34S_cut   = NaN(s,1); % preallocate median vectors
    other_d34S_05 = NaN(s,1); other_d34S_05_cut = NaN(s,1); % preallocate pct 05 vectors
    other_d34S_25 = NaN(s,1); other_d34S_25_cut = NaN(s,1); % preallocate pct 25 vectors
    other_d34S_75 = NaN(s,1); other_d34S_75_cut = NaN(s,1); % preallocate pct 75 vectors
    other_d34S_95 = NaN(s,1); other_d34S_95_cut = NaN(s,1); % preallocate pct 95 vectors

    for zz = 1:s  
        mean_d34S(zz,1)         = nanmean(reshape(otherSO4_d34S(:,zz,:),remaininginstances*n,1));
        median_d34S(zz,1)       = nanmedian(reshape(otherSO4_d34S(:,zz,:),remaininginstances*n,1));
        other_d34S_05(zz,1)     = prctile(reshape(otherSO4_d34S(:,zz,:),remaininginstances*n,1),05);
        other_d34S_25(zz,1)     = prctile(reshape(otherSO4_d34S(:,zz,:),remaininginstances*n,1),25);
        other_d34S_75(zz,1)     = prctile(reshape(otherSO4_d34S(:,zz,:),remaininginstances*n,1),75);
        other_d34S_95(zz,1)     = prctile(reshape(otherSO4_d34S(:,zz,:),remaininginstances*n,1),95);

        mean_d34S_cut(zz,1)     = nanmean(reshape(otherSO4_d34S_cut(:,zz,:),remaininginstances*n,1));
        median_d34S_cut(zz,1)   = nanmedian(reshape(otherSO4_d34S_cut(:,zz,:),remaininginstances*n,1));
        other_d34S_05_cut(zz,1) = prctile(reshape(otherSO4_d34S_cut(:,zz,:),remaininginstances*n,1),05);
        other_d34S_25_cut(zz,1) = prctile(reshape(otherSO4_d34S_cut(:,zz,:),remaininginstances*n,1),25);
        other_d34S_75_cut(zz,1) = prctile(reshape(otherSO4_d34S_cut(:,zz,:),remaininginstances*n,1),75);
        other_d34S_95_cut(zz,1) = prctile(reshape(otherSO4_d34S_cut(:,zz,:),remaininginstances*n,1),95);
    end

    % (6) reshape and save a subset of the data
    all_mean_d34S     = NaN(remaininginstances, s); all_mean_d34S_cut     = NaN(remaininginstances, s); % preallocate mean matrix
    all_median_d34S   = NaN(remaininginstances, s); all_median_d34S_cut   = NaN(remaininginstances, s); % preallocate median matrix
    all_other_d34S_05 = NaN(remaininginstances, s); all_other_d34S_05_cut = NaN(remaininginstances, s); % preallocate pct 05 matrix
    all_other_d34S_25 = NaN(remaininginstances, s); all_other_d34S_25_cut = NaN(remaininginstances, s); % preallocate pct 25 matrix
    all_other_d34S_75 = NaN(remaininginstances, s); all_other_d34S_75_cut = NaN(remaininginstances, s); % preallocate pct 75 matrix
    all_other_d34S_95 = NaN(remaininginstances, s); all_other_d34S_95_cut = NaN(remaininginstances, s); % preallocate pct 95 matrix
    alldist     = {};
    alldist_cut = {};

    for xx = 1:remaininginstances
    for yy = 1:s 

        all_mean_d34S(xx,yy)         = nanmean(otherSO4_d34S(xx,yy,:));
        all_median_d34S(xx,yy)       = nanmedian(otherSO4_d34S(xx,yy,:));
        all_other_d34S_05(xx,yy)     = prctile(otherSO4_d34S(xx,yy,:),05);
        all_other_d34S_25(xx,yy)     = prctile(otherSO4_d34S(xx,yy,:),25);
        all_other_d34S_75(xx,yy)     = prctile(otherSO4_d34S(xx,yy,:),75);
        all_other_d34S_95(xx,yy)     = prctile(otherSO4_d34S(xx,yy,:),95);

        all_mean_d34S_cut(xx,yy)     = nanmean(otherSO4_d34S_cut(xx,yy,:));
        all_median_d34S_cut(xx,yy)   = nanmedian(otherSO4_d34S_cut(xx,yy,:));
        all_other_d34S_05_cut(xx,yy) = prctile(otherSO4_d34S_cut(xx,yy,:),05);
        all_other_d34S_25_cut(xx,yy) = prctile(otherSO4_d34S_cut(xx,yy,:),25);
        all_other_d34S_75_cut(xx,yy) = prctile(otherSO4_d34S_cut(xx,yy,:),75);
        all_other_d34S_95_cut(xx,yy) = prctile(otherSO4_d34S_cut(xx,yy,:),95);

        % if requested, save the full distribution of FeS2 d34S values
        if RecordFullFeS2Distribution == 1
            datfit = reshape(otherSO4_d34S(xx,yy,:),n,1); 
            datfit(datfit == Inf) = NaN;
            datfit(datfit ==-Inf) = NaN;
            if sum(~isnan(datfit))>100
                alldist{xx,yy} = fitdist(datfit,'kernel');
            else
                alldist{xx,yy} = NaN;
            end

            datfit_cut = reshape(otherSO4_d34S_cut(xx,yy,:),n,1); 
            datfit_cut(datfit_cut == Inf) = NaN;
            datfit_cut(datfit_cut ==-Inf) = NaN;
            if sum(~isnan(datfit_cut))>100
                alldist_cut{xx,yy} = fitdist(datfit_cut,'kernel');
            else
                alldist_cut{xx,yy} = NaN;
            end

        else
            alldist{xx,yy}     = NaN;
            alldist_cut{xx,yy} = NaN;
        end
    end
    end

    % (7) save FeS2 d34S into the main output variable 
    river.calculated_other_d34S.uncut.mean       = mean_d34S;              clear mean_d34S;             % uncut mean vector
    river.calculated_other_d34S.uncut.median     = median_d34S;            clear median_d34S;           % uncut median vector
    river.calculated_other_d34S.uncut.pct05      = other_d34S_05;          clear other_d34S_05;         % uncut pct 05 vector
    river.calculated_other_d34S.uncut.pct25      = other_d34S_25;          clear other_d34S_25;         % uncut pct 25 vector
    river.calculated_other_d34S.uncut.pct75      = other_d34S_75;          clear other_d34S_75;         % uncut pct 75 vector
    river.calculated_other_d34S.uncut.pct95      = other_d34S_95;          clear other_d34S_95;         % uncut pct 95 vector

    river.calculated_other_d34S.uncut.all.mean   = all_mean_d34S;          clear all_mean_d34S;         % uncut mean matrix
    river.calculated_other_d34S.uncut.all.median = all_median_d34S;        clear all_median_d34S;       % uncut mean matrix
    river.calculated_other_d34S.uncut.all.pct05  = all_other_d34S_05;      clear all_other_d34S_05;     % uncut pct 05 matrix
    river.calculated_other_d34S.uncut.all.pct25  = all_other_d34S_25;      clear all_other_d34S_25;     % uncut pct 25 matrix
    river.calculated_other_d34S.uncut.all.pct75  = all_other_d34S_75;      clear all_other_d34S_75;     % uncut pct 75 matrix
    river.calculated_other_d34S.uncut.all.pct95  = all_other_d34S_95;      clear all_other_d34S_95;     % uncut pct 95 matrix

    river.calculated_other_d34S.cut.mean         = mean_d34S_cut;          clear mean_d34S_cut;         % cut mean vector
    river.calculated_other_d34S.cut.median       = median_d34S_cut;        clear median_d34S_cut;       % cut median vector
    river.calculated_other_d34S.cut.pct05        = other_d34S_05_cut;      clear other_d34S_05_cut;     % cut pct 05 vector
    river.calculated_other_d34S.cut.pct25        = other_d34S_25_cut;      clear other_d34S_25_cut;     % cut pct 25 vector
    river.calculated_other_d34S.cut.pct75        = other_d34S_75_cut;      clear other_d34S_75_cut;     % cut pct 75 vector
    river.calculated_other_d34S.cut.pct95        = other_d34S_95_cut;      clear other_d34S_95_cut;     % cut pct 95 vector

    river.calculated_other_d34S.cut.all.mean     = all_mean_d34S_cut;      clear all_mean_d34S_cut;     % cut mean matrix
    river.calculated_other_d34S.cut.all.median   = all_median_d34S_cut;    clear all_median_d34S_cut;   % cut median matrix
    river.calculated_other_d34S.cut.all.pct05    = all_other_d34S_05_cut;  clear all_other_d34S_05_cut; % cut pct 05 matrix
    river.calculated_other_d34S.cut.all.pct25    = all_other_d34S_25_cut;  clear all_other_d34S_25_cut; % cut pct 25 matrix
    river.calculated_other_d34S.cut.all.pct75    = all_other_d34S_75_cut;  clear all_other_d34S_75_cut; % cut pct 75 matrix
    river.calculated_other_d34S.cut.all.pct95    = all_other_d34S_95_cut;  clear all_other_d34S_95_cut; % cut pct 95 matrix

    if exist('results','var')
       disp('calculated the d34S of "other" SO4, saved results in field "calculated_other_d34S"');
    end

    % (8) if requested, save the full distributions
    if RecordFullFeS2Distribution == 1
       river.calculated_other_d34S.uncut.all.dist = alldist; 
       river.calculated_other_d34S.cut.all.dist   = alldist_cut; 
       clear alldist;
    end    

end % if EMHaved34S == 1 
end % if sum(ismember(ObsList,'d34S'))==0
end % if sum(~isnan(river.observations.d34S_conc))>0
end % if sum(ismember(fnobs,'d34S_conc'))>0

clear other_d34S;
clear store_den_nonpyrite; 
clear store_num_nonpyrite;