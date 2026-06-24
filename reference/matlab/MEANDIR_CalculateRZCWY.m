% this script calculates net and gross R, Z, C, W, and Y for the inversion results
% this script is only executed if the user-defined variable CalculateRZCWY = 1

if  CalculateRZCWY == 1                 % only proceed if the user requested RZCWY values
    if sum(successfullsimulations>0)>0  % only proceed if there was at least 1 successful result

        % (0) define variables to hold results
        R_gross_median_unscaled = NaN(s,1); R_gross_pct05_unscaled = NaN(s,1); R_gross_pct25_unscaled = NaN(s,1); R_gross_pct75_unscaled = NaN(s,1); R_gross_pct95_unscaled = NaN(s,1); % declare R matricies
        Z_gross_median_unscaled = NaN(s,1); Z_gross_pct05_unscaled = NaN(s,1); Z_gross_pct25_unscaled = NaN(s,1); Z_gross_pct75_unscaled = NaN(s,1); Z_gross_pct95_unscaled = NaN(s,1); % declare Z matricies
        C_gross_median_unscaled = NaN(s,1); C_gross_pct05_unscaled = NaN(s,1); C_gross_pct25_unscaled = NaN(s,1); C_gross_pct75_unscaled = NaN(s,1); C_gross_pct95_unscaled = NaN(s,1); % declare C matricies
        Y_gross_median_unscaled = NaN(s,1); Y_gross_pct05_unscaled = NaN(s,1); Y_gross_pct25_unscaled = NaN(s,1); Y_gross_pct75_unscaled = NaN(s,1); Y_gross_pct95_unscaled = NaN(s,1); % declare Y matricies
        W_gross_median_unscaled = NaN(s,1); W_gross_pct05_unscaled = NaN(s,1); W_gross_pct25_unscaled = NaN(s,1); W_gross_pct75_unscaled = NaN(s,1); W_gross_pct95_unscaled = NaN(s,1); % declare W matricies

        R_gross_median_scaled   = NaN(s,1); R_gross_pct05_scaled   = NaN(s,1); R_gross_pct25_scaled   = NaN(s,1); R_gross_pct75_scaled   = NaN(s,1); R_gross_pct95_scaled   = NaN(s,1); % declare R matricies
        Z_gross_median_scaled   = NaN(s,1); Z_gross_pct05_scaled   = NaN(s,1); Z_gross_pct25_scaled   = NaN(s,1); Z_gross_pct75_scaled   = NaN(s,1); Z_gross_pct95_scaled   = NaN(s,1); % declare Z matricies
        C_gross_median_scaled   = NaN(s,1); C_gross_pct05_scaled   = NaN(s,1); C_gross_pct25_scaled   = NaN(s,1); C_gross_pct75_scaled   = NaN(s,1); C_gross_pct95_scaled   = NaN(s,1); % declare C matricies
        Y_gross_median_scaled   = NaN(s,1); Y_gross_pct05_scaled   = NaN(s,1); Y_gross_pct25_scaled   = NaN(s,1); Y_gross_pct75_scaled   = NaN(s,1); Y_gross_pct95_scaled   = NaN(s,1); % declare Y matricies
        W_gross_median_scaled   = NaN(s,1); W_gross_pct05_scaled   = NaN(s,1); W_gross_pct25_scaled   = NaN(s,1); W_gross_pct75_scaled   = NaN(s,1); W_gross_pct95_scaled   = NaN(s,1); % declare W matricies

        R_net_median_unscaled   = NaN(s,1); R_net_pct05_unscaled   = NaN(s,1); R_net_pct25_unscaled   = NaN(s,1); R_net_pct75_unscaled   = NaN(s,1); R_net_pct95_unscaled   = NaN(s,1); % declare R matricies
        Z_net_median_unscaled   = NaN(s,1); Z_net_pct05_unscaled   = NaN(s,1); Z_net_pct25_unscaled   = NaN(s,1); Z_net_pct75_unscaled   = NaN(s,1); Z_net_pct95_unscaled   = NaN(s,1); % declare Z matricies
        C_net_median_unscaled   = NaN(s,1); C_net_pct05_unscaled   = NaN(s,1); C_net_pct25_unscaled   = NaN(s,1); C_net_pct75_unscaled   = NaN(s,1); C_net_pct95_unscaled   = NaN(s,1); % declare C matricies
        Y_net_median_unscaled   = NaN(s,1); Y_net_pct05_unscaled   = NaN(s,1); Y_net_pct25_unscaled   = NaN(s,1); Y_net_pct75_unscaled   = NaN(s,1); Y_net_pct95_unscaled   = NaN(s,1); % declare Y matricies
        W_net_median_unscaled   = NaN(s,1); W_net_pct05_unscaled   = NaN(s,1); W_net_pct25_unscaled   = NaN(s,1); W_net_pct75_unscaled   = NaN(s,1); W_net_pct95_unscaled   = NaN(s,1); % declare W matricies

        R_net_median_scaled     = NaN(s,1); R_net_pct05_scaled     = NaN(s,1); R_net_pct25_scaled     = NaN(s,1); R_net_pct75_scaled     = NaN(s,1); R_net_pct95_scaled     = NaN(s,1); % declare R matricies
        Z_net_median_scaled     = NaN(s,1); Z_net_pct05_scaled     = NaN(s,1); Z_net_pct25_scaled     = NaN(s,1); Z_net_pct75_scaled     = NaN(s,1); Z_net_pct95_scaled     = NaN(s,1); % declare Z matricies
        C_net_median_scaled     = NaN(s,1); C_net_pct05_scaled     = NaN(s,1); C_net_pct25_scaled     = NaN(s,1); C_net_pct75_scaled     = NaN(s,1); C_net_pct95_scaled     = NaN(s,1); % declare C matricies
        Y_net_median_scaled     = NaN(s,1); Y_net_pct05_scaled     = NaN(s,1); Y_net_pct25_scaled     = NaN(s,1); Y_net_pct75_scaled     = NaN(s,1); Y_net_pct95_scaled     = NaN(s,1); % declare Y matricies
        W_net_median_scaled     = NaN(s,1); W_net_pct05_scaled     = NaN(s,1); W_net_pct25_scaled     = NaN(s,1); W_net_pct75_scaled     = NaN(s,1); W_net_pct95_scaled     = NaN(s,1); % declare W matricies

        agg_RZC_denominator_equi_gross_unscaled = NaN(s,remaininginstances); agg_R_numerator_equi_gross_unscaled = NaN(s,remaininginstances); SO4fromFeS2_gross_unscaled = NaN(s,remaininginstances);
        agg_RZC_denominator_equi_gross_scaled   = NaN(s,remaininginstances); agg_R_numerator_equi_gross_scaled   = NaN(s,remaininginstances); SO4fromFeS2_gross_scaled   = NaN(s,remaininginstances);
        agg_RZC_denominator_equi_net_unscaled   = NaN(s,remaininginstances); agg_R_numerator_equi_net_unscaled   = NaN(s,remaininginstances); SO4fromFeS2_net_unscaled   = NaN(s,remaininginstances);
        agg_RZC_denominator_equi_net_scaled     = NaN(s,remaininginstances); agg_R_numerator_equi_net_scaled     = NaN(s,remaininginstances); SO4fromFeS2_net_scaled     = NaN(s,remaininginstances);

        if isequal(carbonisotopematch,'DIC'); DICfromCorg_gross_unscaled  = NaN(s,remaininginstances); DICfromCorg_gross_scaled  = NaN(s,remaininginstances); DICfromCorg_net_unscaled  = NaN(s,remaininginstances); DICfromCorg_net_scaled  = NaN(s,remaininginstances);
        else;                                 HCO3fromCorg_gross_unscaled = NaN(s,remaininginstances); HCO3fromCorg_gross_scaled = NaN(s,remaininginstances); HCO3fromCorg_net_unscaled = NaN(s,remaininginstances); HCO3fromCorg_net_scaled = NaN(s,remaininginstances);
        end

        % find the index of end-members included in the denominator of R, Z, and C
        RZC_Denominator_EMindx = find(ismember(EMList0,RZC_Denominator_EMList));
        disp('beginning to calculate R, Z, C, W, Y');

        % define conc2equi structure
        conc2equi = MEANDIR_DefineConc2Equi;  

        % iterate over the number of samples   
        for i=1:s 

            % (1) calculate the denominator of R, Z, and C
            RZC_denominator_equi_gross_unscaled = 0; % variable to accumulate gross contributions
            RZC_denominator_equi_gross_scaled   = 0; % variable to accumulate gross contributions
            RZC_denominator_equi_net_unscaled   = 0; % variable to accumulate net contributions
            RZC_denominator_equi_net_scaled     = 0; % variable to accumulate net contributions

            % iterate over the ions of interest
            for j=1:length(RZC_Denominator_IonList)

                % (1.1) get the results for this ion and the river data for this ion
                activeresults = eval(sprintf('results.%s',RZC_Denominator_IonList{j}));
                riverdata     = eval(sprintf('river.InvertedRiverValues.%s.all',RZC_Denominator_IonList{j}));
                rivervalues   = riverdata(i,:); % get the river data just for this variable from just this sample

                % (1.2) set the contributions for this particular ion to begin at zero
                contributionfromthision_gross_unscaled = 0;
                contributionfromthision_gross_scaled   = 0;
                contributionfromthision_net_unscaled   = 0;
                contributionfromthision_net_scaled     = 0;

                % (1.3) get the scaling factor for equivalents/concentration
                if     isequal(EMUnits,'conc'); conc2equimult = eval(sprintf('conc2equi.%s',RZC_Denominator_IonList{j}));
                elseif isequal(EMUnits,'equi'); conc2equimult = 1;
                end

                % (1.4) now iterate over each end-member
                fraction_allsinks               = sum(activeresults(ismember(EMList0,EMsinks),:,i),1);
                fraction_allsources             = sum(activeresults(ismember(EMList0,EMsources),:,i),1);
                fraction_allnet                 = sum(activeresults(:,:,i),1);
                for k = 1:length(RZC_Denominator_EMList)
                    fraction_activeionfromactiveEM         = activeresults(RZC_Denominator_EMindx(k),:,i);
                    contributionfromthision_gross_unscaled = contributionfromthision_gross_unscaled + conc2equimult*(rivervalues.*fraction_activeionfromactiveEM);
                    contributionfromthision_gross_scaled   = contributionfromthision_gross_scaled   + conc2equimult*(rivervalues.*fraction_activeionfromactiveEM./fraction_allnet);
                    contributionfromthision_net_unscaled   = contributionfromthision_net_unscaled   + conc2equimult*(rivervalues.*fraction_activeionfromactiveEM.*(1+fraction_allsinks./fraction_allsources));
                    contributionfromthision_net_scaled     = contributionfromthision_net_scaled     + conc2equimult*(rivervalues.*fraction_activeionfromactiveEM./fraction_allnet.*(1+fraction_allsinks./fraction_allsources));
                end
                RZC_denominator_equi_gross_unscaled = RZC_denominator_equi_gross_unscaled + contributionfromthision_gross_unscaled;
                RZC_denominator_equi_gross_scaled   = RZC_denominator_equi_gross_scaled   + contributionfromthision_gross_scaled;
                RZC_denominator_equi_net_unscaled   = RZC_denominator_equi_net_unscaled   + contributionfromthision_net_unscaled;
                RZC_denominator_equi_net_scaled     = RZC_denominator_equi_net_scaled     + contributionfromthision_net_scaled;
            end

            % store the result for this sample
            agg_RZC_denominator_equi_gross_unscaled(i,:)   = RZC_denominator_equi_gross_unscaled;
            agg_RZC_denominator_equi_gross_scaled(i,:)     = RZC_denominator_equi_gross_scaled;
            agg_RZC_denominator_equi_net_unscaled(i,:)     = RZC_denominator_equi_net_unscaled;
            agg_RZC_denominator_equi_net_scaled(i,:)       = RZC_denominator_equi_net_scaled;

            % (2) calculate the numerator of R
            R_numerator_equi_gross_unscaled = 0; % variable to accumulate gross contributions
            R_numerator_equi_net_unscaled   = 0; % variable to accumulate net contributions
            R_numerator_equi_gross_scaled   = 0; % variable to accumulate gross contributions
            R_numerator_equi_net_scaled     = 0; % variable to accumulate net contributions

            R_Numerator_EMindx = find(ismember(EMList0,R_Numerator_EMList));
            for j=1:length(R_Numerator_IonList)  % iterate over the ions of interest

                % (2.1) get the results for this ion and the river data
                activeresults = eval(sprintf('results.%s',R_Numerator_IonList{j}));
                riverdata     = eval(sprintf('river.InvertedRiverValues.%s.all',R_Numerator_IonList{j}));
                rivervalues   = riverdata(i,:);

                % (2.2) set the contributions for this particular ion to begin at zero
                contributionfromthision_gross_unscaled = 0;
                contributionfromthision_gross_scaled   = 0;
                contributionfromthision_net_unscaled   = 0;
                contributionfromthision_net_scaled     = 0;

                % (2.3) get the scaling factor for equivalents/concentration
                if     isequal(EMUnits,'conc'); conc2equimult = eval(sprintf('conc2equi.%s',R_Numerator_IonList{j}));
                elseif isequal(EMUnits,'equi'); conc2equimult = 1;
                end

                % (2.4) now iterate over each end-member
                fraction_allsinks               = sum(activeresults(ismember(EMList0,EMsinks),:,i),1);
                fraction_allsources             = sum(activeresults(ismember(EMList0,EMsources),:,i),1);
                fraction_allnet                 = sum(activeresults(:,:,i),1);
                for k = 1:length(R_Numerator_EMList)
                    fraction_activeionfromactiveEM         = activeresults(R_Numerator_EMindx(k),:,i);
                    contributionfromthision_gross_unscaled = contributionfromthision_gross_unscaled + conc2equimult*(rivervalues.*fraction_activeionfromactiveEM);
                    contributionfromthision_gross_scaled   = contributionfromthision_gross_scaled   + conc2equimult*(rivervalues.*fraction_activeionfromactiveEM./fraction_allnet);
                    contributionfromthision_net_unscaled   = contributionfromthision_net_unscaled   + conc2equimult*(rivervalues.*fraction_activeionfromactiveEM.*(1+fraction_allsinks./fraction_allsources));
                    contributionfromthision_net_scaled     = contributionfromthision_net_scaled     + conc2equimult*(rivervalues.*fraction_activeionfromactiveEM./fraction_allnet.*(1+fraction_allsinks./fraction_allsources));
                end
                R_numerator_equi_gross_unscaled = R_numerator_equi_gross_unscaled + contributionfromthision_gross_unscaled;
                R_numerator_equi_gross_scaled   = R_numerator_equi_gross_scaled   + contributionfromthision_gross_scaled;
                R_numerator_equi_net_unscaled   = R_numerator_equi_net_unscaled   + contributionfromthision_net_unscaled;
                R_numerator_equi_net_scaled     = R_numerator_equi_net_scaled     + contributionfromthision_net_scaled;
            end

            % store the result for this sample
            agg_R_numerator_equi_gross_unscaled(i,:) = R_numerator_equi_gross_unscaled;
            agg_R_numerator_equi_gross_scaled(i,:)   = R_numerator_equi_gross_scaled;
            agg_R_numerator_equi_net_unscaled(i,:)   = R_numerator_equi_net_unscaled;
            agg_R_numerator_equi_net_scaled(i,:)     = R_numerator_equi_net_scaled;


            % (3) solve for Z, gross and net
            % (3.1) here the end-members either contain pyrite or are pyrite
            if  ZfromEM==1
                   InvRivDatSO4                        = river.InvertedRiverValues.SO4.all;
                   Z_Numerator_EMindx                  = find(ismember(EMList0,Z_Numerator_EMList));
                   Z_numerator_equi_gross_unscaled     = 0;
                   Z_numerator_equi_gross_scaled       = 0;
                   Z_numerator_equi_net_unscaled       = 0;
                   Z_numerator_equi_net_scaled         = 0;

                   fraction_allsinks                   = sum(results.SO4(ismember(EMList0,EMsinks),:,i),1);
                   fraction_allsources                 = sum(results.SO4(ismember(EMList0,EMsources),:,i),1);
                   fraction_allnet                     = sum(results.SO4(:,:,i),1);
                   for k = 1:length(Z_Numerator_EMList)
                       fraction_SO4fromactiveEM        = results.SO4(Z_Numerator_EMindx(k),:,i);
                       Z_numerator_equi_gross_unscaled = Z_numerator_equi_gross_unscaled +  (InvRivDatSO4(i,:).*fraction_SO4fromactiveEM);
                       Z_numerator_equi_gross_scaled   = Z_numerator_equi_gross_scaled   +  (InvRivDatSO4(i,:).*fraction_SO4fromactiveEM./fraction_allnet);
                       Z_numerator_equi_net_unscaled   = Z_numerator_equi_net_unscaled   +  (InvRivDatSO4(i,:).*fraction_SO4fromactiveEM.*(1+fraction_allsinks./fraction_allsources));
                       Z_numerator_equi_net_scaled     = Z_numerator_equi_net_scaled     +  (InvRivDatSO4(i,:).*fraction_SO4fromactiveEM./fraction_allnet.*(1+fraction_allsinks./fraction_allsources));
                   end
                   SO4fromFeS2_gross_unscaled(i,:)     = Z_numerator_equi_gross_unscaled;
                   SO4fromFeS2_gross_scaled(i,:)       = Z_numerator_equi_gross_scaled;
                   SO4fromFeS2_net_unscaled(i,:)       = Z_numerator_equi_net_unscaled;
                   SO4fromFeS2_net_scaled(i,:)         = Z_numerator_equi_net_scaled;

            % (3.2) here the SO4 from pyrite weathering is excess SO4
            elseif ZfromSO4excess==1
                   InvRivDatSO4                        =  repmat(river.model_variable.SO4,1,remaininginstances);
                   SO4fromFeS2_gross_unscaled(i,:)     =  river.model_variable.SO4(i).*(1-sum(results.SO4(ismember(EMList0,EMsources),:,i),1));
                   SO4fromFeS2_gross_scaled(i,:)       =  NaN(size(SO4fromFeS2_gross_unscaled(i,:)));
                   SO4fromFeS2_net_unscaled(i,:)       =  ExcessSO4(i,:);
                   SO4fromFeS2_net_scaled(i,:)         =  NaN(size(ExcessSO4(i,:)));

            % (3.3) here the SO4 from weathering is SO4 in the river
            elseif ZfromriverSO4==1
                   InvRivDatSO4                        = river.InvertedRiverValues.SO4.all;
                   SO4fromFeS2_gross_unscaled(i,:)     = InvRivDatSO4(i,:);
                   SO4fromFeS2_gross_scaled(i,:)       = NaN(size(InvRivDatSO4(i,:)));
                   SO4fromFeS2_net_unscaled(i,:)       = InvRivDatSO4(i,:);
                   SO4fromFeS2_net_scaled(i,:)         = NaN(size(InvRivDatSO4(i,:)));

            % (3.4) in this case, Z is not calculated
            elseif Znotcalculated == 1
                   InvRivDatSO4(i,:)                   = NaN(1,remaininginstances);
            end

            % (6) solve for C - HCO3 from Corg oxidation divided by weathering ions
            if isequal(carbonisotopematch,'DIC')
                if    ~isempty(C_Numerator_EMList)
                      InvRivDatDIC                        = river.InvertedRiverValues.DIC.all;
                      C_Numerator_EMindx                  = find(ismember(EMList0,C_Numerator_EMList));
                      C_numerator_gross_DIC_unscaled      = 0;
                      C_numerator_gross_DIC_scaled        = 0;
                      C_numerator_net_DIC_unscaled        = 0;
                      C_numerator_net_DIC_scaled          = 0;
                      fraction_allsinks                   = sum(results.DIC(ismember(EMList0,EMsinks),:,i),1);
                      fraction_allsources                 = sum(results.DIC(ismember(EMList0,EMsources),:,i),1);
                      fraction_allnet                     = sum(results.DIC(:,:,i),1);
                      for k = 1:length(C_Numerator_EMList)
                          fraction_DICfromactiveEM        = results.DIC(C_Numerator_EMindx(k),:,i);
                          C_numerator_gross_DIC_unscaled  = C_numerator_gross_DIC_unscaled + (InvRivDatDIC(i,:).*fraction_DICfromactiveEM);
                          C_numerator_gross_DIC_scaled    = C_numerator_gross_DIC_scaled   + (InvRivDatDIC(i,:).*fraction_DICfromactiveEM./fraction_allnet);
                          C_numerator_net_DIC_unscaled    = C_numerator_net_DIC_unscaled   + (InvRivDatDIC(i,:).*fraction_DICfromactiveEM.*(1+fraction_allsinks./fraction_allsources));
                          C_numerator_net_DIC_scaled      = C_numerator_net_DIC_scaled     + (InvRivDatDIC(i,:).*fraction_DICfromactiveEM./fraction_allnet.*(1+fraction_allsinks./fraction_allsources));
                      end
                      DICfromCorg_gross_unscaled(i,:)     = C_numerator_gross_DIC_unscaled;
                      DICfromCorg_gross_scaled(i,:)       = C_numerator_gross_DIC_scaled;
                      DICfromCorg_net_unscaled(i,:)       = C_numerator_net_DIC_unscaled; 
                      DICfromCorg_net_scaled(i,:)         = C_numerator_net_DIC_scaled;
                end
            else
                if    ~isempty(C_Numerator_EMList)
                      InvRivDatHCO3                       = river.InvertedRiverValues.HCO3.all;
                      C_Numerator_EMindx                  = find(ismember(EMList0,C_Numerator_EMList));
                      C_numerator_gross_HCO3_unscaled     = 0;
                      C_numerator_gross_HCO3_scaled       = 0;
                      C_numerator_net_HCO3_unscaled       = 0;
                      C_numerator_net_HCO3_scaled         = 0;
                      fraction_allsinks                   = sum(results.HCO3(ismember(EMList0,EMsinks),:,i),1);
                      fraction_allsources                 = sum(results.HCO3(ismember(EMList0,EMsources),:,i),1);
                      fraction_allnet                     = sum(results.HCO3(:,:,i),1);
                      for k = 1:length(C_Numerator_EMList)
                          fraction_HCO3fromactiveEM       = results.HCO3(C_Numerator_EMindx(k),:,i);
                          C_numerator_gross_HCO3_unscaled = C_numerator_gross_HCO3_unscaled + (InvRivDatHCO3(i,:).*fraction_HCO3fromactiveEM);
                          C_numerator_gross_HCO3_scaled   = C_numerator_gross_HCO3_scaled   + (InvRivDatHCO3(i,:).*fraction_HCO3fromactiveEM./fraction_allnet);
                          C_numerator_net_HCO3_unscaled   = C_numerator_net_HCO3_unscaled   + (InvRivDatHCO3(i,:).*fraction_HCO3fromactiveEM.*(1+fraction_allsinks./fraction_allsources));
                          C_numerator_net_HCO3_scaled     = C_numerator_net_HCO3_scaled     + (InvRivDatHCO3(i,:).*fraction_HCO3fromactiveEM./fraction_allnet.*(1+fraction_allsinks./fraction_allsources));
                      end
                      HCO3fromCorg_gross_unscaled(i,:)    = C_numerator_gross_HCO3_unscaled;
                      HCO3fromCorg_gross_scaled(i,:)      = C_numerator_gross_HCO3_scaled;
                      HCO3fromCorg_net_unscaled(i,:)      = C_numerator_net_HCO3_unscaled;
                      HCO3fromCorg_net_scaled(i,:)        = C_numerator_net_HCO3_scaled;
                end  % whether C_Numerator_EMList is empty
            end      % whether carbonisotopematch equals DIC

        end          % end main loop (for i==1:s)

        % now actually solve for the variables of interest
        % solve for R
        R_gross_unscaled = agg_R_numerator_equi_gross_unscaled./agg_RZC_denominator_equi_gross_unscaled;
        R_gross_scaled   = agg_R_numerator_equi_gross_scaled./agg_RZC_denominator_equi_gross_scaled;
        R_net_unscaled   = agg_R_numerator_equi_net_unscaled./agg_RZC_denominator_equi_net_unscaled;
        R_net_scaled     = agg_R_numerator_equi_net_scaled./agg_RZC_denominator_equi_net_scaled;

        % solve for Z
        if     isequal(EMUnits,'conc'); SO4conc2equimult = conc2equi.SO4;
        elseif isequal(EMUnits,'equi'); SO4conc2equimult = 1;
        end
        Z_gross_unscaled = SO4conc2equimult*SO4fromFeS2_gross_unscaled./agg_RZC_denominator_equi_gross_unscaled;
        Z_gross_scaled   = SO4conc2equimult*SO4fromFeS2_gross_scaled./agg_RZC_denominator_equi_gross_scaled;
        Z_net_unscaled   = SO4conc2equimult*SO4fromFeS2_net_unscaled./agg_RZC_denominator_equi_net_unscaled;
        Z_net_scaled     = SO4conc2equimult*SO4fromFeS2_net_scaled./agg_RZC_denominator_equi_net_scaled;

        % solve for Y
        Y_gross_unscaled = SO4fromFeS2_gross_unscaled./InvRivDatSO4;
        Y_gross_scaled   = SO4fromFeS2_gross_scaled./InvRivDatSO4;
        Y_net_unscaled   = SO4fromFeS2_net_unscaled./InvRivDatSO4;
        Y_net_scaled     = SO4fromFeS2_net_scaled./InvRivDatSO4;

        % solve for W    
        if     isequal(EMUnits,'equi')
               SO4conc2equimult  = 1; 
               InvRivDatnorm     = river.InvertedRiverValues.norm.all; 
        elseif isequal(EMUnits,'conc')
               SO4conc2equimult  = conc2equi.SO4; 
               InvRivDatnorm     = 0; 
               for lll=1:length(ObsInNormalization)
                    conc2equimult = eval(sprintf('conc2equi.%s',ObsInNormalization{lll}));
                    InvRivDatnorm = InvRivDatnorm + conc2equimult.*eval(sprintf('river.InvertedRiverValues.%s.all',ObsInNormalization{lll}));
               end
        end
        W_gross_unscaled = SO4conc2equimult.*SO4fromFeS2_gross_unscaled./InvRivDatnorm;
        W_gross_scaled   = SO4conc2equimult.*SO4fromFeS2_gross_scaled./InvRivDatnorm;
        W_net_unscaled   = SO4conc2equimult.*SO4fromFeS2_net_unscaled./InvRivDatnorm;
        W_net_scaled     = SO4conc2equimult.*SO4fromFeS2_net_scaled./InvRivDatnorm;  

        % solve for C
        if   isequal(carbonisotopematch,'DIC')
             C_gross_unscaled = DICfromCorg_gross_unscaled./agg_RZC_denominator_equi_gross_unscaled;
             C_gross_scaled   = DICfromCorg_gross_scaled./agg_RZC_denominator_equi_gross_scaled;
             C_net_unscaled   = DICfromCorg_net_unscaled./agg_RZC_denominator_equi_net_unscaled;
             C_net_scaled     = DICfromCorg_net_scaled./agg_RZC_denominator_equi_net_scaled;
        else
             C_gross_unscaled = HCO3fromCorg_gross_unscaled./agg_RZC_denominator_equi_gross_unscaled;
             C_gross_scaled   = HCO3fromCorg_gross_scaled./agg_RZC_denominator_equi_gross_scaled;
             C_net_unscaled   = HCO3fromCorg_net_unscaled./agg_RZC_denominator_equi_net_unscaled;
             C_net_scaled     = HCO3fromCorg_net_scaled./agg_RZC_denominator_equi_net_scaled;
        end

        % calculate parameters of interest
        R_gross_median_unscaled(1:s,1) = nanmedian(R_gross_unscaled',1)'; R_net_median_unscaled(1:s,1) = nanmedian(R_net_unscaled',1)'; R_gross_median_scaled(1:s,1) = nanmedian(R_gross_scaled',1)'; R_net_median_scaled(1:s,1)  = nanmedian(R_net_scaled',1)';
        R_gross_pct05_unscaled(1:s,1)  = prctile(R_gross_unscaled',05)';  R_net_pct05_unscaled(1:s,1)  = prctile(R_net_unscaled',05)';  R_gross_pct05_scaled(1:s,1)  = prctile(R_gross_scaled',05)';  R_net_pct05_scaled(1:s,1)   = prctile(R_net_scaled',05)'; 
        R_gross_pct25_unscaled(1:s,1)  = prctile(R_gross_unscaled',25)';  R_net_pct25_unscaled(1:s,1)  = prctile(R_net_unscaled',25)';  R_gross_pct25_scaled(1:s,1)  = prctile(R_gross_scaled',25)';  R_net_pct25_scaled(1:s,1)   = prctile(R_net_scaled',25)';
        R_gross_pct75_unscaled(1:s,1)  = prctile(R_gross_unscaled',75)';  R_net_pct75_unscaled(1:s,1)  = prctile(R_net_unscaled',75)';  R_gross_pct75_scaled(1:s,1)  = prctile(R_gross_scaled',75)';  R_net_pct75_scaled(1:s,1)   = prctile(R_net_scaled',75)';
        R_gross_pct95_unscaled(1:s,1)  = prctile(R_gross_unscaled',95)';  R_net_pct95_unscaled(1:s,1)  = prctile(R_net_unscaled',95)';  R_gross_pct95_scaled(1:s,1)  = prctile(R_gross_scaled',95)';  R_net_pct95_scaled(1:s,1)   = prctile(R_net_scaled',95)';

        Z_gross_median_unscaled(1:s,1) = nanmedian(Z_gross_unscaled',1)'; Z_net_median_unscaled(1:s,1) = nanmedian(Z_net_unscaled',1)'; Z_gross_median_scaled(1:s,1) = nanmedian(Z_gross_scaled',1)'; Z_net_median_scaled(1:s,1) = nanmedian(Z_net_scaled',1)';
        Z_gross_pct05_unscaled(1:s,1)  = prctile(Z_gross_unscaled',05)';  Z_net_pct05_unscaled(1:s,1)  = prctile(Z_net_unscaled',05)';  Z_gross_pct05_scaled(1:s,1)  = prctile(Z_gross_scaled',05)';  Z_net_pct05_scaled(1:s,1)  = prctile(Z_net_scaled',05)'; 
        Z_gross_pct25_unscaled(1:s,1)  = prctile(Z_gross_unscaled',25)';  Z_net_pct25_unscaled(1:s,1)  = prctile(Z_net_unscaled',25)';  Z_gross_pct25_scaled(1:s,1)  = prctile(Z_gross_scaled',25)';  Z_net_pct25_scaled(1:s,1)  = prctile(Z_net_scaled',25)'; 
        Z_gross_pct75_unscaled(1:s,1)  = prctile(Z_gross_unscaled',75)';  Z_net_pct75_unscaled(1:s,1)  = prctile(Z_net_unscaled',75)';  Z_gross_pct75_scaled(1:s,1)  = prctile(Z_gross_scaled',75)';  Z_net_pct75_scaled(1:s,1)  = prctile(Z_net_scaled',75)';
        Z_gross_pct95_unscaled(1:s,1)  = prctile(Z_gross_unscaled',95)';  Z_net_pct95_unscaled(1:s,1)  = prctile(Z_net_unscaled',95)';  Z_gross_pct95_scaled(1:s,1)  = prctile(Z_gross_scaled',95)';  Z_net_pct95_scaled(1:s,1)  = prctile(Z_net_scaled',95)';

        C_gross_median_unscaled(1:s,1) = nanmedian(C_gross_unscaled',1)'; C_net_median_unscaled(1:s,1) = nanmedian(C_net_unscaled',1)'; C_gross_median_scaled(1:s,1) = nanmedian(C_gross_scaled',1)'; C_net_median_scaled(1:s,1) = nanmedian(C_net_scaled',1)';
        C_gross_pct05_unscaled(1:s,1)  = prctile(C_gross_unscaled',05)';  C_net_pct05_unscaled(1:s,1)  = prctile(C_net_unscaled',05)';  C_gross_pct05_scaled(1:s,1)  = prctile(C_gross_scaled',05)';  C_net_pct05_scaled(1:s,1)  = prctile(C_net_scaled',05)'; 
        C_gross_pct25_unscaled(1:s,1)  = prctile(C_gross_unscaled',25)';  C_net_pct25_unscaled(1:s,1)  = prctile(C_net_unscaled',25)';  C_gross_pct25_scaled(1:s,1)  = prctile(C_gross_scaled',25)';  C_net_pct25_scaled(1:s,1)  = prctile(C_net_scaled',25)'; 
        C_gross_pct75_unscaled(1:s,1)  = prctile(C_gross_unscaled',75)';  C_net_pct75_unscaled(1:s,1)  = prctile(C_net_unscaled',75)';  C_gross_pct75_scaled(1:s,1)  = prctile(C_gross_scaled',75)';  C_net_pct75_scaled(1:s,1)  = prctile(C_net_scaled',75)';
        C_gross_pct95_unscaled(1:s,1)  = prctile(C_gross_unscaled',95)';  C_net_pct95_unscaled(1:s,1)  = prctile(C_net_unscaled',95)';  C_gross_pct95_scaled(1:s,1)  = prctile(C_gross_scaled',95)';  C_net_pct95_scaled(1:s,1)  = prctile(C_net_scaled',95)';

        W_gross_median_unscaled(1:s,1) = nanmedian(W_gross_unscaled',1)'; W_net_median_unscaled(1:s,1) = nanmedian(W_net_unscaled',1)'; W_gross_median_scaled(1:s,1) = nanmedian(W_gross_scaled',1)'; W_net_median_scaled(1:s,1) = nanmedian(W_net_scaled',1)';
        W_gross_pct05_unscaled(1:s,1)  = prctile(W_gross_unscaled',05)';  W_net_pct05_unscaled(1:s,1)  = prctile(W_net_unscaled',05)';  W_gross_pct05_scaled(1:s,1)  = prctile(W_gross_scaled',05)';  W_net_pct05_scaled(1:s,1)  = prctile(W_net_scaled',05)';
        W_gross_pct25_unscaled(1:s,1)  = prctile(W_gross_unscaled',25)';  W_net_pct25_unscaled(1:s,1)  = prctile(W_net_unscaled',25)';  W_gross_pct25_scaled(1:s,1)  = prctile(W_gross_scaled',25)';  W_net_pct25_scaled(1:s,1)  = prctile(W_net_scaled',25)';
        W_gross_pct75_unscaled(1:s,1)  = prctile(W_gross_unscaled',75)';  W_net_pct75_unscaled(1:s,1)  = prctile(W_net_unscaled',75)';  W_gross_pct75_scaled(1:s,1)  = prctile(W_gross_scaled',75)';  W_net_pct75_scaled(1:s,1)  = prctile(W_net_scaled',75)';
        W_gross_pct95_unscaled(1:s,1)  = prctile(W_gross_unscaled',95)';  W_net_pct95_unscaled(1:s,1)  = prctile(W_net_unscaled',95)';  W_gross_pct95_scaled(1:s,1)  = prctile(W_gross_scaled',95)';  W_net_pct95_scaled(1:s,1)  = prctile(W_net_scaled',95)';

        Y_gross_median_unscaled(1:s,1) = nanmedian(Y_gross_unscaled',1)'; Y_net_median_unscaled(1:s,1) = nanmedian(Y_net_unscaled',1)'; Y_gross_median_scaled(1:s,1) = nanmedian(Y_gross_scaled',1)'; Y_net_median_scaled(1:s,1) = nanmedian(Y_net_scaled',1)';
        Y_gross_pct05_unscaled(1:s,1)  = prctile(Y_gross_unscaled',05)';  Y_net_pct05_unscaled(1:s,1)  = prctile(Y_net_unscaled',05)';  Y_gross_pct05_scaled(1:s,1)  = prctile(Y_gross_scaled',05)';  Y_net_pct05_scaled(1:s,1)  = prctile(Y_net_scaled',05)'; 
        Y_gross_pct25_unscaled(1:s,1)  = prctile(Y_gross_unscaled',25)';  Y_net_pct25_unscaled(1:s,1)  = prctile(Y_net_unscaled',25)';  Y_gross_pct25_scaled(1:s,1)  = prctile(Y_gross_scaled',25)';  Y_net_pct25_scaled(1:s,1)  = prctile(Y_net_scaled',25)'; 
        Y_gross_pct75_unscaled(1:s,1)  = prctile(Y_gross_unscaled',75)';  Y_net_pct75_unscaled(1:s,1)  = prctile(Y_net_unscaled',75)';  Y_gross_pct75_scaled(1:s,1)  = prctile(Y_gross_scaled',75)';  Y_net_pct75_scaled(1:s,1)  = prctile(Y_net_scaled',75)'; 
        Y_gross_pct95_unscaled(1:s,1)  = prctile(Y_gross_unscaled',95)';  Y_net_pct95_unscaled(1:s,1)  = prctile(Y_net_unscaled',95)';  Y_gross_pct95_scaled(1:s,1)  = prctile(Y_gross_scaled',95)';  Y_net_pct95_scaled(1:s,1)  = prctile(Y_net_scaled',95)'; 

        % save all of those calculations and remove clutter from the workspace
        river.RZCWY.R.gross.unscaled.all    = R_gross_unscaled;        clear R_gross_unscaled;        river.RZCWY.R.net.unscaled.all    = R_net_unscaled;        clear R_net_unscaled;        river.RZCWY.R.gross.scaled.all    = R_gross_scaled;        clear R_gross_scaled;        river.RZCWY.R.net.scaled.all    = R_net_scaled;        clear R_net_scaled;  
        river.RZCWY.R.gross.unscaled.median = R_gross_median_unscaled; clear R_gross_median_unscaled; river.RZCWY.R.net.unscaled.median = R_net_median_unscaled; clear R_net_median_unscaled; river.RZCWY.R.gross.scaled.median = R_gross_median_scaled; clear R_gross_median_scaled; river.RZCWY.R.net.scaled.median = R_net_median_scaled; clear R_net_median_scaled;
        river.RZCWY.R.gross.unscaled.pct05  = R_gross_pct05_unscaled;  clear R_gross_pct05_unscaled;  river.RZCWY.R.net.unscaled.pct05  = R_net_pct05_unscaled;  clear R_net_pct05_unscaled;  river.RZCWY.R.gross.scaled.pct05  = R_gross_pct05_scaled;  clear R_gross_pct05_scaled;  river.RZCWY.R.net.scaled.pct05  = R_net_pct05_scaled;  clear R_net_pct05_scaled;
        river.RZCWY.R.gross.unscaled.pct25  = R_gross_pct25_unscaled;  clear R_gross_pct25_unscaled;  river.RZCWY.R.net.unscaled.pct25  = R_net_pct25_unscaled;  clear R_net_pct25_unscaled;  river.RZCWY.R.gross.scaled.pct25  = R_gross_pct25_scaled;  clear R_gross_pct25_scaled;  river.RZCWY.R.net.scaled.pct25  = R_net_pct25_scaled;  clear R_net_pct25_scaled;
        river.RZCWY.R.gross.unscaled.pct75  = R_gross_pct75_unscaled;  clear R_gross_pct75_unscaled;  river.RZCWY.R.net.unscaled.pct75  = R_net_pct75_unscaled;  clear R_net_pct75_unscaled;  river.RZCWY.R.gross.scaled.pct75  = R_gross_pct75_scaled;  clear R_gross_pct75_scaled;  river.RZCWY.R.net.scaled.pct75  = R_net_pct75_scaled;  clear R_net_pct75_scaled;
        river.RZCWY.R.gross.unscaled.pct95  = R_gross_pct95_unscaled;  clear R_gross_pct95_unscaled;  river.RZCWY.R.net.unscaled.pct95  = R_net_pct95_unscaled;  clear R_net_pct95_unscaled;  river.RZCWY.R.gross.scaled.pct95  = R_gross_pct95_scaled;  clear R_gross_pct95_scaled;  river.RZCWY.R.net.scaled.pct95  = R_net_pct95_scaled;  clear R_net_pct95_scaled;

        river.RZCWY.Z.gross.unscaled.all    = Z_gross_unscaled;        clear Z_gross_unscaled;        river.RZCWY.Z.net.unscaled.all    = Z_net_unscaled;        clear Z_net_unscaled;        river.RZCWY.Z.gross.scaled.all    = Z_gross_scaled;        clear Z_gross_scaled;        river.RZCWY.Z.net.scaled.all    = Z_net_scaled;        clear Z_net_scaled;
        river.RZCWY.Z.gross.unscaled.median = Z_gross_median_unscaled; clear Z_gross_median_unscaled; river.RZCWY.Z.net.unscaled.median = Z_net_median_unscaled; clear Z_net_median_unscaled; river.RZCWY.Z.gross.scaled.median = Z_gross_median_scaled; clear Z_gross_median_scaled; river.RZCWY.Z.net.scaled.median = Z_net_median_scaled; clear Z_net_median_scaled;
        river.RZCWY.Z.gross.unscaled.pct05  = Z_gross_pct05_unscaled;  clear Z_gross_pct05_unscaled;  river.RZCWY.Z.net.unscaled.pct05  = Z_net_pct05_unscaled;  clear Z_net_pct05_unscaled;  river.RZCWY.Z.gross.scaled.pct05  = Z_gross_pct05_scaled;  clear Z_gross_pct05_scaled;  river.RZCWY.Z.net.scaled.pct05  = Z_net_pct05_scaled;  clear Z_net_pct05_scaled;
        river.RZCWY.Z.gross.unscaled.pct25  = Z_gross_pct25_unscaled;  clear Z_gross_pct25_unscaled;  river.RZCWY.Z.net.unscaled.pct25  = Z_net_pct25_unscaled;  clear Z_net_pct25_unscaled;  river.RZCWY.Z.gross.scaled.pct25  = Z_gross_pct25_scaled;  clear Z_gross_pct25_scaled;  river.RZCWY.Z.net.scaled.pct25  = Z_net_pct25_scaled;  clear Z_net_pct25_scaled;
        river.RZCWY.Z.gross.unscaled.pct75  = Z_gross_pct75_unscaled;  clear Z_gross_pct75_unscaled;  river.RZCWY.Z.net.unscaled.pct75  = Z_net_pct75_unscaled;  clear Z_net_pct75_unscaled;  river.RZCWY.Z.gross.scaled.pct75  = Z_gross_pct75_scaled;  clear Z_gross_pct75_scaled;  river.RZCWY.Z.net.scaled.pct75  = Z_net_pct75_scaled;  clear Z_net_pct75_scaled;
        river.RZCWY.Z.gross.unscaled.pct95  = Z_gross_pct95_unscaled;  clear Z_gross_pct95_unscaled;  river.RZCWY.Z.net.unscaled.pct95  = Z_net_pct95_unscaled;  clear Z_net_pct95_unscaled;  river.RZCWY.Z.gross.scaled.pct95  = Z_gross_pct95_scaled;  clear Z_gross_pct95_scaled;  river.RZCWY.Z.net.scaled.pct95  = Z_net_pct95_scaled;  clear Z_net_pct95_scaled;

        river.RZCWY.C.gross.unscaled.all    = C_gross_unscaled;        clear C_gross_unscaled;        river.RZCWY.C.net.unscaled.all    = C_net_unscaled;        clear C_net_unscaled;        river.RZCWY.C.gross.scaled.all    = C_gross_scaled;        clear C_gross_scaled;        river.RZCWY.C.net.scaled.all    = C_net_scaled;        clear C_net_scaled;
        river.RZCWY.C.gross.unscaled.median = C_gross_median_unscaled; clear C_gross_median_unscaled; river.RZCWY.C.net.unscaled.median = C_net_median_unscaled; clear C_net_median_unscaled; river.RZCWY.C.gross.scaled.median = C_gross_median_scaled; clear C_gross_median_scaled; river.RZCWY.C.net.scaled.median = C_net_median_scaled; clear C_net_median_scaled;
        river.RZCWY.C.gross.unscaled.pct05  = C_gross_pct05_unscaled;  clear C_gross_pct05_unscaled;  river.RZCWY.C.net.unscaled.pct05  = C_net_pct05_unscaled;  clear C_net_pct05_unscaled;  river.RZCWY.C.gross.scaled.pct05  = C_gross_pct05_scaled;  clear C_gross_pct05_scaled;  river.RZCWY.C.net.scaled.pct05  = C_net_pct05_scaled;  clear C_net_pct05_scaled;
        river.RZCWY.C.gross.unscaled.pct25  = C_gross_pct25_unscaled;  clear C_gross_pct25_unscaled;  river.RZCWY.C.net.unscaled.pct25  = C_net_pct25_unscaled;  clear C_net_pct25_unscaled;  river.RZCWY.C.gross.scaled.pct25  = C_gross_pct25_scaled;  clear C_gross_pct25_scaled;  river.RZCWY.C.net.scaled.pct25  = C_net_pct25_scaled;  clear C_net_pct25_scaled;
        river.RZCWY.C.gross.unscaled.pct75  = C_gross_pct75_unscaled;  clear C_gross_pct75_unscaled;  river.RZCWY.C.net.unscaled.pct75  = C_net_pct75_unscaled;  clear C_net_pct75_unscaled;  river.RZCWY.C.gross.scaled.pct75  = C_gross_pct75_scaled;  clear C_gross_pct75_scaled;  river.RZCWY.C.net.scaled.pct75  = C_net_pct75_scaled;  clear C_net_pct75_scaled;
        river.RZCWY.C.gross.unscaled.pct95  = C_gross_pct95_unscaled;  clear C_gross_pct95_unscaled;  river.RZCWY.C.net.unscaled.pct95  = C_net_pct95_unscaled;  clear C_net_pct95_unscaled;  river.RZCWY.C.gross.scaled.pct95  = C_gross_pct95_scaled;  clear C_gross_pct95_scaled;  river.RZCWY.C.net.scaled.pct95  = C_net_pct95_scaled;  clear C_net_pct95_scaled;

        river.RZCWY.W.gross.unscaled.all    = W_gross_unscaled;        clear W_gross_unscaled;        river.RZCWY.W.net.unscaled.all    = W_net_unscaled;        clear W_net_unscaled;        river.RZCWY.W.gross.scaled.all    = W_gross_scaled;        clear W_gross_scaled;        river.RZCWY.W.net.scaled.all    = W_net_scaled;        clear W_net_scaled;
        river.RZCWY.W.gross.unscaled.median = W_gross_median_unscaled; clear W_gross_median_unscaled; river.RZCWY.W.net.unscaled.median = W_net_median_unscaled; clear W_net_median_unscaled; river.RZCWY.W.gross.scaled.median = W_gross_median_scaled; clear W_gross_median_scaled; river.RZCWY.W.net.scaled.median = W_net_median_scaled; clear W_net_median_scaled;
        river.RZCWY.W.gross.unscaled.pct05  = W_gross_pct05_unscaled;  clear W_gross_pct05_unscaled;  river.RZCWY.W.net.unscaled.pct05  = W_net_pct05_unscaled;  clear W_net_pct05_unscaled;  river.RZCWY.W.gross.scaled.pct05  = W_gross_pct05_scaled;  clear W_gross_pct05_scaled;  river.RZCWY.W.net.scaled.pct05  = W_net_pct05_scaled;  clear W_net_pct05_scaled;
        river.RZCWY.W.gross.unscaled.pct25  = W_gross_pct25_unscaled;  clear W_gross_pct25_unscaled;  river.RZCWY.W.net.unscaled.pct25  = W_net_pct25_unscaled;  clear W_net_pct25_unscaled;  river.RZCWY.W.gross.scaled.pct25  = W_gross_pct25_scaled;  clear W_gross_pct25_scaled;  river.RZCWY.W.net.scaled.pct25  = W_net_pct25_scaled;  clear W_net_pct25_scaled;
        river.RZCWY.W.gross.unscaled.pct75  = W_gross_pct75_unscaled;  clear W_gross_pct75_unscaled;  river.RZCWY.W.net.unscaled.pct75  = W_net_pct75_unscaled;  clear W_net_pct75_unscaled;  river.RZCWY.W.gross.scaled.pct75  = W_gross_pct75_scaled;  clear W_gross_pct75_scaled;  river.RZCWY.W.net.scaled.pct75  = W_net_pct75_scaled;  clear W_net_pct75_scaled;
        river.RZCWY.W.gross.unscaled.pct95  = W_gross_pct95_unscaled;  clear W_gross_pct95_unscaled;  river.RZCWY.W.net.unscaled.pct95  = W_net_pct95_unscaled;  clear W_net_pct95_unscaled;  river.RZCWY.W.gross.scaled.pct95  = W_gross_pct95_scaled;  clear W_gross_pct95_scaled;  river.RZCWY.W.net.scaled.pct95  = W_net_pct95_scaled;  clear W_net_pct95_scaled;

        river.RZCWY.Y.gross.unscaled.all    = Y_gross_unscaled;        clear Y_gross_unscaled;        river.RZCWY.Y.net.unscaled.all    = Y_net_unscaled;        clear Y_net_unscaled;        river.RZCWY.Y.gross.scaled.all    = Y_gross_scaled;        clear Y_gross_scaled;        river.RZCWY.Y.net.scaled.all    = Y_net_scaled;        clear Y_net_scaled;
        river.RZCWY.Y.gross.unscaled.median = Y_gross_median_unscaled; clear Y_gross_median_unscaled; river.RZCWY.Y.net.unscaled.median = Y_net_median_unscaled; clear Y_net_median_unscaled; river.RZCWY.Y.gross.scaled.median = Y_gross_median_scaled; clear Y_gross_median_scaled; river.RZCWY.Y.net.scaled.median = Y_net_median_scaled; clear Y_net_median_scaled;
        river.RZCWY.Y.gross.unscaled.pct05  = Y_gross_pct05_unscaled;  clear Y_gross_pct05_unscaled;  river.RZCWY.Y.net.unscaled.pct05  = Y_net_pct05_unscaled;  clear Y_net_pct05_unscaled;  river.RZCWY.Y.gross.scaled.pct05  = Y_gross_pct05_scaled;  clear Y_gross_pct05_scaled;  river.RZCWY.Y.net.scaled.pct05  = Y_net_pct05_scaled;  clear Y_net_pct05_scaled;
        river.RZCWY.Y.gross.unscaled.pct25  = Y_gross_pct25_unscaled;  clear Y_gross_pct25_unscaled;  river.RZCWY.Y.net.unscaled.pct25  = Y_net_pct25_unscaled;  clear Y_net_pct25_unscaled;  river.RZCWY.Y.gross.scaled.pct25  = Y_gross_pct25_scaled;  clear Y_gross_pct25_scaled;  river.RZCWY.Y.net.scaled.pct25  = Y_net_pct25_scaled;  clear Y_net_pct25_scaled;
        river.RZCWY.Y.gross.unscaled.pct75  = Y_gross_pct75_unscaled;  clear Y_gross_pct75_unscaled;  river.RZCWY.Y.net.unscaled.pct75  = Y_net_pct75_unscaled;  clear Y_net_pct75_unscaled;  river.RZCWY.Y.gross.scaled.pct75  = Y_gross_pct75_scaled;  clear Y_gross_pct75_scaled;  river.RZCWY.Y.net.scaled.pct75  = Y_net_pct75_scaled;  clear Y_net_pct75_scaled;
        river.RZCWY.Y.gross.unscaled.pct95  = Y_gross_pct95_unscaled;  clear Y_gross_pct95_unscaled;  river.RZCWY.Y.net.unscaled.pct95  = Y_net_pct95_unscaled;  clear Y_net_pct95_unscaled;  river.RZCWY.Y.gross.scaled.pct95  = Y_gross_pct95_scaled;  clear Y_gross_pct95_scaled;  river.RZCWY.Y.net.scaled.pct95  = Y_net_pct95_scaled;  clear Y_net_pct95_scaled;

        % tell the user that R, Z, C, W, and Y are calculated and saved in river
        disp('calculated R, Z, C, W, Y saved in the field "RZCWY"');
    end

else
    % tell the user that R, Z, C, W, and Y were not calculated
    disp('did not calculate R, Z, C, W, Y because CalculateRZCWY is not 1');

end % end of function