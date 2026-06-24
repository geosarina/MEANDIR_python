function  [MinFractionalContribution, MaxFractionalContribution] = MEANDIR_resetDICcont(ObsList, EMList0, RiverColumn, DegasDICContributionMin, DegasDICContributionMax, MinFractionalContribution, MaxFractionalContribution)          

          % this script resets values of MinFractionalContribution and MaxFractionalContribution to limit  
          % the amount of DIC that can be lost through degassing. the logic of the equations is as follows.
          %
          % if 100% of DIC were sourced through the Corg oxidation end-member and that end-member did not source any additional
          % elements to the normalization variable,  the maximum contribution of the Corg oxidation end-member to the normalization
          % variable would be DIC/normalization. for example, if the normalization were Ca + DIC in equal parts, 
          % the maximum contribution of Corg oxidation to the normalization variable would be 0.5. however, we normally constrain this value
          % as 100% within the range given to the optimizer just to keep things simple. now assume for the same case that I want to limit the
          % contribution to 35% of DIC. in this case, I would just have to multiply the desired fraction by the maximum fractional contribution, 
          % 0.35*DIC./(DIC + Ca) = 0.35*0.5 = 0.175. 
          %
          % stated slightly differently, the normalized river value gives the contribution to the normalization variable 
          % if the end-member sourced 100% of the DIC and nothing else. by multiplying that ratio by the desired maximum fraction of the DIC,
          % we calculate the maximum contribution to the normalization variable to generate the desired maximum contribution to the
          % DIC budget. 
                   
          river_DIC_norm = RiverColumn(ismember(ObsList,'DIC'));             % get the river DIC/normalization ratio
          FracDICMin     = DegasDICContributionMin*river_DIC_norm;           % multiply the river DIC/Norm by the minimum fractional contribution
          FracDICMax     = DegasDICContributionMax*river_DIC_norm;           % multiply the river DIC/Norm by the maximum fractional contribution
          MinFractionalContribution(ismember(EMList0,'degas')) = FracDICMin; % reset the fraction of the normalization that can be supplied by degassing
          MaxFractionalContribution(ismember(EMList0,'degas')) = FracDICMax; % reset the fraction of the normalization that can be supplied by degassing
end