function [OrEMdist, TrEMdist, distEMdist0, distEMdist0min, distEMdist0max, smpadddist_mean, smpadddist_std, smpadddist_min, smpadddist_max, lgu_sign] = MEANDIR_makeEMdistributions(nOL, ObsList, nEM, EMList0, EM, EMdatasource, EndMembersWithNegativeRatios, CoupleFeS2SO4intoEM, NormalizationType, ObsInNormalization, IterateOver, indxisopos, ConvertDelta2RList, ConvertDelta2R)

        % this function calculates the end-member chemical distributions that are then used in the Monte Carlo inversion. additionally, 
        % this script generates variables that encode information on which variables are defined relative to samples and which are defined 
        % as fractionation relative to all other weathering inputs

        % define cell arrays to hold information on the distributions of end-member chemistry. 
        OrEMdist         = {};           % OrEMdist holds the original distributions, based on user-given properties
        TrEMdist         = {};           % TrEMdist holds the truncated distributions, after removing unphysical values of the original distributions
        distEMdist0      = {};           % distEMdist0 contains information on the type of distribution for each element in each end-member
        distEMdist0min   = {};           % when there are uniform distributions with at least one end-member relative to the sample, this matrix records which end-member minimum values are floating and which are fixed (minimum values)
        distEMdist0max   = {};           % when there are uniform distributions with at least one end-member relative to the sample, this matrix records which end-member maximum values are floating and which are fixed (maximum values)
        smpadddist_mean  = NaN(nOL,nEM); % in cases of a normal distribution relative to a sample, this matrix records the offset
        smpadddist_std   = NaN(nOL,nEM); % in cases of a normal distribution relative to a sample, this matrix records the standard deviation
        smpadddist_min   = NaN(nOL,nEM); % in cases of a uniform distribution relative to a sample, this matrix records the minimum
        smpadddist_max   = NaN(nOL,nEM); % in cases of a uniform distribution relative to a sample, this matrix records the maximum
        lgu_sign         = NaN(nOL,nEM); % in case of a log-uniform distribution with negative values, this matrix records the sign of the distribution

        for ii=1:nOL % iterate over each each dissolved constituent
        for jj=1:nEM % iterate over each end-member       

              % clear the current probability distribution and the truncated distribution
              clear pd;
              clear td; 

              % find and record the distribution type for the active observation and active end-member
              % this can change within this script depending on the user-selected values
              dist_type          = eval(sprintf('EM.%s.disttype.%s',EMdatasource,ObsList{ii}));
              distEMdist0{ii,jj} = dist_type; % save the type of the distribution (this will often be overwritten in the code below)

              % (1) if the distribution is normal
              if   isequal(dist_type,'NOR')

                   % (1.1) pull the initial mean and standard deviation of the normal distribution
                   meanentry = eval(sprintf('EM.%s.%s.%s_Men%s%s',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));
                   stdentry  = eval(sprintf('EM.%s.%s.%s_Sig%s%s',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));
                   % if the entry is not 123456789i, this is not a normal distribution relative to sample values
                   % if the entry is not 314159265i, this is not a normal distribution about fractionation

                   % (1.2) this first set of code is going to be used in the majority of cases for normal distributions
                   if     meanentry ~= 123456789i & meanentry ~= 314159265i

                          % (1.2.1) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                          if  sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                              conv        = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                              record_mean = (meanentry/1000+1)*conv;
                              record_std  = (stdentry/1000)*conv;
                          else
                              record_mean = meanentry;
                              record_std  = stdentry;
                          end                

                          % (1.2.2) save the distribution
                          pd = makedist('Normal','mu',record_mean,'sigma',record_std);

                   % (1.3) address the case where values should be replaced by sample values
                   elseif meanentry==123456789i & isequal(IterateOver,'Samples')

                          % (1.3.1) reset the distribution to be 0 +/- 0
                          pd = makedist('Normal','mu',0,'sigma',0);   
                          distEMdist0{ii,jj} = 'SMP-NOR';

                          % (1.3.2) pull the values for the offset from the sample and for the standard deviation
                          pull_mean = eval(sprintf('EM.%s.%s.%s_Men%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));
                          pull_std  = stdentry;

                          % (1.3.3) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                          if   sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                               conv        = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                               record_mean = (pull_mean/1000)*conv;
                               record_std  = (pull_std/1000)*conv;
                          else
                               record_mean = pull_mean;
                               record_std  = pull_std;
                          end     

                          % (1.3.4) save the values which will later be used to construct the end-member matrix
                          smpadddist_mean(ii,jj) = record_mean;
                          smpadddist_std(ii,jj)  = record_std;

                   % (1.4) address the case where we are considering a distribution of fractionation factors.
                   elseif meanentry==314159265i & isequal(IterateOver,'Samples')

                          % (1.4.1) pull the values for the fractionation and standard deviation
                          pull_mean = eval(sprintf('EM.%s.%s.%s_Men%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));
                          pull_std  = stdentry;

                          % (1.4.2) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                          if   sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                               conv        = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                               record_mean = (pull_mean/1000)*conv;
                               record_std  = (pull_std/1000)*conv;
                          else
                               record_mean = pull_mean;
                               record_std  = pull_std;
                          end

                          % (1.4.3) save the distribution
                          pd      = makedist('Normal','mu',record_mean,'sigma',record_std);   
                          distEMdist0{ii,jj} = 'EPS-NOR';

                   % (1.5) check for any issues with the distribution choices
                   elseif (meanentry==123456789i | meanentry==314159265i) & isequal(IterateOver,'End-members')
                          disp('warning! MEANDIR cannot match end-members chemistry to samples or use fractionation factors when IterateOver is "End-members"');
                          disp(sprintf('either change IterateOver to "Samples" or define ranges for observation: %s in end-member: %s',ObsList{ii},EMList0{jj}))
                   else
                          disp(sprintf('warning! could not determine normal distribution values for "%s" in "%s"',ObsList{ii},EMList0{jj}));
                   end

               % (2) if the distribution is uniform    
               elseif isequal(dist_type,'UNI')           

                      % (2.1) pull the initial minimum and maximum of the uniform distribution
                      minentry = eval(sprintf('EM.%s.%s.%s_Min%s%s',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType)); % Get the minimum
                      maxentry = eval(sprintf('EM.%s.%s.%s_Max%s%s',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType)); % Get the maximum   
                      % if the entry is not +/- 123456789, it is not replaced by sample values
                      % if the entry is not +/- 314159265, this entry is not about fractionation

                      % (2.2) address the case where both values are end-points of the distribution (most cases)
                      if      minentry~=-123456789 & maxentry~=123456789 & minentry~=-314159265 & maxentry~=314159265

                              % (2.2.1) when the minimum and maximum entry are equal
                              if    minentry==maxentry % if min equals max, set the distribution to be normal with zero error

                                    % (2.2.1.1) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                                    if  sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                        conv       = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                        record_min = (minentry/1000+1)*conv;
                                        record_max = (maxentry/1000+1)*conv;
                                    else
                                        record_min = minentry;
                                        record_max = maxentry;
                                    end       

                                    % (2.2.1.2) save the uniform distribution
                                    pd = makedist('Normal','mu',(record_min+record_max)/2,'sigma',0);   
                                    distEMdist0{ii,jj} = 'NOR';
                                    stdentry = 0; % the use of this line is to allow this end-member variable to skip the truncation stage

                              % (2.2.2) when the minimum and maximum entry are not equal 
                              else         
                                    % (2.2.2.1) check if the minimum entry is less than the maximum entry
                                    if    minentry<maxentry

                                        % (2.2.2.1.1) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                                          if  sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                              conv       = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                              record_min = (minentry/1000+1)*conv;
                                              record_max = (maxentry/1000+1)*conv;
                                          else
                                              record_min = minentry;
                                              record_max = maxentry;
                                          end    

                                          % (2.2.2.1.2) save the uniform distribution
                                          pd = makedist('Uniform','Lower',record_min,'Upper',record_max);

                                    % (2.2.2.2) give warnings if the minimum entry is less than the maximum entry
                                    else
                                          disp(sprintf('\nproblem with end-members! the lower normalized ratio for %s in %s is higher than its higher normalized ratio',ObsList{ii},EMList0{jj}));
                                          clear pd; clear t; % this is just to make sure no values are accidentally carried through the calculation
                                    end
                              end  

                      % now address the cases of uniform distribution but where the end-member selections are relative to the samples
                      % (2.3) first, address the case where the minimum is equal to the sample but the maximum is not. set the value of 
                      % smpadddist_min to be the desired offset from the sample, and set the value of smpadddist_max to be the maximum. 
                      elseif  minentry==-123456789 & maxentry~=123456789 & isequal(IterateOver,'Samples') & minentry~=-314159265 & maxentry~=314159265

                              % (2.3.1) reset the distributions
                              distEMdist0min{ii,jj} = 'SMP';
                              distEMdist0max{ii,jj} = 'UNI';  

                              % (2.3.2) pull the values
                              pull_min = eval(sprintf('EM.%s.%s.%s_Min%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));
                              pull_max = maxentry; 

                              % (2.3.3) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                              if   sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                   conv       = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                   record_min = (pull_min/1000)*conv;
                                   record_max = (pull_max/1000+1)*conv;
                              else
                                   record_min = pull_min;
                                   record_max = pull_max;
                              end

                              % (2.3.4) save the distribution values
                              smpadddist_min(ii,jj) = record_min;
                              smpadddist_max(ii,jj) = record_max;

                      % (2.4) second, address the case where the minimum is not equal to the sample but the maximum is. set the value of 
                      % smpadddist_min to be the minimum value and set the value of smpadddist_max to be the desired offset from the sample. 
                      elseif  minentry~=-123456789 & maxentry==123456789 & isequal(IterateOver,'Samples') & minentry~=-314159265 & maxentry~=314159265

                              % (2.4.1) reset the distributions
                              distEMdist0min{ii,jj} = 'UNI';
                              distEMdist0max{ii,jj} = 'SMP';

                              % (2.4.2) pull the values
                              pull_min = minentry;
                              pull_max = eval(sprintf('EM.%s.%s.%s_Max%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));

                              % (2.4.3) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                              if   sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                   conv       = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                   record_min = (pull_min/1000+1)*conv;
                                   record_max = (pull_max/1000)*conv;
                              else
                                   record_min = pull_min;
                                   record_max = pull_max;
                              end     

                              % (2.4.4) save the distribution values
                              smpadddist_min(ii,jj) = record_min;
                              smpadddist_max(ii,jj) = record_max;

                      % (2.5) third, address the case where the minimum and maximum are floating relative to each sample. 
                      % set the values of smpadddist_min and smpadddist_max to be the desired offsets from the sample.
                      elseif  minentry==-123456789 & maxentry==123456789 & isequal(IterateOver,'Samples') & minentry~=-314159265 & maxentry~=314159265

                              % (2.5.1) reset the distributions
                              distEMdist0min{ii,jj} = 'SMP';
                              distEMdist0max{ii,jj} = 'SMP';

                              % (2.5.2) pull the values
                              pull_min = eval(sprintf('EM.%s.%s.%s_Min%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));
                              pull_max = eval(sprintf('EM.%s.%s.%s_Max%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));

                              % (2.5.3) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                              if   sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                   conv       = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                   record_min = (pull_min/1000)*conv;
                                   record_max = (pull_max/1000)*conv;
                              else
                                   record_min = pull_min;
                                   record_max = pull_max;
                              end     

                              % (2.5.4) save the distribution values
                              smpadddist_min(ii,jj) = record_min;
                              smpadddist_max(ii,jj) = record_max;

                      % (2.6) fourth, address the case where the minimum and maximum refer to isotope effects to be applied later
                      elseif  minentry==-314159265 & maxentry==314159265

                              % (2.6.1) pull the minimum and maximum fractionation factors
                              pull_min = eval(sprintf('EM.%s.%s.%s_Min%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));
                              pull_max = eval(sprintf('EM.%s.%s.%s_Max%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));

                              % (2.6.2) address the case where the fractionation factor minimum and maximum are equal
                              if    pull_min==pull_max % if min equals max, set the distribution to be normal with zero error

                                    % (2.6.2.1) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                                    if   sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                         conv       = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                         record_min = (pull_min/1000)*conv;
                                         record_max = (pull_max/1000)*conv;
                                    else
                                         record_min = pull_min;
                                         record_max = pull_max;
                                    end                                   

                                    % (2.6.2.2) save the distribution
                                    pd = makedist('Normal','mu',(record_min+record_max)/2,'sigma',0);   
                                    distEMdist0{ii,jj} = 'EPS-NOR';

                              % (2.6.3) address the case where the fractionation factor minimum and maximum are not equal
                              else     
                                    % (2.6.3.1) check that the minimum is less than the maximum
                                    if    pull_min<pull_max

                                          % (2.6.3.1.1) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                                          if   sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                               conv       = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                               record_min = (pull_min/1000)*conv;
                                               record_max = (pull_max/1000)*conv;
                                          else
                                               record_min = pull_min;
                                               record_max = pull_max;
                                          end

                                          % (2.6.3.1.2) save the distribution 
                                          pd = makedist('Uniform','Lower',record_min,'Upper',record_max);
                                          distEMdist0{ii,jj} = 'EPS-UNI';

                                    % (2.6.3.2) if the minimum is greater than the maximum, display an erorr
                                    else
                                          disp(sprintf('\nproblem with end-members! the lower normalized ratio for "%s" in "%s" is higher than its higher normalized ratio',ObsList{ii},EMList0{jj}));
                                          clear pd; clear t; % this is just to make sure no values are accidentally carried through the calculation
                                    end
                              end   
                      else;   disp(sprintf('warning! could not determine uniform distribution values for "%s" in "%s"',ObsList{ii},EMList0{jj}));
                      end      

                      % (2.7) if either uniform end-member is floating relative to the sample value, set 
                      % the entry of distEMdist0 to indicate that and set pd equal to 0+/-0
                      if     (minentry==-123456789 | maxentry==123456789)
                             pd = makedist('Normal','mu',0,'sigma',0);   
                             distEMdist0{ii,jj} = 'SMP-UNI';
                      end       
                      % (2.8) check for warnings to display
                      if (minentry==-123456789 | maxentry==123456789 | minentry==-314159265 | maxentry==314159265) & isequal(IterateOver,'End-members')
                         disp('warning! MEANDIR cannot match end-members chemistry to samples when IterateOver is "End-members"');
                         disp(sprintf('either change IterateOver to "Samples" or define ranges for observation: %s in end-member: %s',ObsList{ii},EMList0{jj}))
                      end                   

               % (3) if the distribution is log uniform   
               elseif isequal(dist_type,'LGU') 

                      % (3.1) pull the initial minimum and maximum of the log-uniform distribution
                      minentry = eval(sprintf('EM.%s.%s.%s_Min%s%s',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));
                      maxentry = eval(sprintf('EM.%s.%s.%s_Max%s%s',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));

                      % (3.2) address the case where both values are end-points of the distribution (most cases)
                      % acccess this first portion of the code if none of the end-member values should be relative to the sample chemistry
                      if    minentry~=-123456789 & maxentry~=123456789 & minentry~=-314159265 & maxentry~=314159265

                                % (3.2.1) when the minimum and maximum entry are equal
                                if   minentry==maxentry    

                                     % (3.2.1.1) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                                     if  sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                         conv       = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                         record_min = (minentry/1000+1)*conv;
                                         record_max = (maxentry/1000+1)*conv;
                                     else
                                         record_min = minentry;
                                         record_max = maxentry;
                                     end     

                                     % (3.2.1.2) save the uniform distribution value
                                     pd = makedist('Normal','mu',(record_min+record_max)/2,'sigma',0);
                                     distEMdist0{ii,jj} = 'NOR';
                                     stdentry = 0; % the use of this line is to allow this end-member variable to skip the truncation stage

                                % (3.2.2) when the minimum and maximum entry are not equal
                                else

                                     % (3.2.2.1) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                                     if  sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                         conv        = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                         record_min  = (minentry/1000+1)*conv;
                                         record_max  = (maxentry/1000+1)*conv;
                                     else
                                         record_min  = minentry;
                                         record_max  = maxentry;
                                     end     

                                     % (3.2.2.2) consider if either entry is 0, which will not work for taking a log
                                     if  record_min==0 | record_max==0
                                          disp('  '); 
                                          disp('error! cannot define a log-normal distribution when only one of the values is zero'); 
                                          disp(sprintf('check variable "%s" in end-member "%s" within the EMgroup "%s"',ObsList{ii},EMList0{jj},EMdatasource));
                                          disp('  '); 
                                          break;

                                     % (3.2.2.3) access this portion of the code if there is a log-normal distribution with different values and neither is equal to zero
                                     else    
                                          % (3.2.2.3.1) take the natural log of the end-member distribution values
                                          if  record_min<record_max
                                              if   record_min<0 | record_max < 0; isneg = 1; lgu_sign(ii,jj) = -1;
                                              else;                               isneg = 0; lgu_sign(ii,jj) = 1;
                                              end
                                              if isneg == 1
                                                  v1 = log(-record_min);
                                                  v2 = log(-record_max);
                                                  if v1<v2; lowval = v1; highval = v2; end
                                                  if v2<v1; lowval = v2; highval = v1; end
                                              else
                                                  v1 = log(record_min); 
                                                  v2 = log(record_max);
                                                  if v1<v2; lowval = v1; highval = v2; end
                                                  if v2<v1; lowval = v2; highval = v1; end
                                              end
                                              pd = makedist('Uniform','Lower',lowval,'Upper',highval);  

                                          % (3.2.2.3.2) generate an error if the minimum is greater than the maximum
                                          else
                                              disp(sprintf('\nproblem with end-members! the lower normalized ratio for %s in %s is higher than its higher normalized ratio',ObsList{ii},EMList0{jj}));
                                              clear pd; clear t; % this is just to make sure no values are accidentally carried through the calculation
                                          end
                                     end
                                end

                      % now address the cases of log-uniform distribution but where the end-member selections are relative to the samples
                      % (3.3) first, address the case where the minimum is equal to the sample but the maximum is not. set the value of 
                      % smpadddist_min to be the desired offset from the sample, and set the value of smpadddist_max to be the maximum.
                       elseif   minentry==-123456789 & maxentry~=123456789 & isequal(IterateOver,'Samples') & minentry~=-314159265 & maxentry~=314159265

                                % (3.3.1) reset the distributions
                                distEMdist0min{ii,jj} = 'SMP';
                                distEMdist0max{ii,jj} = 'LGU';

                                % (3.3.2) pull the values
                                pull_min = eval(sprintf('EM.%s.%s.%s_Min%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));
                                pull_max = maxentry; 

                                % (3.3.3) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                                if   sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                     conv       = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                     record_min = (pull_min/1000)*conv;
                                     record_max = (pull_max/1000+1)*conv;
                                else
                                     record_min = pull_min;
                                     record_max = pull_max;
                                end     

                                % (3.3.4) save the distribution values
                                smpadddist_min(ii,jj) = record_min;
                                smpadddist_max(ii,jj) = record_max;

                      % (3.4) second, address the case where the minimum is not equal to the sample but the maximum is. set the value of 
                      % smpadddist_min to be the minimum value and set the value of smpadddist_max to be the desired offset from the sample
                       elseif   minentry~=-123456789 & maxentry==123456789 & isequal(IterateOver,'Samples') & minentry~=-314159265 & maxentry~=314159265

                                % (3.4.1) reset the distributions
                                distEMdist0min{ii,jj} = 'LGU';
                                distEMdist0max{ii,jj} = 'SMP';  

                                % (3.4.2) pull the values
                                pull_min = minentry;
                                pull_max = eval(sprintf('EM.%s.%s.%s_Max%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType)); 

                                % (3.4.3) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                                if   sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                     conv       = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                     record_min = (pull_min/1000+1)*conv;
                                     record_max = (pull_max/1000)*conv;
                                else
                                     record_min = pull_min;
                                     record_max = pull_max;
                                end

                                % (3.4.4) save the distribution values
                                smpadddist_min(ii,jj) = record_min;
                                smpadddist_max(ii,jj) = record_max;

                      % (3.5) third, address the case where the minimum and maximum are floating relative to each sample 
                      % set the values of smpadddist_min and smpadddist_max to be the desired offsets from the sample
                       elseif   minentry==-123456789 & maxentry==123456789 & isequal(IterateOver,'Samples') & minentry~=-314159265 & maxentry~=314159265

                                % (3.5.1) reset the distributions
                                distEMdist0min{ii,jj} = 'SMP';
                                distEMdist0max{ii,jj} = 'SMP';

                                % (3.5.2) pull the values
                                pull_min = eval(sprintf('EM.%s.%s.%s_Min%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));
                                pull_max = eval(sprintf('EM.%s.%s.%s_Max%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));

                                % (3.5.3) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                                if   sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                     conv      = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                     record_min   = (pull_min/1000)*conv;
                                     record_max   = (pull_max/1000)*conv;
                                else
                                     record_min  = pull_min;
                                     record_max  = pull_max;
                                end    

                                % (3.5.4) save the distribution values
                                smpadddist_min(ii,jj) = record_min;
                                smpadddist_max(ii,jj) = record_max;

                      % (3.6) fourth, address the case where the minimum and maximum refer to isotope effects to be applied later
                      elseif   minentry==-314159265 & maxentry==314159265

                               % (3.6.1) pull the minimum and maximum fractionation factors
                               pull_min = eval(sprintf('EM.%s.%s.%s_Min%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));
                               pull_max = eval(sprintf('EM.%s.%s.%s_Max%s%s_addvalue',EMdatasource,EMList0{jj},dist_type,ObsList{ii},NormalizationType));

                               % (3.6.2) address the case where the fractionation factor minimum and maximum are equal
                                if   pull_min==pull_max

                                    % (3.6.2.1) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                                    if   sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                         conv       = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                         record_min = (pull_min/1000)*conv;
                                         record_max = (pull_max/1000)*conv;
                                    else
                                         record_min = pull_min;
                                         record_max = pull_max;
                                    end

                                    % (3.6.2.2) save the distribution
                                    pd = makedist('Normal','mu',(record_min+record_max)/2,'sigma',0);
                                    distEMdist0{ii,jj} = 'EPS-NOR';

                                % (3.6.3) address the case where the fractionation factor minimum and maximum are not equal     
                                else           

                                     % (3.6.3.1) evaluate whether the current observation is delta notation that should be converted to isotopic ratios
                                     if   sum(ismember(ConvertDelta2RList,ObsList{ii}))>0
                                          conv       = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,ObsList{ii}));
                                          record_min = (pull_min/1000)*conv;
                                          record_max = (pull_max/1000)*conv;
                                     else
                                          record_min = pull_min;
                                          record_max = pull_max;
                                     end

                                    % (3.6.3.2) consider if either entry is 0, which will not work for taking a log
                                     if record_min==0 | record_max==0
                                          disp('  '); 
                                          disp('error! cannot define a log-normal distribution when only one of the values is zero'); 
                                          disp(sprintf('check variable "%s" in end-member "%s" within the EMgroup "%s"',ObsList{ii},EMList0{jj},EMdatasource));
                                          disp('  '); 
                                          break;

                                     % (3.2.2.3) access this portion of the code if there is a log-normal distribution with different values and neither is equal to zero     
                                     else      
                                          if   record_min<record_max
                                               if   record_min<0 | record_max < 0; isneg = 1; lgu_sign(ii,jj) = -1;
                                               else;                               isneg = 0; lgu_sign(ii,jj) = 1;
                                               end         
                                               if isneg==1
                                                   v1 = log(-record_min); 
                                                   v2 = log(-record_max);
                                                   if v1<v2; lowval = v1; highval = v2; end
                                                   if v2<v1; lowval = v2; highval = v1; end
                                               else
                                                   v1 = log(record_min); 
                                                   v2 = log(record_max);
                                                   if v1<v2; lowval = v1; highval = v2; end
                                                   if v2<v1; lowval = v2; highval = v1; end
                                               end                                           
                                               pd = makedist('Uniform','Lower',lowval,'Upper',highval);  
                                               distEMdist0{ii,jj} = 'EPS-LGU';
                                          else                    
                                               disp(sprintf('\nproblem with end-members! the lower normalized ratio for %s in %s is higher than its higher normalized ratio',ObsList{ii},EMList0{jj}));
                                               clear pd; clear t; % this is just to make sure no values are accidentally carried through the calculation
                                          end
                                     end
                                end      
                      else;     disp(sprintf('warning! could not determine log-uniform distribution values for "%s" in "%s"',ObsList{ii},EMList0{jj}));
                      end           
                      % (3.7) if either log-uniform end-member is floating, set the entry
                      % of distEMdist0 to indicate that and set pd equal to zero
                      if  (minentry==-123456789 | maxentry==123456789)
                           distEMdist0{ii,jj} = 'SMP-LGU';   
                           pd = makedist('Normal','mu',0,'sigma',0);   
                      end 
                      % (3.8) check for warnings to display
                      if (minentry==-123456789 | maxentry==123456789 | minentry==-314159265 | maxentry==314159265) & isequal(IterateOver,'End-members')
                             disp('warning! MEANDIR cannot match end-members chemistry to samples when IterateOver is "End-members"');
                             disp(sprintf('either change IterateOver to "Samples" or define ranges for observation: %s in end-member: %s',ObsList{ii},EMList0{jj}))
                      end           
               end

               % evaluate whether normal distributions need to be truncated
               if isequal(dist_type,'NOR') & stdentry~=0                         
                  % first consider the lower bound
                  if      sum(indxisopos==ii)>0 % this indicates that the active observation is isotopic information
                          if isequal(ObsList{ii},'d7Li')  & sum(ismember(ConvertDelta2RList,ObsList{ii}))==0;  lowbound = -1000; end % delta values (in per mil) are defined between -1000 and positive infinity
                          if isequal(ObsList{ii},'d13C')  & sum(ismember(ConvertDelta2RList,ObsList{ii}))==0;  lowbound = -1000; end % delta values (in per mil) are defined between -1000 and positive infinity
                          if isequal(ObsList{ii},'d18O')  & sum(ismember(ConvertDelta2RList,ObsList{ii}))==0;  lowbound = -1000; end % delta values (in per mil) are defined between -1000 and positive infinity
                          if isequal(ObsList{ii},'d26Mg') & sum(ismember(ConvertDelta2RList,ObsList{ii}))==0;  lowbound = -1000; end % delta values (in per mil) are defined between -1000 and positive infinity
                          if isequal(ObsList{ii},'d30Si') & sum(ismember(ConvertDelta2RList,ObsList{ii}))==0;  lowbound = -1000; end % delta values (in per mil) are defined between -1000 and positive infinity
                          if isequal(ObsList{ii},'d34S')  & sum(ismember(ConvertDelta2RList,ObsList{ii}))==0;  lowbound = -1000; end % delta values (in per mil) are defined between -1000 and positive infinity
                          if isequal(ObsList{ii},'d42Ca') & sum(ismember(ConvertDelta2RList,ObsList{ii}))==0;  lowbound = -1000; end % delta values (in per mil) are defined between -1000 and positive infinity
                          if isequal(ObsList{ii},'d44Ca') & sum(ismember(ConvertDelta2RList,ObsList{ii}))==0;  lowbound = -1000; end % delta values (in per mil) are defined between -1000 and positive infinity
                          if isequal(ObsList{ii},'d56Fe') & sum(ismember(ConvertDelta2RList,ObsList{ii}))==0;  lowbound = -1000; end % delta values (in per mil) are defined between -1000 and positive infinity
                          if isequal(ObsList{ii},'d98Mo') & sum(ismember(ConvertDelta2RList,ObsList{ii}))==0;  lowbound = -1000; end % delta values (in per mil) are defined between -1000 and positive infinity

                          if isequal(ObsList{ii},'d7Li')  & sum(ismember(ConvertDelta2RList,ObsList{ii}))==1;  lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity
                          if isequal(ObsList{ii},'d13C')  & sum(ismember(ConvertDelta2RList,ObsList{ii}))==1;  lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity
                          if isequal(ObsList{ii},'d18O')  & sum(ismember(ConvertDelta2RList,ObsList{ii}))==1;  lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity
                          if isequal(ObsList{ii},'d26Mg') & sum(ismember(ConvertDelta2RList,ObsList{ii}))==1;  lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity
                          if isequal(ObsList{ii},'d30Si') & sum(ismember(ConvertDelta2RList,ObsList{ii}))==1;  lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity
                          if isequal(ObsList{ii},'d34S')  & sum(ismember(ConvertDelta2RList,ObsList{ii}))==1;  lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity
                          if isequal(ObsList{ii},'d42Ca') & sum(ismember(ConvertDelta2RList,ObsList{ii}))==1;  lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity
                          if isequal(ObsList{ii},'d44Ca') & sum(ismember(ConvertDelta2RList,ObsList{ii}))==1;  lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity
                          if isequal(ObsList{ii},'d56Fe') & sum(ismember(ConvertDelta2RList,ObsList{ii}))==1;  lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity
                          if isequal(ObsList{ii},'d98Mo') & sum(ismember(ConvertDelta2RList,ObsList{ii}))==1;  lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity
                          if isequal(ObsList{ii},'Fmod');                                                      lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity
                          if isequal(ObsList{ii},'Sr8786');                                                    lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity
                          if isequal(ObsList{ii},'Os8788');                                                    lowbound = 00000; end % isotopic ratios are defined between 0 and positive infinity

                  elseif  sum(ismember(EndMembersWithNegativeRatios,EMList0{jj}))>0
                          lowbound = -inf;   
                  else
                          lowbound = 0; % initial lower bounding ranges
                  end      
                  % now consider the upper bound
                  % if normalizing with SumObs and the ratio numerator contributes to the normalization, the upper bound is 1
                  if     isequal(NormalizationType,'SumObs') & sum(ismember(ObsInNormalization,ObsList{ii}))>0
                         highbound = 1;
                  else
                         highbound = inf;
                  end  

                  % having found the upper and lower bounds, truncate the original normal distribution
                   td = truncate(pd,lowbound,highbound);
               else
                   td = pd; 
               end

               OrEMdist{ii,jj} = pd; % store the original distribution
               TrEMdist{ii,jj} = td; % store the truncated/original distrubtion

               % accounting for the SO4 of FeS2 oxidation when forming end-member distrubtions
               % for now just set the SO4 values to equal the pyrite distribution
               % (the d34S of these end-members is set later, when actually pulling values and we 
               % don't have to worry about isotope conversion here because its just SO4, not d34S)
               if isequal(ObsList{ii},'SO4') == 1 & sum(ismember(CoupleFeS2SO4intoEM,EMList0{jj}))>0
                   disttype = eval(sprintf('EM.%s.disttype.FeS2SO4',EMdatasource));
                   if        isequal(disttype,'NOR')
                             pyr_SO4_mean = eval(sprintf('EM.%s.pyri.%s_MenFeS2SO4%s',EMdatasource,disttype,NormalizationType));
                             pyr_SO4_std  = eval(sprintf('EM.%s.pyri.%s_SigFeS2SO4%s',EMdatasource,disttype,NormalizationType));
                             pd = makedist('Normal','mu',pyr_SO4_mean,'sigma',pyr_SO4_std);
                   elseif   isequal(disttype,'UNI')
                             pyr_SO4_min = eval(sprintf('EM.%s.pyri.%s_MinFeS2SO4%s',EMdatasource,disttype,NormalizationType));
                             pyr_SO4_max = eval(sprintf('EM.%s.pyri.%s_MaxFeS2SO4%s',EMdatasource,disttype,NormalizationType));
                             if pyr_SO4_min==pyr_SO4_max; pd = makedist('Normal','mu',(pyr_SO4_min+pyr_SO4_max)/2,'sigma',0);
                             else;                        pd = makedist('Uniform','Lower',pyr_SO4_min,'Upper',pyr_SO4_max); 
                             end  
                   end
                   TrEMdist{ii,jj} = pd; 
               end            
        end % end iteration over end-members (indexed on j)
        end % end iteration over ions (indexed on i)

        [a, b] = size(TrEMdist); 
        if a == 1 & b ==1
            disp('problem! TrEMdist is only a single value, so something went wrong trying to construct the chemical distributions of end-members. check the spreadsheet of end-member information, the matricies EMMatrix and EMUncMatrix, and the function MEANDIR_makeEMdistributions'); 
        end

end % end of function