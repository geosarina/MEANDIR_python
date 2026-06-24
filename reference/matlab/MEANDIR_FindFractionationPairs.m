function fractionation = MEANDIR_FindFractionationPairs(distEMdist0, distEMdist_r, ObsList, ObsList_r, EMList0, EMList_r, carbonisotopematch, indxisopos)

         % this script calculates the structure "fractionation", which contains information
         % on the relationship of isotopes and elements in the inversion
         % similar calculations are done twice in this script, once on the
         % reduced inversion variables and once on the full inversion variables

         clear fractionation; 
         fractionation = struct;

         % (1) first, consider the occurence of fractionation in the 'reduced'
         %  data, which is the set of results that will actually be inverted
         fractionation.red.fracdist_r = [ismember(distEMdist_r,'EPS-NOR') | ismember(distEMdist_r,'EPS-UNI') | ismember(distEMdist_r,'EPS-LGU')];
         fractionation.red.fracloc_r  = find(fractionation.red.fracdist_r);
         fractionation.red.n_r        = length(fractionation.red.fracloc_r);

         % calculate isotopeposition_r and indxisopos_r based on the presence of isotopes in ObsList_r
         isotopeposition_r = (ismember(ObsList_r,'d34S')   | ismember(ObsList_r,'d18O')   | ismember(ObsList_r,'Sr8786') | ...
                              ismember(ObsList_r,'d7Li')   | ismember(ObsList_r,'d26Mg')  | ismember(ObsList_r,'d30Si')  | ...
                              ismember(ObsList_r,'d42Ca')  | ismember(ObsList_r,'d44Ca')  | ismember(ObsList_r,'d56Fe')  | ...
                              ismember(ObsList_r,'d98Mo')  | ismember(ObsList_r,'Os8788') | ismember(ObsList_r,'d13C')  | ismember(ObsList_r,'Fmod'));
         indxisopos_r      = find(isotopeposition_r); % Find the index vector for isotope data           
         
         % iterate over each fractionation that is included in the reduced data
         for jj=1:fractionation.red.n_r

                % pull the active location within the matrix
                activeloc = fractionation.red.fracloc_r(jj);

                % the following loop finds the x, y position of a '1' in a matrix of zero and ones
                % after running the calculation, activeloc is the row position and count is the column position
                count = 0;
                while activeloc>0
                      activeloc = activeloc-size(fractionation.red.fracdist_r,1);
                      count = count+1; 
                end
                activeloc = activeloc+size(fractionation.red.fracdist_r,1);

                % for this instance of fractionation, calculate the isotope, the end-member, and the position of the isotope
                % and of the end-member in both the full inversion matrix and in the reduced inversion matrix
                fractionation.red.iso{jj,1}      = ObsList_r{activeloc};
                fractionation.red.EM{jj,1}       = EMList_r{count};                
                fractionation.red.isopos_r(jj,1) = find(ismember(ObsList_r,ObsList_r{activeloc}));
                fractionation.red.empos_r(jj,1)  = find(ismember(EMList_r,EMList_r{count}));
                fractionation.red.isopos0(jj,1)  = find(ismember(ObsList,ObsList_r{activeloc}));
                fractionation.red.empos0(jj,1)   = find(ismember(EMList0,EMList_r{count}));
                if isequal(ObsList_r{activeloc},'d7Li');                                       fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'Li'));   fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'Li'));   end % if d7Li, find Li
                if isequal(ObsList_r{activeloc},'d18O');                                       fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'SO4'));  fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'SO4'));  end % if d18O, find SO4
                if isequal(ObsList_r{activeloc},'d26Mg');                                      fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'Mg'));   fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'Mg'));   end % if d26Mg, find Mg
                if isequal(ObsList_r{activeloc},'d30Si');                                      fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'Si'));   fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'Si'));   end % if d30Si, find Si
                if isequal(ObsList_r{activeloc},'d34S');                                       fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'SO4'));  fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'SO4'));  end % if d34S, find SO4
                if isequal(ObsList_r{activeloc},'d42Ca');                                      fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'Ca'));   fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'Ca'));   end % if d42Ca, find Ca
                if isequal(ObsList_r{activeloc},'d44Ca');                                      fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'Ca'));   fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'Ca'));   end % if d44Ca, find Ca
                if isequal(ObsList_r{activeloc},'d56Fe');                                      fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'Fe'));   fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'Fe'));   end % if d56Fe, find Fe
                if isequal(ObsList_r{activeloc},'Sr8786');                                     fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'Sr'));   fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'Sr'));   end % if Sr87/86, find Sr
                if isequal(ObsList_r{activeloc},'Os8788');                                     fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'Os'));   fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'Os'));   end % if Os8788, find Os
                if isequal(ObsList_r{activeloc},'d98Mo');                                      fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'Mo'));   fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'Mo'));   end % if d98Mo, find Mo
                if isequal(ObsList_r{activeloc},'d13C') & isequal(carbonisotopematch,'HCO3');  fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'HCO3')); fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'HCO3')); end % if d13C and carbon variable is HCO3, find HCO3
                if isequal(ObsList_r{activeloc},'Fmod') & isequal(carbonisotopematch,'HCO3');  fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'HCO3')); fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'HCO3')); end % if Fmod and carbon variable is HCO3, find HCO3
                if isequal(ObsList_r{activeloc},'d13C') & isequal(carbonisotopematch,'DIC');   fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'DIC'));  fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'DIC'));  end % if d13C and carbon variable is DIC, find DIC
                if isequal(ObsList_r{activeloc},'Fmod') & isequal(carbonisotopematch,'DIC');   fractionation.red.ionpos0(jj,1) = find(ismember(ObsList,'DIC'));  fractionation.red.ionpos_r(jj,1) = find(ismember(ObsList_r,'DIC'));  end % if Fmod and carbon variable is DIC, find DIC

                 % check that fractionation has only been requested for isotopic data
                 if    sum(ismember(fractionation.red.isopos_r,indxisopos_r))~=length(fractionation.red.isopos_r)
                       problemobs = fractionation.red.iso(~ismember(fractionation.red.isopos_r,indxisopos_r));
                       problemEMs = fractionation.red.EM(~ismember(fractionation.red.isopos_r,indxisopos_r));
                       for j=1:length(problemobs)          
                           disp(sprintf('warning! MEANDIR is attempting to apply fractionation to "%s" in "%s" (not to isotopes, to concentration). this warning was generated in MEANDIR_FindFractionationPairs.',problemobs{j},problemEMs{j}));
                       end
                 end   % end of checking for isotopic data
         end           % end of iterating over instances of fractionation


         % (2) second, consider the occurence of fractionation in the 'all' data, which is the full inversion matrix 
         % all instances of fractionation in the reduced matrix will appear in the full inversion matrix, but not all
         % instances in the full matrix must appear in the reduced matrix
         fractionation.all.fracdist0    = [ismember(distEMdist0,'EPS-NOR')  | ismember(distEMdist0,'EPS-UNI')  | ismember(distEMdist0,'EPS-LGU')];
         fractionation.all.fracloc0     = find(fractionation.all.fracdist0);
         fractionation.all.n0           = length(fractionation.all.fracloc0);

         % iterate over each fractionation that is included in the full data
         for jj=1:length(fractionation.all.fracloc0)

                % pull the active location within the matrix
                activeloc = fractionation.all.fracloc0(jj);

                % the following loop finds the x, y position of a '1' in a matrix of zero and ones 
                % after running the calculation, activeloc is the row position and count is the column position
                count = 0;
                while activeloc>0
                      activeloc = activeloc-size(fractionation.all.fracdist0,1);
                      count = count+1; 
                end
                activeloc = activeloc+size(fractionation.all.fracdist0,1);

                % for this instance of fractionation, calculate the isotope, the end-member, and the position of the isotope
                % and of the end-member in both the full inversion matrix and in the reduced inversion matrix
                fractionation.all.iso{jj,1}      = ObsList{activeloc};
                fractionation.all.EM{jj,1}       = EMList0{count};
                fractionation.all.isopos0(jj,1)  = find(ismember(ObsList,ObsList{activeloc}));
                fractionation.all.empos0(jj,1)   = find(ismember(EMList0,EMList0{count}));
                testisopos_r = find(ismember(ObsList_r,ObsList{activeloc}));
                testtempos_r = find(ismember(EMList_r,EMList0{count}));
                if   ~isempty(testisopos_r); fractionation.all.isopos_r(jj,1) = testisopos_r;
                else;                        fractionation.all.isopos_r(jj,1) = NaN;
                end
                if   ~isempty(testtempos_r); fractionation.all.empos_r(jj,1) = testtempos_r;
                else;                        fractionation.all.empos_r(jj,1) = NaN;
                end
                if isequal(ObsList{activeloc},'d7Li');                                       fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'Li'));   rtest = find(ismember(ObsList_r,'Li'));   if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if d7Li, find Li
                if isequal(ObsList{activeloc},'d18O');                                       fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'SO4'));  rtest = find(ismember(ObsList_r,'SO4'));  if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if d18O, find SO4
                if isequal(ObsList{activeloc},'d26Mg');                                      fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'Mg'));   rtest = find(ismember(ObsList_r,'Mg'));   if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if d26Mg, find Mg
                if isequal(ObsList{activeloc},'d30Si');                                      fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'Si'));   rtest = find(ismember(ObsList_r,'Si'));   if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if d30Si, find Si
                if isequal(ObsList{activeloc},'d34S');                                       fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'SO4'));  rtest = find(ismember(ObsList_r,'SO4'));  if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if d34S, find SO4
                if isequal(ObsList{activeloc},'d42Ca');                                      fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'Ca'));   rtest = find(ismember(ObsList_r,'Ca'));   if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if d42Ca, find Ca
                if isequal(ObsList{activeloc},'d44Ca');                                      fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'Ca'));   rtest = find(ismember(ObsList_r,'Ca'));   if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if d44Ca, find Ca                
                if isequal(ObsList{activeloc},'d56Fe');                                      fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'Fe'));   rtest = find(ismember(ObsList_r,'Fe'));   if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if d56Fe, find Fe                
                if isequal(ObsList{activeloc},'Sr8786');                                     fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'Sr'));   rtest = find(ismember(ObsList_r,'Sr'));   if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if Sr87/86, find Sr
                if isequal(ObsList{activeloc},'Os8788');                                     fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'Os'));   rtest = find(ismember(ObsList_r,'Os'));   if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if Os8788, find Os
                if isequal(ObsList{activeloc},'d98Mo');                                      fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'Mo'));   rtest = find(ismember(ObsList_r,'Mo'));   if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if d98Mo, find Mo
                if isequal(ObsList{activeloc},'d13C') & isequal(carbonisotopematch,'HCO3');  fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'HCO3')); rtest = find(ismember(ObsList_r,'HCO3')); if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if d13C and carbon variable is HCO3, find HCO3
                if isequal(ObsList{activeloc},'Fmod') & isequal(carbonisotopematch,'HCO3');  fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'HCO3')); rtest = find(ismember(ObsList_r,'HCO3')); if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if Fmod and carbon variable is HCO3, find HCO3
                if isequal(ObsList{activeloc},'d13C') & isequal(carbonisotopematch,'DIC');   fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'DIC'));  rtest = find(ismember(ObsList_r,'DIC'));  if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if d13C and carbon variable is DIC, find DIC
                if isequal(ObsList{activeloc},'Fmod') & isequal(carbonisotopematch,'DIC');   fractionation.all.ionpos0(jj,1) = find(ismember(ObsList,'DIC'));  rtest = find(ismember(ObsList_r,'DIC'));  if ~isempty(rtest); fractionation.all.ionpos_r(jj,1) = rtest; else; fractionation.all.ionpos_r(jj,1) = NaN; end; end % if Fmod and carbon variable is DIC, find DIC

                 % check that fractionation has only been requested for isotopic data
                 if    sum(ismember(fractionation.all.isopos0,indxisopos))~=length(fractionation.all.isopos0)
                       problemobs = fractionation.all.iso(~ismember(fractionation.all.isopos0,indxisopos));
                       problemEMs = fractionation.all.EM(~ismember(fractionation.all.isopos0,indxisopos));
                       for j=1:length(problemobs)
                           disp(sprintf('warning! MEANDIR is attempting to apply fractionation to "%s" in "%s" (not to isotopes, to concentration). this warning was generated in MEANDIR_FindFractionationPairs.',problemobs{j},problemEMs{j}));
                       end
                 end   % end of checking for isotopic data
         end           % end of iterating over instances of fractionation
end                    % end of function