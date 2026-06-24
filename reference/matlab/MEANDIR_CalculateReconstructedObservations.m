% this script calculates the reconstruction of normalized river water chemistry. 
% for isotopic data, it calculates the actual inversion-constrained values. 

% (1) first define the structure 'recon' which will hold the river observations reconstructed from model results
clear recon; recon = struct;
if sum(ismember(ObsList,'ALK'))    > 0;  recon.ALK    = NaN(s,remaininginstances);  end  %  ALK
if sum(ismember(ObsList,'DIC'))    > 0;  recon.DIC    = NaN(s,remaininginstances);  end  %  DIC
if sum(ismember(ObsList,'Ca'))     > 0;  recon.Ca     = NaN(s,remaininginstances);  end  %  Ca
if sum(ismember(ObsList,'Mg'))     > 0;  recon.Mg     = NaN(s,remaininginstances);  end  %  Mg
if sum(ismember(ObsList,'Na'))     > 0;  recon.Na     = NaN(s,remaininginstances);  end  %  Na
if sum(ismember(ObsList,'K'))      > 0;  recon.K      = NaN(s,remaininginstances);  end  %  K
if sum(ismember(ObsList,'Sr'))     > 0;  recon.Sr     = NaN(s,remaininginstances);  end  %  Sr
if sum(ismember(ObsList,'Fe'))     > 0;  recon.Fe     = NaN(s,remaininginstances);  end  %  Fe
if sum(ismember(ObsList,'Cl'))     > 0;  recon.Cl     = NaN(s,remaininginstances);  end  %  Cl
if sum(ismember(ObsList,'SO4'))    > 0;  recon.SO4    = NaN(s,remaininginstances);  end  %  SO4
if sum(ismember(ObsList,'NO3'))    > 0;  recon.NO3    = NaN(s,remaininginstances);  end  %  NO3
if sum(ismember(ObsList,'PO4'))    > 0;  recon.PO4    = NaN(s,remaininginstances);  end  %  PO4
if sum(ismember(ObsList,'Si'))     > 0;  recon.Si     = NaN(s,remaininginstances);  end  %  Si
if sum(ismember(ObsList,'Ge'))     > 0;  recon.Ge     = NaN(s,remaininginstances);  end  %  Ge
if sum(ismember(ObsList,'Li'))     > 0;  recon.Li     = NaN(s,remaininginstances);  end  %  Li
if sum(ismember(ObsList,'F'))      > 0;  recon.F      = NaN(s,remaininginstances);  end  %  F
if sum(ismember(ObsList,'B'))      > 0;  recon.B      = NaN(s,remaininginstances);  end  %  B
if sum(ismember(ObsList,'Re'))     > 0;  recon.Re     = NaN(s,remaininginstances);  end  %  Re
if sum(ismember(ObsList,'Mo'))     > 0;  recon.Mo     = NaN(s,remaininginstances);  end  %  Mo
if sum(ismember(ObsList,'Os'))     > 0;  recon.Os     = NaN(s,remaininginstances);  end  %  Os
if sum(ismember(ObsList,'HCO3'))   > 0;  recon.HCO3   = NaN(s,remaininginstances);  end  %  HCO3
if sum(ismember(ObsList,'d7Li'))   > 0;  recon.d7Li   = NaN(s,remaininginstances);  end  %  d7Li
if sum(ismember(ObsList,'d18O'))   > 0;  recon.d18O   = NaN(s,remaininginstances);  end  %  d18O
if sum(ismember(ObsList,'d26Mg'))  > 0;  recon.d26Mg  = NaN(s,remaininginstances);  end  %  d26Mg
if sum(ismember(ObsList,'d30Si'))  > 0;  recon.d30Si  = NaN(s,remaininginstances);  end  %  d30Si
if sum(ismember(ObsList,'d34S'))   > 0;  recon.d34S   = NaN(s,remaininginstances);  end  %  d34S
if sum(ismember(ObsList,'d42Ca'))  > 0;  recon.d42Ca  = NaN(s,remaininginstances);  end  %  d42Ca
if sum(ismember(ObsList,'d44Ca'))  > 0;  recon.d44Ca  = NaN(s,remaininginstances);  end  %  d44Ca
if sum(ismember(ObsList,'d56Fe'))  > 0;  recon.d56Fe  = NaN(s,remaininginstances);  end  %  d56Fe
if sum(ismember(ObsList,'Sr8786')) > 0;  recon.Sr8786 = NaN(s,remaininginstances);  end  %  Sr8786
if sum(ismember(ObsList,'d98Mo'))  > 0;  recon.d98Mo  = NaN(s,remaininginstances);  end  %  d98Mo
if sum(ismember(ObsList,'Os8788')) > 0;  recon.Os8788 = NaN(s,remaininginstances);  end  %  Os8788
if sum(ismember(ObsList,'d13C'))   > 0;  recon.d13C   = NaN(s,remaininginstances);  end  %  d13C
if sum(ismember(ObsList,'Fmod'))   > 0;  recon.Fmod   = NaN(s,remaininginstances);  end  %  Fmod

% (2) iterate over each sample and iteration. calculate MEANDIR's
% prediction of the river chemistry and save that into the field recon.
for i=1:s
for j=1:remaininginstances      
    reconriver = EMSucessIterMatrixAll(:,:,j,i)*results.norm(:,j,i);
    for ll = 1:length(ObsList)    
        evalin('base',[sprintf('recon.%s(i,j)',ObsList{ll}) '= reconriver(ll);' ]);
    end
end    
end

% (3a) for isotope data, scale the data by normalized abundance to recover 
% either the delta value or the raw isotope ratio. 
if sum(ismember(ObsList,'d7Li'))  >0;                                      recon.d7Li   = recon.d7Li./recon.Li;   end  %  d7Li
if sum(ismember(ObsList,'d18O'))  >0;                                      recon.d18O   = recon.d18O./recon.SO4;  end  %  d18O
if sum(ismember(ObsList,'d26Mg')) >0;                                      recon.d26Mg  = recon.d26Mg./recon.Mg;  end  %  d26Mg
if sum(ismember(ObsList,'d30Si')) >0;                                      recon.d30Si  = recon.d30Si./recon.Si;  end  %  d30Si
if sum(ismember(ObsList,'d34S'))  >0;                                      recon.d34S   = recon.d34S./recon.SO4;  end  %  d34S
if sum(ismember(ObsList,'d42Ca')) >0;                                      recon.d42Ca  = recon.d42Ca./recon.Ca;  end  %  d42Ca
if sum(ismember(ObsList,'d44Ca')) >0;                                      recon.d44Ca  = recon.d44Ca./recon.Ca;  end  %  d44Ca
if sum(ismember(ObsList,'d56Fe')) >0;                                      recon.d56Fe  = recon.d56Fe./recon.Fe;  end  %  d56Fe
if sum(ismember(ObsList,'Sr8786'))>0;                                      recon.Sr8786 = recon.Sr8786./recon.Sr; end  %  Sr8786
if sum(ismember(ObsList,'d98Mo')) >0;                                      recon.d98Mo  = recon.d98Mo./recon.Mo;  end  %  d98Mo
if sum(ismember(ObsList,'Os8788'))>0;                                      recon.Os8788 = recon.Os8788./recon.Os; end  %  Os8788
if sum(ismember(ObsList,'d13C'))  >0 & isequal(carbonisotopematch,'HCO3'); recon.d13C   = recon.d13C./recon.HCO3; end  %  d13C for HCO3
if sum(ismember(ObsList,'Fmod'))  >0 & isequal(carbonisotopematch,'HCO3'); recon.Fmod   = recon.Fmod./recon.HCO3; end  %  Fmod for HCO3
if sum(ismember(ObsList,'d13C'))  >0 & isequal(carbonisotopematch,'DIC');  recon.d13C   = recon.d13C./recon.DIC;  end  %  d13C for DIC
if sum(ismember(ObsList,'Fmod'))  >0 & isequal(carbonisotopematch,'DIC');  recon.Fmod   = recon.Fmod./recon.DIC;  end  %  Fmod for DIC
    
% (3b) for isotope data, scale isotopic ratios back to delta values (when they were converted within MEANDIR)
if sum(ismember(ObsList,'d7Li'))  >0 & sum(ismember(ConvertDelta2RList,'d7Li')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d7Li'));  recon.d7Li  = (recon.d7Li/conv-1)*1000;  end  %  d7Li
if sum(ismember(ObsList,'d18O'))  >0 & sum(ismember(ConvertDelta2RList,'d18O')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d18O'));  recon.d18O  = (recon.d18O/conv-1)*1000;  end  %  d18O
if sum(ismember(ObsList,'d26Mg')) >0 & sum(ismember(ConvertDelta2RList,'d26Mg'))>0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d26Mg')); recon.d26Mg = (recon.d26Mg/conv-1)*1000; end  %  d26Mg
if sum(ismember(ObsList,'d30Si')) >0 & sum(ismember(ConvertDelta2RList,'d30Si'))>0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d30Si')); recon.d30Si = (recon.d30Si/conv-1)*1000; end  %  d30Si
if sum(ismember(ObsList,'d34S'))  >0 & sum(ismember(ConvertDelta2RList,'d34S')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d34S'));  recon.d34S  = (recon.d34S/conv-1)*1000;  end  %  d34S
if sum(ismember(ObsList,'d42Ca')) >0 & sum(ismember(ConvertDelta2RList,'d42Ca'))>0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d42Ca')); recon.d42Ca = (recon.d42Ca/conv-1)*1000; end  %  d42Ca
if sum(ismember(ObsList,'d44Ca')) >0 & sum(ismember(ConvertDelta2RList,'d44Ca'))>0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d44Ca')); recon.d44Ca = (recon.d44Ca/conv-1)*1000; end  %  d44Ca
if sum(ismember(ObsList,'d56Fe')) >0 & sum(ismember(ConvertDelta2RList,'d56Fe'))>0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d56Fe')); recon.d56Fe = (recon.d56Fe/conv-1)*1000; end  %  d56Fe
if sum(ismember(ObsList,'d98Mo')) >0 & sum(ismember(ConvertDelta2RList,'d98Mo'))>0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d98Mo')); recon.d98Mo = (recon.d98Mo/conv-1)*1000; end  %  d98Mo
if sum(ismember(ObsList,'d13C'))  >0 & sum(ismember(ConvertDelta2RList,'d13C')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d13C'));  recon.d13C  = (recon.d13C/conv-1)*1000;  end  %  d13C

% (4) save the data into the 'river' structure
for ll = 1:length(ObsList)      
    saveObs = ObsList; 
    isotopeindx = find(isotopeposition);
    for llname = 1:length(saveObs) 
        if sum(isotopeindx==llname)==0
            saveObs{llname} = cat(2,saveObs{llname} ,'_norm');
        end
    end    
    evalin('base',[sprintf('river.reconstructed.%s.all',   saveObs{ll}) '=' sprintf('recon.%s;',ObsList{ll})]);     
    evalin('base',[sprintf('river.reconstructed.%s.median',saveObs{ll}) '=' sprintf('nanmedian(recon.%s,2);',ObsList{ll})]);  
    evalin('base',[sprintf('river.reconstructed.%s.mean',  saveObs{ll}) '=' sprintf('nanmean(recon.%s,2);',ObsList{ll})]);       
    evalin('base',[sprintf('river.reconstructed.%s.pct05', saveObs{ll}) '=' sprintf('prctile(recon.%s%s,05)%s;',ObsList{ll},'''','''')]);  
    evalin('base',[sprintf('river.reconstructed.%s.pct25', saveObs{ll}) '=' sprintf('prctile(recon.%s%s,25)%s;',ObsList{ll},'''','''')]);  
    evalin('base',[sprintf('river.reconstructed.%s.pct75', saveObs{ll}) '=' sprintf('prctile(recon.%s%s,75)%s;',ObsList{ll},'''','''')]);  
    evalin('base',[sprintf('river.reconstructed.%s.pct95', saveObs{ll}) '=' sprintf('prctile(recon.%s%s,95)%s;',ObsList{ll},'''','''')]);  
end

% (5) clear excess variables (keep the workspace neat)
clear recon;