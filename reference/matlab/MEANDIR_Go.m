function go = MEANDIR_Go(ObsList)

         % this function evaluates which variables are included in "ObsList" and
         % in the inversion. if a variable is included, the corresponding field
         % of the structure "go" is given a value of 1. if the variable is not included,
         % the corresponding field of the structure "go" is given a value of 0.

         go = struct; % define the structure go
         if sum(ismember(ObsList,'ALK'))   >0; go.ALK    = 1; else; go.ALK    = 0; end % check if ALK is included in the inversion
         if sum(ismember(ObsList,'DIC'))   >0; go.DIC    = 1; else; go.DIC    = 0; end % check if DIC is included in the inversion
         if sum(ismember(ObsList,'Ca'))    >0; go.Ca     = 1; else; go.Ca     = 0; end % check if Ca is included in the inversion
         if sum(ismember(ObsList,'Mg'))    >0; go.Mg     = 1; else; go.Mg     = 0; end % check if Mg is included in the inversion
         if sum(ismember(ObsList,'Na'))    >0; go.Na     = 1; else; go.Na     = 0; end % check if Na is included in the inversion
         if sum(ismember(ObsList,'K'))     >0; go.K      = 1; else; go.K      = 0; end % check if K is included in the inversion
         if sum(ismember(ObsList,'Sr'))    >0; go.Sr     = 1; else; go.Sr     = 0; end % check if Sr is included in the inversion
         if sum(ismember(ObsList,'Fe'))    >0; go.Fe     = 1; else; go.Fe     = 0; end % check if Fe is included in the inversion
         if sum(ismember(ObsList,'Cl'))    >0; go.Cl     = 1; else; go.Cl     = 0; end % check if Cl is included in the inversion
         if sum(ismember(ObsList,'SO4'))   >0; go.SO4    = 1; else; go.SO4    = 0; end % check if SO4 is included in the inversion
         if sum(ismember(ObsList,'NO3'))   >0; go.NO3    = 1; else; go.NO3    = 0; end % check if NO3 is included in the inversion
         if sum(ismember(ObsList,'PO4'))   >0; go.PO4    = 1; else; go.PO4    = 0; end % check if PO4 is included in the inversion
         if sum(ismember(ObsList,'Si'))    >0; go.Si     = 1; else; go.Si     = 0; end % check if Si is included in the inversion
         if sum(ismember(ObsList,'Ge'))    >0; go.Ge     = 1; else; go.Ge     = 0; end % check if Ge is included in the inversion
         if sum(ismember(ObsList,'Li'))    >0; go.Li     = 1; else; go.Li     = 0; end % check if Li is included in the inversion
         if sum(ismember(ObsList,'F'))     >0; go.F      = 1; else; go.F      = 0; end % check if F is included in the inversion
         if sum(ismember(ObsList,'B'))     >0; go.B      = 1; else; go.B      = 0; end % check if B is included in the inversion
         if sum(ismember(ObsList,'Re'))    >0; go.Re     = 1; else; go.Re     = 0; end % check if Re is included in the inversion
         if sum(ismember(ObsList,'Mo'))    >0; go.Mo     = 1; else; go.Mo     = 0; end % check if Mo is included in the inversion
         if sum(ismember(ObsList,'Os'))    >0; go.Os     = 1; else; go.Os     = 0; end % check if Os is included in the inversion
         if sum(ismember(ObsList,'HCO3'))  >0; go.HCO3   = 1; else; go.HCO3   = 0; end % check if HCO3 is included in the inversion
         if sum(ismember(ObsList,'d7Li'))  >0; go.d7Li   = 1; else; go.d7Li   = 0; end % check if d7Li is included in the inversion
         if sum(ismember(ObsList,'d13C'))  >0; go.d13C   = 1; else; go.d13C   = 0; end % check if d13C is included in the inversion
         if sum(ismember(ObsList,'d18O'))  >0; go.d18O   = 1; else; go.d18O   = 0; end % check if d18O is included in the inversion
         if sum(ismember(ObsList,'d26Mg')) >0; go.d26Mg  = 1; else; go.d26Mg  = 0; end % check if d26Mg is included in the inversion
         if sum(ismember(ObsList,'d30Si')) >0; go.d30Si  = 1; else; go.d30Si  = 0; end % check if d30Si is included in the inversion
         if sum(ismember(ObsList,'d34S'))  >0; go.d34S   = 1; else; go.d34S   = 0; end % check if d34S is included in the inversion
         if sum(ismember(ObsList,'d42Ca')) >0; go.d42Ca  = 1; else; go.d42Ca  = 0; end % check if d42Ca is included in the inversion
         if sum(ismember(ObsList,'d44Ca')) >0; go.d44Ca  = 1; else; go.d44Ca  = 0; end % check if d44Ca is included in the inversion
         if sum(ismember(ObsList,'d56Fe')) >0; go.d56Fe  = 1; else; go.d56Fe  = 0; end % check if d56Fe is included in the inversion
         if sum(ismember(ObsList,'Sr8786'))>0; go.Sr8786 = 1; else; go.Sr8786 = 0; end % check if 87Sr/86Sr is included in the inversion
         if sum(ismember(ObsList,'d98Mo')) >0; go.d98Mo  = 1; else; go.d98Mo  = 0; end % check if d98Mo is included in the inversion
         if sum(ismember(ObsList,'Os8788'))>0; go.Os8788 = 1; else; go.Os8788 = 0; end % check if 187Os/188Os is included in the inversion
         if sum(ismember(ObsList,'Fmod'))  >0; go.Fmod   = 1; else; go.Fmod   = 0; end % check if Fmod is included in the inversion

end % end of function