function IterEM_smpfrac = MEANDIR_updateIterEM_smpfrac(smpres,IterEM_smpfrac,i, go)

         % this function saves the results in "smpres" into the structure variable IterEM_smpfrac
         % this function is required for removing results from the parfor loop without generating transparency or classification errors

                            IterEM_smpfrac.smpresultsnorm(:,i)   = smpres.smpresultsnorm;       % normalization
                            IterEM_smpfrac.smpresultsSO4(:,i)    = smpres.smpresultsSO4;        % SO4 
                            IterEM_smpfrac.smpresultsd34S(:,i)   = smpres.smpresultsd34S;       % d34S 
         if go.ALK    == 1; IterEM_smpfrac.smpresultsALK(:,i)    = smpres.smpresultsALK;    end % ALK
         if go.DIC    == 1; IterEM_smpfrac.smpresultsDIC(:,i)    = smpres.smpresultsDIC;    end % DIC
         if go.Ca     == 1; IterEM_smpfrac.smpresultsCa(:,i)     = smpres.smpresultsCa;     end % Ca 
         if go.Mg     == 1; IterEM_smpfrac.smpresultsMg(:,i)     = smpres.smpresultsMg;     end % Mg
         if go.Na     == 1; IterEM_smpfrac.smpresultsNa(:,i)     = smpres.smpresultsNa;     end % Na
         if go.K      == 1; IterEM_smpfrac.smpresultsK(:,i)      = smpres.smpresultsK;      end % K
         if go.Sr     == 1; IterEM_smpfrac.smpresultsSr(:,i)     = smpres.smpresultsSr;     end % Sr
         if go.Fe     == 1; IterEM_smpfrac.smpresultsFe(:,i)     = smpres.smpresultsFe;     end % Fe
         if go.Cl     == 1; IterEM_smpfrac.smpresultsCl(:,i)     = smpres.smpresultsCl;     end % Cl
         if go.NO3    == 1; IterEM_smpfrac.smpresultsNO3(:,i)    = smpres.smpresultsNO3;    end % NO3 
         if go.PO4    == 1; IterEM_smpfrac.smpresultsPO4(:,i)    = smpres.smpresultsPO4;    end % PO4 
         if go.Si     == 1; IterEM_smpfrac.smpresultsSi(:,i)     = smpres.smpresultsSi;     end % Si
         if go.Ge     == 1; IterEM_smpfrac.smpresultsGe(:,i)     = smpres.smpresultsGe;     end % Ge
         if go.Li     == 1; IterEM_smpfrac.smpresultsLi(:,i)     = smpres.smpresultsLi;     end % Li
         if go.F      == 1; IterEM_smpfrac.smpresultsF(:,i)      = smpres.smpresultsF;      end % F
         if go.B      == 1; IterEM_smpfrac.smpresultsB(:,i)      = smpres.smpresultsB;      end % B
         if go.Re     == 1; IterEM_smpfrac.smpresultsRe(:,i)     = smpres.smpresultsRe;     end % Re
         if go.Mo     == 1; IterEM_smpfrac.smpresultsMo(:,i)     = smpres.smpresultsMo;     end % Mo
         if go.Os     == 1; IterEM_smpfrac.smpresultsOs(:,i)     = smpres.smpresultsOs;     end % Os
         if go.HCO3   == 1; IterEM_smpfrac.smpresultsHCO3(:,i)   = smpres.smpresultsHCO3;   end % HCO3
         if go.d7Li   == 1; IterEM_smpfrac.smpresultsd7Li(:,i)   = smpres.smpresultsd7Li;   end % d7Li
         if go.d13C   == 1; IterEM_smpfrac.smpresultsd13C(:,i)   = smpres.smpresultsd13C;   end % d13C
         if go.d18O   == 1; IterEM_smpfrac.smpresultsd18O(:,i)   = smpres.smpresultsd18O;   end % d18O
         if go.d26Mg  == 1; IterEM_smpfrac.smpresultsd26Mg(:,i)  = smpres.smpresultsd26Mg;  end % d26Mg
         if go.d30Si  == 1; IterEM_smpfrac.smpresultsd30Si(:,i)  = smpres.smpresultsd30Si;  end % d30Si
         if go.d42Ca  == 1; IterEM_smpfrac.smpresultsd42Ca(:,i)  = smpres.smpresultsd42Ca;  end % d42Ca
         if go.d44Ca  == 1; IterEM_smpfrac.smpresultsd44Ca(:,i)  = smpres.smpresultsd44Ca;  end % d44Ca
         if go.d56Fe  == 1; IterEM_smpfrac.smpresultsd56Fe(:,i)  = smpres.smpresultsd56Fe;  end % d56Fe
         if go.Sr8786 == 1; IterEM_smpfrac.smpresultsSr8786(:,i) = smpres.smpresultsSr8786; end % 87Sr/86Sr
         if go.d98Mo  == 1; IterEM_smpfrac.smpresultsd98Mo(:,i)  = smpres.smpresultsd98Mo;  end % d98Mo
         if go.Os8788 == 1; IterEM_smpfrac.smpresultsOs8788(:,i) = smpres.smpresultsOs8788; end % 188Os
         if go.Fmod   == 1; IterEM_smpfrac.smpresultsFmod(:,i)   = smpres.smpresultsFmod;   end % Fmod

end % end of function