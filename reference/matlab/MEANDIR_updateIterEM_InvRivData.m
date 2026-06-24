function IterEM_InvRivData = MEANDIR_updateIterEM_InvRivData(InvRivDatsave,IterEM_InvRivData,i, go)
       
         % this function save the results in InvRivDatsave into the structure variable IterEM_InvRivData
         % this function is required for removing results from the parfor loop without generating transparency or classification errors
         % this function is called only when 'IterateOver' equals 'End-members'

                            IterEM_InvRivData.norm(i,1)   = InvRivDatsave.norm;       % norm
                            IterEM_InvRivData.SO4(i,1)    = InvRivDatsave.SO4;        % SO4 
                            IterEM_InvRivData.d34S(i,1)   = InvRivDatsave.d34S;       % d34S
         if go.ALK    == 1; IterEM_InvRivData.ALK(i,1)    = InvRivDatsave.ALK;    end % ALK
         if go.DIC    == 1; IterEM_InvRivData.DIC(i,1)    = InvRivDatsave.DIC;    end % DIC
         if go.Ca     == 1; IterEM_InvRivData.Ca(i,1)     = InvRivDatsave.Ca;     end % Ca
         if go.Mg     == 1; IterEM_InvRivData.Mg(i,1)     = InvRivDatsave.Mg;     end % Mg
         if go.Na     == 1; IterEM_InvRivData.Na(i,1)     = InvRivDatsave.Na;     end % Na
         if go.K      == 1; IterEM_InvRivData.K(i,1)      = InvRivDatsave.K;      end % K
         if go.Sr     == 1; IterEM_InvRivData.Sr(i,1)     = InvRivDatsave.Sr;     end % Sr
         if go.Fe     == 1; IterEM_InvRivData.Fe(i,1)     = InvRivDatsave.Fe;     end % Fe
         if go.Cl     == 1; IterEM_InvRivData.Cl(i,1)     = InvRivDatsave.Cl;     end % Cl
         if go.NO3    == 1; IterEM_InvRivData.NO3(i,1)    = InvRivDatsave.NO3;    end % NO3
         if go.PO4    == 1; IterEM_InvRivData.PO4(i,1)    = InvRivDatsave.PO4;    end % PO4
         if go.Si     == 1; IterEM_InvRivData.Si(i,1)     = InvRivDatsave.Si;     end % Si
         if go.Ge     == 1; IterEM_InvRivData.Ge(i,1)     = InvRivDatsave.Ge;     end % Ge
         if go.Li     == 1; IterEM_InvRivData.Li(i,1)     = InvRivDatsave.Li;     end % Li
         if go.F      == 1; IterEM_InvRivData.F(i,1)      = InvRivDatsave.F;      end % F
         if go.B      == 1; IterEM_InvRivData.B(i,1)      = InvRivDatsave.B;      end % B
         if go.Re     == 1; IterEM_InvRivData.Re(i,1)     = InvRivDatsave.Re;     end % Re
         if go.Mo     == 1; IterEM_InvRivData.Mo(i,1)     = InvRivDatsave.Mo;     end % Mo 
         if go.Os     == 1; IterEM_InvRivData.Os(i,1)     = InvRivDatsave.Os;     end % Os
         if go.HCO3   == 1; IterEM_InvRivData.HCO3(i,1)   = InvRivDatsave.HCO3;   end % HCO3 
         if go.d7Li   == 1; IterEM_InvRivData.d7Li(i,1)   = InvRivDatsave.d7Li;   end % d7Li
         if go.d13C   == 1; IterEM_InvRivData.d13C(i,1)   = InvRivDatsave.d13C;   end % d13C
         if go.d18O   == 1; IterEM_InvRivData.d18O(i,1)   = InvRivDatsave.d18O;   end % d18O
         if go.d26Mg  == 1; IterEM_InvRivData.d26Mg(i,1)  = InvRivDatsave.d26Mg;  end % d26Mg
         if go.d30Si  == 1; IterEM_InvRivData.d30Si(i,1)  = InvRivDatsave.d30Si;  end % d30Si
         if go.d42Ca  == 1; IterEM_InvRivData.d42Ca(i,1)  = InvRivDatsave.d42Ca;  end % d42Ca
         if go.d44Ca  == 1; IterEM_InvRivData.d44Ca(i,1)  = InvRivDatsave.d44Ca;  end % d44Ca
         if go.d56Fe  == 1; IterEM_InvRivData.d56Fe(i,1)  = InvRivDatsave.d56Fe;  end % d56Fe
         if go.Sr8786 == 1; IterEM_InvRivData.Sr8786(i,1) = InvRivDatsave.Sr8786; end % Sr8786
         if go.d98Mo  == 1; IterEM_InvRivData.d98Mo(i,1)  = InvRivDatsave.d98Mo;  end % d98Mo
         if go.Os8788 == 1; IterEM_InvRivData.Os8788(i,1) = InvRivDatsave.Os8788; end % Os8788
         if go.Fmod   == 1; IterEM_InvRivData.Fmod(i,1)   = InvRivDatsave.Fmod;   end % Fmod
         
end % end of function