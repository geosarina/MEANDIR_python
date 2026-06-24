function [IterEM_smpfrac, IterEM_InvRivData] = MEANDIR_DeclareInternalLoopIterEMStructures(go, nEM, s)

         % this script defines two structure variables to be populated with
         % inversion results when Iterover equals "End-members"
         %   IterEM_smpfrac stores sets of fractional contributions
         %   IterEM_InvRivData stores the actual river observations used in the inversion

                             IterEM_smpfrac = struct;                       IterEM_InvRivData        = struct;
                             IterEM_smpfrac.smpresultsnorm   = NaN(nEM,s);  IterEM_InvRivData.norm   = NaN(s,1);       %  norm
                             IterEM_smpfrac.smpresultsSO4    = NaN(nEM,s);  IterEM_InvRivData.SO4    = NaN(s,1);       %  SO4 
                             IterEM_smpfrac.smpresultsd34S   = NaN(nEM,s);  IterEM_InvRivData.d34S   = NaN(s,1);       %  d34S
         if go.ALK    == 1;  IterEM_smpfrac.smpresultsALK    = NaN(nEM,s);  IterEM_InvRivData.ALK    = NaN(s,1);  end  %  ALK
         if go.DIC    == 1;  IterEM_smpfrac.smpresultsDIC    = NaN(nEM,s);  IterEM_InvRivData.DIC    = NaN(s,1);  end  %  DIC
         if go.Ca     == 1;  IterEM_smpfrac.smpresultsCa     = NaN(nEM,s);  IterEM_InvRivData.Ca     = NaN(s,1);  end  %  Ca
         if go.Mg     == 1;  IterEM_smpfrac.smpresultsMg     = NaN(nEM,s);  IterEM_InvRivData.Mg     = NaN(s,1);  end  %  Mg
         if go.Na     == 1;  IterEM_smpfrac.smpresultsNa     = NaN(nEM,s);  IterEM_InvRivData.Na     = NaN(s,1);  end  %  Na
         if go.K      == 1;  IterEM_smpfrac.smpresultsK      = NaN(nEM,s);  IterEM_InvRivData.K      = NaN(s,1);  end  %  K
         if go.Sr     == 1;  IterEM_smpfrac.smpresultsSr     = NaN(nEM,s);  IterEM_InvRivData.Sr     = NaN(s,1);  end  %  Sr
         if go.Fe     == 1;  IterEM_smpfrac.smpresultsFe     = NaN(nEM,s);  IterEM_InvRivData.Fe     = NaN(s,1);  end  %  Fe
         if go.Cl     == 1;  IterEM_smpfrac.smpresultsCl     = NaN(nEM,s);  IterEM_InvRivData.Cl     = NaN(s,1);  end  %  Cl
         if go.NO3    == 1;  IterEM_smpfrac.smpresultsNO3    = NaN(nEM,s);  IterEM_InvRivData.NO3    = NaN(s,1);  end  %  NO3
         if go.PO4    == 1;  IterEM_smpfrac.smpresultsPO4    = NaN(nEM,s);  IterEM_InvRivData.PO4    = NaN(s,1);  end  %  PO4
         if go.Si     == 1;  IterEM_smpfrac.smpresultsSi     = NaN(nEM,s);  IterEM_InvRivData.Si     = NaN(s,1);  end  %  Si
         if go.Ge     == 1;  IterEM_smpfrac.smpresultsGe     = NaN(nEM,s);  IterEM_InvRivData.Ge     = NaN(s,1);  end  %  Ge
         if go.Li     == 1;  IterEM_smpfrac.smpresultsLi     = NaN(nEM,s);  IterEM_InvRivData.Li     = NaN(s,1);  end  %  Li
         if go.F      == 1;  IterEM_smpfrac.smpresultsF      = NaN(nEM,s);  IterEM_InvRivData.F      = NaN(s,1);  end  %  F
         if go.B      == 1;  IterEM_smpfrac.smpresultsB      = NaN(nEM,s);  IterEM_InvRivData.B      = NaN(s,1);  end  %  B
         if go.Re     == 1;  IterEM_smpfrac.smpresultsRe     = NaN(nEM,s);  IterEM_InvRivData.Re     = NaN(s,1);  end  %  Re
         if go.Mo     == 1;  IterEM_smpfrac.smpresultsMo     = NaN(nEM,s);  IterEM_InvRivData.Mo     = NaN(s,1);  end  %  Mo
         if go.Os     == 1;  IterEM_smpfrac.smpresultsOs     = NaN(nEM,s);  IterEM_InvRivData.Os     = NaN(s,1);  end  %  Os
         if go.HCO3   == 1;  IterEM_smpfrac.smpresultsHCO3   = NaN(nEM,s);  IterEM_InvRivData.HCO3   = NaN(s,1);  end  %  HCO3
         if go.d7Li   == 1;  IterEM_smpfrac.smpresultsd7Li   = NaN(nEM,s);  IterEM_InvRivData.d7Li   = NaN(s,1);  end  %  d7Li
         if go.d13C   == 1;  IterEM_smpfrac.smpresultsd13C   = NaN(nEM,s);  IterEM_InvRivData.d13C   = NaN(s,1);  end  %  d13C
         if go.d18O   == 1;  IterEM_smpfrac.smpresultsd18O   = NaN(nEM,s);  IterEM_InvRivData.d18O   = NaN(s,1);  end  %  d18O
         if go.d26Mg  == 1;  IterEM_smpfrac.smpresultsd26Mg  = NaN(nEM,s);  IterEM_InvRivData.d26Mg  = NaN(s,1);  end  %  d26Mg
         if go.d30Si  == 1;  IterEM_smpfrac.smpresultsd30Si  = NaN(nEM,s);  IterEM_InvRivData.d30Si  = NaN(s,1);  end  %  d30Si
         if go.d42Ca  == 1;  IterEM_smpfrac.smpresultsd42Ca  = NaN(nEM,s);  IterEM_InvRivData.d42Ca  = NaN(s,1);  end  %  d42Ca
         if go.d44Ca  == 1;  IterEM_smpfrac.smpresultsd44Ca  = NaN(nEM,s);  IterEM_InvRivData.d44Ca  = NaN(s,1);  end  %  d44Ca
         if go.d56Fe  == 1;  IterEM_smpfrac.smpresultsd56Fe  = NaN(nEM,s);  IterEM_InvRivData.d56Fe  = NaN(s,1);  end  %  d56Fe
         if go.Sr8786 == 1;  IterEM_smpfrac.smpresultsSr8786 = NaN(nEM,s);  IterEM_InvRivData.Sr8786 = NaN(s,1);  end  %  Sr8786
         if go.d98Mo  == 1;  IterEM_smpfrac.smpresultsd98Mo  = NaN(nEM,s);  IterEM_InvRivData.d98Mo  = NaN(s,1);  end  %  d98Mo
         if go.Os8788 == 1;  IterEM_smpfrac.smpresultsOs8788 = NaN(nEM,s);  IterEM_InvRivData.Os8788 = NaN(s,1);  end  %  Os8788
         if go.Fmod   == 1;  IterEM_smpfrac.smpresultsFmod   = NaN(nEM,s);  IterEM_InvRivData.Fmod   = NaN(s,1);  end  %  Fmod

end      % end of function