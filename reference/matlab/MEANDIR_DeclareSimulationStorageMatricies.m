function [InvRivDatsave, InvRivDat, smpres] = MEANDIR_DeclareSimulationStorageMatricies(maxsuccess, nEM, go)

          % this file defines matricies that will be used to store inversion results
          % "InvRivDatsave" will store information on river chemistry
          % "InvRivDat" will also store information on river chemistry, but only one value per subfield
          % "smpres" will store results for the fractional contributions of end-members to the normalization variable

          InvRivDatsave = struct;
          InvRivDat     = struct; 
          smpres        = struct;

                             InvRivDatsave.norm   = NaN(maxsuccess,1); InvRivDat.norm   = NaN;  smpres.smpresultsnorm   = NaN(nEM,maxsuccess);     % normalization
                             InvRivDatsave.SO4    = NaN(maxsuccess,1); InvRivDat.SO4    = NaN;  smpres.smpresultsSO4    = NaN(nEM,maxsuccess);     % SO4
                             InvRivDatsave.d34S   = NaN(maxsuccess,1); InvRivDat.d34S   = NaN;  smpres.smpresultsd34S   = NaN(nEM,maxsuccess);     % d34S
          if go.ALK    == 1; InvRivDatsave.ALK    = NaN(maxsuccess,1); InvRivDat.ALK    = NaN;  smpres.smpresultsALK    = NaN(nEM,maxsuccess); end % ALK
          if go.DIC    == 1; InvRivDatsave.DIC    = NaN(maxsuccess,1); InvRivDat.DIC    = NaN;  smpres.smpresultsDIC    = NaN(nEM,maxsuccess); end % DIC
          if go.Ca     == 1; InvRivDatsave.Ca     = NaN(maxsuccess,1); InvRivDat.Ca     = NaN;  smpres.smpresultsCa     = NaN(nEM,maxsuccess); end % Ca
          if go.Mg     == 1; InvRivDatsave.Mg     = NaN(maxsuccess,1); InvRivDat.Mg     = NaN;  smpres.smpresultsMg     = NaN(nEM,maxsuccess); end % Mg
          if go.Na     == 1; InvRivDatsave.Na     = NaN(maxsuccess,1); InvRivDat.Na     = NaN;  smpres.smpresultsNa     = NaN(nEM,maxsuccess); end % Na
          if go.K      == 1; InvRivDatsave.K      = NaN(maxsuccess,1); InvRivDat.K      = NaN;  smpres.smpresultsK      = NaN(nEM,maxsuccess); end % K
          if go.Sr     == 1; InvRivDatsave.Sr     = NaN(maxsuccess,1); InvRivDat.Sr     = NaN;  smpres.smpresultsSr     = NaN(nEM,maxsuccess); end % Sr
          if go.Fe     == 1; InvRivDatsave.Fe     = NaN(maxsuccess,1); InvRivDat.Fe     = NaN;  smpres.smpresultsFe     = NaN(nEM,maxsuccess); end % Fe
          if go.Cl     == 1; InvRivDatsave.Cl     = NaN(maxsuccess,1); InvRivDat.Cl     = NaN;  smpres.smpresultsCl     = NaN(nEM,maxsuccess); end % Cl
          if go.NO3    == 1; InvRivDatsave.NO3    = NaN(maxsuccess,1); InvRivDat.NO3    = NaN;  smpres.smpresultsNO3    = NaN(nEM,maxsuccess); end % NO3
          if go.PO4    == 1; InvRivDatsave.PO4    = NaN(maxsuccess,1); InvRivDat.PO4    = NaN;  smpres.smpresultsPO4    = NaN(nEM,maxsuccess); end % PO4
          if go.Si     == 1; InvRivDatsave.Si     = NaN(maxsuccess,1); InvRivDat.Si     = NaN;  smpres.smpresultsSi     = NaN(nEM,maxsuccess); end % Si
          if go.Ge     == 1; InvRivDatsave.Ge     = NaN(maxsuccess,1); InvRivDat.Ge     = NaN;  smpres.smpresultsGe     = NaN(nEM,maxsuccess); end % Ge
          if go.Li     == 1; InvRivDatsave.Li     = NaN(maxsuccess,1); InvRivDat.Li     = NaN;  smpres.smpresultsLi     = NaN(nEM,maxsuccess); end % Li
          if go.F      == 1; InvRivDatsave.F      = NaN(maxsuccess,1); InvRivDat.F      = NaN;  smpres.smpresultsF      = NaN(nEM,maxsuccess); end % F
          if go.B      == 1; InvRivDatsave.B      = NaN(maxsuccess,1); InvRivDat.B      = NaN;  smpres.smpresultsB      = NaN(nEM,maxsuccess); end % B
          if go.Re     == 1; InvRivDatsave.Re     = NaN(maxsuccess,1); InvRivDat.Re     = NaN;  smpres.smpresultsRe     = NaN(nEM,maxsuccess); end % Re
          if go.Mo     == 1; InvRivDatsave.Mo     = NaN(maxsuccess,1); InvRivDat.Mo     = NaN;  smpres.smpresultsMo     = NaN(nEM,maxsuccess); end % Mo
          if go.Os     == 1; InvRivDatsave.Os     = NaN(maxsuccess,1); InvRivDat.Os     = NaN;  smpres.smpresultsOs     = NaN(nEM,maxsuccess); end % Os
          if go.HCO3   == 1; InvRivDatsave.HCO3   = NaN(maxsuccess,1); InvRivDat.HCO3   = NaN;  smpres.smpresultsHCO3   = NaN(nEM,maxsuccess); end % HCO3
          if go.d7Li   == 1; InvRivDatsave.d7Li   = NaN(maxsuccess,1); InvRivDat.d7Li   = NaN;  smpres.smpresultsd7Li   = NaN(nEM,maxsuccess); end % d7Li
          if go.d13C   == 1; InvRivDatsave.d13C   = NaN(maxsuccess,1); InvRivDat.d13C   = NaN;  smpres.smpresultsd13C   = NaN(nEM,maxsuccess); end % d13C
          if go.d18O   == 1; InvRivDatsave.d18O   = NaN(maxsuccess,1); InvRivDat.d18O   = NaN;  smpres.smpresultsd18O   = NaN(nEM,maxsuccess); end % d18O
          if go.d26Mg  == 1; InvRivDatsave.d26Mg  = NaN(maxsuccess,1); InvRivDat.d26Mg  = NaN;  smpres.smpresultsd26Mg  = NaN(nEM,maxsuccess); end % d26Mg
          if go.d30Si  == 1; InvRivDatsave.d30Si  = NaN(maxsuccess,1); InvRivDat.d30Si  = NaN;  smpres.smpresultsd30Si  = NaN(nEM,maxsuccess); end % d30Si
          if go.d42Ca  == 1; InvRivDatsave.d42Ca  = NaN(maxsuccess,1); InvRivDat.d42Ca  = NaN;  smpres.smpresultsd42Ca  = NaN(nEM,maxsuccess); end % d42Ca
          if go.d44Ca  == 1; InvRivDatsave.d44Ca  = NaN(maxsuccess,1); InvRivDat.d44Ca  = NaN;  smpres.smpresultsd44Ca  = NaN(nEM,maxsuccess); end % d44Ca
          if go.d56Fe  == 1; InvRivDatsave.d56Fe  = NaN(maxsuccess,1); InvRivDat.d56Fe  = NaN;  smpres.smpresultsd56Fe  = NaN(nEM,maxsuccess); end % d56Fe
          if go.Sr8786 == 1; InvRivDatsave.Sr8786 = NaN(maxsuccess,1); InvRivDat.Sr8786 = NaN;  smpres.smpresultsSr8786 = NaN(nEM,maxsuccess); end % 87Sr/86Sr
          if go.d98Mo  == 1; InvRivDatsave.d98Mo  = NaN(maxsuccess,1); InvRivDat.d98Mo  = NaN;  smpres.smpresultsd98Mo  = NaN(nEM,maxsuccess); end % d98Mo
          if go.Os8788 == 1; InvRivDatsave.Os8788 = NaN(maxsuccess,1); InvRivDat.Os8788 = NaN;  smpres.smpresultsOs8788 = NaN(nEM,maxsuccess); end % 187Os/188Os
          if go.Fmod   == 1; InvRivDatsave.Fmod   = NaN(maxsuccess,1); InvRivDat.Fmod   = NaN;  smpres.smpresultsFmod   = NaN(nEM,maxsuccess); end % Fmod

end       % end of function