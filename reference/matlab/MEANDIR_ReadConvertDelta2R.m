function ConvertDelta2R = MEANDIR_ReadConvertDelta2R

         % this function reads the conversion factors for moving from delta notation to 
         % isotopic ratios from the sheet "MEANDIR_DeltaNotationToR" of the workbook "MEANDIR_UserEntries"

         % open the spreadsheet containing the conversion factors
         [~, ~, all] = xlsread('MEANDIR_UserEntries','MEANDIR_DeltaNotationToR');
         header     = all(1,:);
         for i=1:length(header) 
             if isnan(header{i})
                header{i} = 'NaN';
             end
         end

         % find the list of elements, isotopes, ratios, and standards
         isotope = all(2:end,ismember(header,'IsotopeSystem'));
         ratio   = all(2:end,ismember(header,'Values'));

         % check the entries and remove unknowns
         for i=1:length(ratio);   if ~isnumeric(ratio{i}); ratio{i}   =  NaN;  end; end
         for i=1:length(isotope); if isnan(isotope{i});    isotope{i} = 'NaN'; end; end

         % convert the ratios to numbers
         ratio = cell2mat(ratio);

         % save the conversions
         clear Delta2RConversion; 
         ConvertDelta2R         = struct; 
         ConvertDelta2R.isotope = isotope;
         ConvertDelta2R.ratio   = ratio;

end % end of function