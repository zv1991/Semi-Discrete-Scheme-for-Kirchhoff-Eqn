## Author: Zurab Vashakidze <zurab.vashakidze@gmail.com>
## Created: 2023-06-16

function save_as_csv (table, filecsv, dirname, cHeader)
  ### Generate CSV file for SLR parameters ###
  csvfile     = fullfile(dirname,filecsv);
  commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; # insert commaas
  commaHeader = commaHeader(:)';
  textHeader  = strjoin(cHeader, ','); # cHeader in text with commas
  clear ('ans', 'cHeader', 'commaHeader');

  ### write header to file ###
  fid = fopen(csvfile,'w');
  fprintf(fid,'%s\n',textHeader);
  fclose(fid);
  clear ('ans', 'fid', 'textHeader');

  ### write data to end of file ###
  dlmwrite(csvfile,table,'-append','precision','%.15e');
  clear ('filecsv', 'csvfile');
endfunction
