data4       = [n1,max_CondN2];

%%% Generate CSV file for program execution time %%%
filecsv4    = sprintf('max_CondN2_Test%d_osc=%d.csv',problem,osc);
csvfile4    = fullfile(dirname,filecsv4);
cHeader     = {'div_numb','max_CondN2'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader  = strjoin(cHeader, ','); %cHeader in text with commas
clear ans cHeader commaHeader

%write header to file
fid = fopen(csvfile4,'w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
clear ans fid textHeader

%write data to end of file
dlmwrite(csvfile4,data4,'-append','precision','%.15f');
clear filecsv4 csvfile4
