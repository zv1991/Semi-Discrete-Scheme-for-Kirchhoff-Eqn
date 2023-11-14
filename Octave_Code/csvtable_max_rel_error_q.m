data5       = [n1,max_rel_error_q];

%%% Generate CSV file for program execution time %%%
filecsv5    = sprintf('max_rel_error_q_Test%d_osc=%d.csv',problem,osc);
csvfile5    = fullfile(dirname,filecsv5);
cHeader     = {'div_numb','max_rel_error_q'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader  = strjoin(cHeader, ','); %cHeader in text with commas
clear ans cHeader commaHeader

%write header to file
fid = fopen(csvfile5,'w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
clear ans fid textHeader

%write data to end of file
dlmwrite(csvfile5,data5,'-append','precision','%.15f');
clear filecsv5 csvfile5
