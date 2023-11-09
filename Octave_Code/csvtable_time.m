data3       = [n1,time];

%%% Generate CSV file for program execution time %%%
filecsv3    = sprintf('Time_Test%d_osc=%d.csv',problem,osc);
csvfile3    = fullfile(dirname,filecsv3);
cHeader     = {'div_numb','time'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader  = strjoin(cHeader, ','); %cHeader in text with commas
clear ans cHeader commaHeader

%write header to file
fid = fopen(csvfile3,'w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
clear ans fid textHeader

%write data to end of file
dlmwrite(csvfile3,data3,'-append','precision','%.15f');
clear filecsv3 csvfile3
