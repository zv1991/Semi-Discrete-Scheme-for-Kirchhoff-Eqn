data2       = [tau1,error];

%%% Generate CSV file for Log-log-scale graph data %%%
filecsv2    = sprintf('Log_graph_Test%d_osc=%d.csv',problem,osc);
csvfile2    = fullfile(dirname,filecsv2);
cHeader     = {'tau','error'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader  = strjoin(cHeader, ','); %cHeader in text with commas
clear ans cHeader commaHeader

%write header to file
fid = fopen(csvfile2,'w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
clear ans fid textHeader

%write data to end of file
dlmwrite(csvfile2,data2,'-append','precision','%.15f');
clear filecsv2 csvfile2
