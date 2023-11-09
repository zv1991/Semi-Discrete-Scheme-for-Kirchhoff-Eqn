data1       = [osc,beta1,beta0,min_res,max_res,r_squared];

%%% Generate CSV file for SLR parameters %%%
filecsv1    = sprintf('Table_Test%d_osc=%d.csv',problem,osc);
csvfile1    = fullfile(dirname,filecsv1);
cHeader     = {'osc','slope','yintercept','min_res','max_res','r_squared'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader  = strjoin(cHeader, ','); %cHeader in text with commas
clear ans cHeader commaHeader

%write header to file
fid = fopen(csvfile1,'w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
clear ans fid textHeader

%write data to end of file
dlmwrite(csvfile1,data1,'-append','precision','%.15f');
clear filecsv1 csvfile1
