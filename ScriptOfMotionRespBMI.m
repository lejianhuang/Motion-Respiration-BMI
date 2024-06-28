%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physio.dcm and confounds_timeseries.tsv were provided in the same folder.
dcm_DIR='physio.dcm';
MC_DIR='confounds_timeseries.tsv';
out_path='./';

% ---------------------------------------------------------------------
% Participants were scanned on a clinical 3 Tesla Siemens Magnetom Prisma
% with a 64 channel-head/neck coil. 
% Rs-fMRI images were acquired with the following parameters: 
% TR/TE = 555/22 ms;  number of volumes = 1110; multiband accelerator = 8. 
% During rs-fMRI scanning, a respiratory belt and pulse oximeter were 
% sampled at 400 Hz to measure respiration.

TR=0.555;
num_vol=1110;
slices=64;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adatped from extractCMRRPhysio.m (E. Auerbach, CMRR, 2016) to extract 
% and save RESP.log and Info.log files. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dcmInfo = dicominfo(dcm_DIR);
if (~isempty(dcmInfo) && isfield(dcmInfo,'ImageType') && strcmp(dcmInfo.ImageType,'ORIGINAL\PRIMARY\RAWDATA\PHYSIO') ...
        && isfield(dcmInfo,'Private_7fe1_10xx_Creator') && strcmp(deblank(char(dcmInfo.Private_7fe1_10xx_Creator)),'SIEMENS CSA NON-IMAGE'))
    np = size(dcmInfo.Private_7fe1_1010,1);
    rows = dcmInfo.AcquisitionNumber;
    columns = np/rows;
    numFiles = columns/1024;
    if (rem(np,rows) || rem(columns,1024)), error('Invalid image size (%dx%d)!', columns, rows); end
    dcmData = reshape(dcmInfo.Private_7fe1_1010,[],numFiles)';
  
    [~,~,endian] = computer;
    needswap = ~strcmp(endian,'L');
    for idx=1:numFiles
        datalen = typecast(dcmData(idx,1:4),'uint32');
        
        if needswap, datalen = swapbytes(datalen); end
        filenamelen = typecast(dcmData(idx,5:8),'uint32');
        if needswap, filenamelen = swapbytes(filenamelen); end
        filename = char(dcmData(idx,9:9+filenamelen-1));   
        if(contains(filename, '_Info.log'))
            Info_DIR = fullfile(out_path, 'Info.log');
            Info_Data = dcmData(idx,1025:1025+datalen-1);
            fp_Info = fopen(Info_DIR,'w');
            fwrite(fp_Info, char(Info_Data));
            fclose(fp_Info);
        end
        
        if(contains(filename, '_RESP.log'))
            RESP_DIR = fullfile(out_path, 'RESP.log');
            RESP_Data = dcmData(idx,1025:1025+datalen-1);
            fp_RESP = fopen(RESP_DIR,'w');
            fwrite(fp_RESP, char(RESP_Data));
            fclose(fp_RESP);
        end      
    end
else
    error('%s is not a valid encoded physio DICOM format file!', fn);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process Info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Info_fid=fopen(Info_DIR,'r');
Info_line = fgetl(Info_fid);

line_num=1;
vol_num=1;
while ischar(Info_line)
    if(line_num>=11&&line_num<11+num_vol*slices)
        Info_split=split(Info_line);
        Info_value=cellfun(@str2num,Info_split(2:6),'un',0);
        Info_matrix(vol_num,:)=Info_value';
        vol_num=vol_num+1;
    end

    if(line_num==11+num_vol*slices)
        Info_split=split(Info_line);
        first_echo=cell2mat(cellfun(@str2num,Info_split(3),'un',0));
    end

    if(line_num==11+num_vol*slices+1)
        Info_split=split(Info_line);
        last_echo=cell2mat(cellfun(@str2num,Info_split(3),'un',0));
    end

    line_num=line_num+1;
    Info_line = fgetl(Info_fid);
end
Info_table=cell2table(Info_matrix, 'VariableName', {'Vol', 'Slice', 'Start_Echo', 'End_Echo', 'Unknown'});

Info_table_DIR=fullfile(out_path, 'Info_table.csv');
writetable(Info_table, Info_table_DIR, 'Delimiter', ' ')

Echo_distance=Info_table.Start_Echo(65)-Info_table.Start_Echo(1);

for i=1:num_vol-1
    Echo(i,1)=i;
    Echo(i,2)=Info_table.Start_Echo((i-1)*64+1);
    Echo(i,3)=Info_table.Start_Echo(i*64+1)-1;
end

Echo(num_vol,1)=num_vol;
Echo(num_vol,2)=Info_table.Start_Echo((num_vol-1)*64+1);
Echo(num_vol,3)=last_echo;

Echo_table=table(Echo(:,1), Echo(:,2), Echo(:,3), 'VariableNames', {'Vol', 'Start_Echo', 'End_Echo'});

Echo_DIR=fullfile(out_path, 'Echo_table.csv');
writetable(Echo_table, Echo_DIR, 'Delimiter', ' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process RESP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RESP_fid=fopen(RESP_DIR,'r');
RESP_line = fgetl(RESP_fid);
line_num=1;
vol_num=1;

while ischar(RESP_line)
    if(line_num>=9)
        RESP_split=split(RESP_line);
        RESP_value=cellfun(@str2num,[RESP_split(2) RESP_split(4)],'un',0);
        RESP_matrix(vol_num,:)=RESP_value';
        vol_num=vol_num+1;
    end

    line_num=line_num+1;
    RESP_line = fgetl(RESP_fid);
end
RESP_table=cell2table(RESP_matrix, 'VariableName', {'Echo_Time', 'RESP'});

RESP_table_DIR=fullfile(out_path, 'RESP_table.csv');
writetable(RESP_table, RESP_table_DIR, 'Delimiter', ' ')

tick_all_start=Echo(1,2);
tick_all_end=Echo(end,3);

RESP_all_idx=find(RESP_table.Echo_Time>=tick_all_start & RESP_table.Echo_Time<=tick_all_end);
RESP_all=RESP_table.RESP(RESP_all_idx);
RESP_all_DIR=[out_path 'RESP_all.txt'];
dlmwrite(RESP_all_DIR, RESP_all);

for vol=1:1110
    tick_start=Echo(vol,2);
    tick_end=Echo(vol,3);
    RESP_idx=find(RESP_table.Echo_Time>=tick_start & RESP_table.Echo_Time<=tick_end);
    RESP_vol=RESP_table.RESP(RESP_idx);
    RESP(vol,1)=mean(RESP_vol);
end
mean_RESP=mean(RESP,'omitnan'); % calcu. nonnan mean
RESP(isnan(RESP))=mean_RESP; % replace nan with mean


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each subject, two features were derived from the z-scored 
% respiration signal: mean time between differentiated inhalations (mTdI) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F=1/TR;
index=0:num_vol/2-1;
FrequencyResolution=F/num_vol;
Frequency=index.*FrequencyResolution;
     

% find the maximum frequency of RESP signal
RESP_fft=fft(detrend(RESP));
RESP_fft_magnitude=abs(RESP_fft(index+1));
max_Fre=Frequency(find(RESP_fft_magnitude==max(RESP_fft_magnitude)));
Tmax=1/max_Fre;
Tmax_num=ceil(Tmax/0.555);
half_Tmax_num=round(Tmax_num/2);

% differene of RESP signal to be scaled to N(0, 1)
dRESP=zscore([0;diff(RESP)]);

% Inverse dRESP
DataInv = 1.01*max(dRESP) - dRESP; 
% find the starting time of difference of each breath
[Minima,MinIdx] = findpeaks(DataInv);
% set -0.5 as threshold to select the most representive 
New_MinIdx=MinIdx(find(dRESP(MinIdx)<-0.5));

% locate the peak time of each breath cycle
% and calculate mean time between differentiated inhalations (mTdI)
iter=0;
for t=2:length(New_MinIdx)
    t2=New_MinIdx(t);
    t1=New_MinIdx(t-1);
    dRESP_t1t2=dRESP(t1:1:t2);
    if(((t2-t1)>Tmax_num-half_Tmax_num) && ((t2-t1)<Tmax_num+half_Tmax_num))               
        peak_Idx=find(dRESP_t1t2==max(dRESP_t1t2));
        if(length(peak_Idx)==1) % only accounts for those that has one peak
            iter=iter+1;
            TdI(iter,1)=(peak_Idx-1)*0.555; % idx starts from 1  
        end
    end
end

mTdI=mean(TdI);

% Calculate maximum absolute values of Pearson’s correlation coefficient between
% the differentiated respiration signal (dRESP) and dMCdim across 
% 6 motion dimensions (|r|max).  
confounds = tdfread(MC_DIR); %loading the abs motion file
MC_all = [confounds.rot_x confounds.rot_y confounds.rot_z ...
    confounds.trans_x confounds.trans_y confounds.trans_z];

dMC=[zeros(1,6);diff(MC_all)];
rot_fd=50 *(dMC(:,1:3)); % 50 mm times radian of rotation displacement
trans_fd=dMC(:,4:6);
fd_all=[rot_fd trans_fd];
dRESP_fd=[dRESP,fd_all];
R_dRESP_MC=corr(dRESP_fd);
R_RESP_MC=R_dRESP_MC(1,2:7);

R_RESP_MC_max=max(abs(R_RESP_MC));




