function cadwell2mat(FNameIn)

% Converts TXT files from the clinical cadwell system to MAT files and
% saves as *.mat within the same directory as the source. The function
% has a single input that is the full file name (FNameIn) to the TXT file
% to be converted.
%
% Example: cadwell2mat('I:\sample.txt');
% Example: cadwell2mat;
%
% Version Date: 20200113
% Author: Tyler Davis

if ~nargin
    [path,name,ext] = lastPath('\*.txt','Choose txt data file...');
    FNameIn = fullfile(path,[name,ext]);
end

fid = fopen(FNameIn);

% chartype = 'int16=>char';
chartype = 'int8=>char';
switch chartype
    case 'int8=>char'
        charfactor = 1;
    case 'int16=>char'
        charfactor = 2; %used to convert from #characters to #bytes
end

hdr_str = fread(fid,[1,10000],chartype);

bod = regexp(hdr_str,'Date.Time','once');
if isempty(bod) %if header data does not exist, just populate with nans (this occurs if exported without comments in NeuroWorks)
    bod = 0;
    
    Header.Fs = nan;
    Header.Duration = nan;
    hdr_cell = regexp(hdr_str,'\r\n','split')';
    labels = regexp(hdr_cell{1},'\t','split');
    Header.ChanLabels = repmat({''},length(labels)-3,1); % used to be -3
    
else
    bod = regexp(hdr_str(bod:end),'%\r\n','once') + bod + 1;
    
    hdr_cell = regexp(hdr_str(1:bod),'\r\n','split')';
    patient_info = regexp(hdr_cell{8},':','split');
    Header.PatientName = strtrim(patient_info{2});
    
    labels = regexp(hdr_cell{12},'\t','split');
    Header.DataHeaderLabels = [];
    Header.ChanLabels = unique(regexprep(labels(2:end-1),'\s+','')); %removing any extra whitespace characters and saving channel labels % used to be 3:end-1
    recording_info = regexp(hdr_cell{5},'\d+/\d+/\d+ \d+:\d+:\d+','match');
    Header.StartTime = recording_info{1};
    Header.Duration = etime(datevec(recording_info{2},'mm/dd/yyyy HH:MM:SS'),datevec(recording_info{1},'mm/dd/yyyy HH:MM:SS'));
    
    recording_info = regexp(hdr_cell{6},':','split');
    Header.Units = strtrim(recording_info{2});
    
    Header.Fs = str2double(regexp(hdr_cell{9},'[\d+\.]+','match'));
    Header.ChannelCount = str2double(regexp(hdr_cell{10},'[\d+\.]+','match'));    
end
data_hdr_length = 11; %number of characters before voltage values in data segments

% if isfield(Header,'StartTime') % used to be FNameOut
%     FNameOut = fullfile(fileparts(FNameIn),[datestr(Header.StartTime,'yyyymmdd-HHMMSS'),'.mat']);
% else
%     FNameOut = regexprep(FNameIn,'.txt$','.mat');
% end
[filepath,name] = fileparts(FNameIn);
FNameOut = fullfile(filepath,strcat(name,'.mat'));

% Checking if csvmat file already exists
if exist(FNameOut,'file')
    disp('mat file already exists!')
    return
end

fseek(fid,0,'eof');
eof = ftell(fid); %file size in bytes
fseek(fid,bod*charfactor,'bof'); %bod is number of characters (not bytes)

cnt = 0; %cumulative number of samples read
chunksize = 1e8; %number of 8 (or 16) bit characters to read in at each loop iteration (not bytes) ~ default: 1e8
mobj = matfile(FNameOut);
stop_flag = false;
starttime = tic;
elapsedtime = 0;
while ~stop_flag
    clc; fprintf('%0.1f%%, %0.0fs elapsed\n',cnt/Header.Fs/Header.Duration*100,elapsedtime)
    if (eof-ftell(fid))>chunksize*charfactor %convert from character to bytes using charfactor
        data_str = fread(fid,[1,chunksize],chartype);   
    else
        data_str = fread(fid,[1,(eof-ftell(fid))./charfactor],chartype);
        stop_flag = true;
    end
    idx = regexp(data_str,'\r\n'); idx = idx(end)+1;   
    data = cell2mat(textscan(data_str(1:idx),'%f','delimiter','\t','endofline','\r\n','treatasempty',{'/',':','\s','AM','PM','SHORT','OFF','ON','AMPSAT'}));
    if stop_flag
        data_rem = rem(length(data),length(Header.ChanLabels)+data_hdr_length);
        if data_rem>0
            disp('Data has been truncated!')
        end
    else
        data_rem = 0;
    end
    data = reshape(data(1:length(data)-data_rem),length(Header.ChanLabels)+data_hdr_length,[])'; %1st 11 values are a data header (something?,hour,nan,min,nan,sec,sample)    
    
    n = size(data,1);
    mobj.data_hdr((1:n)+cnt,1) = datenum(data(:,[5,1,3,6,8,10])); %[datenum,samplenum]
    mobj.D((1:n)+cnt,1:length(Header.ChanLabels)) = single(data(:,data_hdr_length+1:end));
    if ~stop_flag
        cnt = cnt+n;
        fseek(fid,(idx-chunksize)*charfactor,'cof'); %convert from character to bytes using charfactor      
    end
    elapsedtime = toc(starttime);
end
% data_samples = diff(mobj.data_hdr([1,end],2))+1; %number of samples based on data header
% data_seconds = diff(mobj.data_hdr([1,end],1))*86400; %total number of seconds based on data header (may not be exact)
Header.ChannelSamples = size(mobj,'D'); Header.ChannelSamples(end) = [];  %number of samples per channel
mobj.Header = Header;
% if Header.Duration~=round(data_seconds)    
%     disp('Header durations and data duration do not match!')
% end
% fs = floor(data_samples/data_seconds);
% if ~any(Header.Fs==[fs-1,fs,fs+1])
%     disp('Header sampling frequency and data sampling frequency do not match!')
% end

fclose(fid);
