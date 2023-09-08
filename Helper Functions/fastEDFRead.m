function [Header,varargout] = fastEDFRead(varargin)

% Reads EDF files over a specified range
%
% Example: [Header, D] = fastEDFRead('File', NS5file, 'Range', SampleRange);
%
% SampleRange = [startSample, endSample];
% Samples are with respect to sample rate of file
%
% Version Date: 20171009
% Author: Tyler Davis



% Parsing input
p = inputParser;
defaultRange = [];
defaultFile = '';
addParameter(p,'Range',defaultRange,@(x)length(x)==2);
addParameter(p,'File',defaultFile,@(x)exist(x,'file'));

parse(p,varargin{:});

% Defining variables
Range = p.Results.Range;
FNameIn = p.Results.File;

if isempty(FNameIn)
    [PathIn,NameIn,ExtIn] = lastPath('\*.edf','Choose edf file...');
    FNameIn = fullfile(PathIn,[NameIn,ExtIn]);
end

% Open file
[fid,msg] = fopen(FNameIn,'r');
if fid == -1
    error(msg)
end

% HEADER
Header.ver        = str2double(char(fread(fid,8)'));
Header.patientID  = fread(fid,80,'*char')';
Header.recordID   = fread(fid,80,'*char')';
Header.startdate  = fread(fid,8,'*char')';% (dd.mm.yy)
Header.starttime  = fread(fid,8,'*char')';% (hh.mm.ss)
Header.bytes      = str2double(fread(fid,8,'*char')'); %header bytes
reserved          = fread(fid,44);
Header.records    = str2double(fread(fid,8,'*char')'); %number of records per channel
Header.duration   = str2double(fread(fid,8,'*char')'); %duration per record

% Number of signals (channels)
Header.ChannelCount = str2double(fread(fid,4,'*char')');
for ii = 1:Header.ChannelCount
    Header.label{ii} = regexprep(fread(fid,16,'*char')','\W','');
end
for ii = 1:Header.ChannelCount
    Header.transducer{ii} = fread(fid,80,'*char')';
end

% Physical dimension
for ii = 1:Header.ChannelCount
    Header.units{ii} = fread(fid,8,'*char')';
end

% Physical minimum
for ii = 1:Header.ChannelCount
    Header.physicalMin(ii) = str2double(fread(fid,8,'*char')');
end

% Physical maximum
for ii = 1:Header.ChannelCount
    Header.physicalMax(ii) = str2double(fread(fid,8,'*char')');
end

% Digital minimum
for ii = 1:Header.ChannelCount
    Header.digitalMin(ii) = str2double(fread(fid,8,'*char')');
end

% Digital maximum
for ii = 1:Header.ChannelCount
    Header.digitalMax(ii) = str2double(fread(fid,8,'*char')');
end
for ii = 1:Header.ChannelCount
    Header.prefilter{ii} = fread(fid,80,'*char')';
end
for ii = 1:Header.ChannelCount
    Header.samples(ii) = str2double(fread(fid,8,'*char')'); %samples per record per channel
end
for ii = 1:Header.ChannelCount
    reserved    = fread(fid,32,'*char')';
end

Header.label = regexprep(Header.label,'\W','');
Header.units = regexprep(Header.units,'\W','');

BegOfData = ftell(fid);
fseek(fid,0,'eof');
EndOfFile = ftell(fid);
fseek(fid,BegOfData,'bof');

Header.DataBytes = EndOfFile - BegOfData;
Header.ChannelSamples = Header.DataBytes/Header.ChannelCount/2;
if Header.records==-1 %if file not closed properly
    Header.records = floor(Header.ChannelSamples/Header.samples(1));
    Header.ChannelSamples = Header.records*Header.samples(1);
    Header.DataBytes = Header.ChannelSamples*Header.ChannelCount*2;
end

Header.Sf = abs((Header.physicalMax - Header.physicalMin)./(Header.digitalMax - Header.digitalMin));
Header.Dc = Header.physicalMax - Header.Sf .* Header.digitalMax;
Header.Fs = ceil(Header.samples./Header.duration);

%Accounting for the one extra channel problem; may be wrong (ALEX)
%Header.samples(Header.ChannelCount) = Header.samples(Header.ChannelCount-1)
%Header.ChannelCount = Header.ChannelCount-1

% Reading data
if nargout>1
    if all(Header.samples(1)==Header.samples)
        
        if isempty(Range)
            % Determining system memory to maximize data segments
            if ispc
                SystemMemory = regexp(evalc('feature memstats'),'\d*(?= MB)','match');
                SystemMemory = str2double(SystemMemory{2})*1e6; % Units bytes
            else
                [r,w] = unix('free | grep Mem');
                stats = str2double(regexp(w, '[0-9]*', 'match'));
                memsize = stats(1)/1e6;
                SystemMemory = (stats(3)+stats(end))/1e6;
            end
            
            
            MaxSamples = floor((0.75*SystemMemory)/Header.DataBytes*Header.ChannelSamples);
            Range = double([1,min(MaxSamples,Header.ChannelSamples)]);
            if MaxSamples<Header.ChannelSamples
                fprintf('Warning!! Data is too large to load. Only the 1st %0.0f samples will be loaded.\n',MaxSamples)
            end
        end
        
        % Seeking to beginning of data segment
        RecordSamples = Header.samples(1)*Header.ChannelCount;
        StartBytes = floor(Range(1)/Header.samples(1))*RecordSamples*2;
        EndBytes = ceil(Range(2)/Header.samples(1))*RecordSamples*2;
        TotalSamples = (EndBytes-StartBytes)/2;
        ChannelSamples = TotalSamples/Header.ChannelCount;
        Records = TotalSamples/RecordSamples;
        
        fseek(fid,BegOfData+StartBytes,'bof');
        
        sl = Header.samples(1); %samples per record
        cc = Header.ChannelCount;
        varargout{1} = zeros(ChannelSamples,cc,'int16');
        for k = 1:Records
            if ~rem(k,100)
                clc; disp(round(k/Records*100))
            end
            varargout{1}((1:sl)+sl*(k-1),1:cc) = fread(fid,[sl,cc],'*int16');
        end
        TrimStart = rem(Range(1),Header.samples(1))-1;
        varargout{1}(1:TrimStart,:) = [];
        if rem(Range(2),Header.samples(1))
            TrimEnd = Header.samples(1)-rem(Range(2),Header.samples(1))-1;
            varargout{1}(end-TrimEnd:end,:) = [];
        end
        if (diff(Range)+1)==size(varargout{1},1)
            Header.ChannelSamples = diff(Range)+1;
        else
            disp('Channel samples does not equal specified input size!')
        end
        
    else
        disp('sample lengths not equal')
        disp (Header.samples)
    end
end
fclose(fid);
