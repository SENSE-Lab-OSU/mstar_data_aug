function RT_IMAGE_STRUCT = MSTAR_LOAD_IMAGE(PR_IMAGE_FILE,PR_HEADER_ONLY)
%RT_IMAGE_STRUCT = MSTAR_LOAD_IMAGE(PR_IMAGE_FILE,PR_HEADER_ONLY)
%
% This file opens and parses the data in the Phoenix format SAR data file
% specified in PR_IMAGE_FILE. The data is returned in a matlab structure.
% If PR_HEADER_ONLY evaluates to TRUE, only the header information will be
% read, but the image data itself will not be read. 
%

if nargin < 2
   PR_HEADER_ONLY = 0;
end

RT_IMAGE_STRUCT = struct('NumberOfRows', [],...       % number of rows in the image
                         'NumberOfColumns', [],...    % number of columns in the image
                         'TargetAz', [],...           % azimuth angle of target, if available
                         'MeasuredDepression', [],... % measured depression angle of radar platform
                         'MeasuredRange', [],...      % measured range, if available
                         'Polarization', [],...       % radar polarization
                         'TargetType',[],...          % class of the target, if available
                         'TargetSerNum',[],...        % serial number of the target, if available
                         'ImageData', [],...
                         'Bandwidth',[],...
                         'CrossRangeRes',[],...
                         'RangeRes',[],...
                         'CrossRangePix',[],...
                         'RangePix',[]);            % complex array of pixel data
                         

LV_FID = fopen(PR_IMAGE_FILE,'rt');
if LV_FID < 0
   error('MSTAR:PHOENIX:FILEOPEN','Could not open file %s',PR_IMAGE_FILE);
end

% Find the start of the Phoenix header...

LV_FOUND_IND = 0;
while (~LV_FOUND_IND)
   if feof(LV_FID)
      fclose(LV_FID);
      error('MSTAR:PHOENIX:NOHDR','Could not find header start for %s',PR_IMAGE_FILE);
   end
   
   LV_RECORD = fgets(LV_FID);
   
   LV_TEMP=sscanf(LV_RECORD,'[PhoenixHeaderVer%f]');
   if ~isempty(LV_TEMP)
      LV_FOUND_IND = 1;
   end
end

% Search remaining components looking for the header end...

LV_PHOENIX_LENGTH = 0;
LV_NATIVE_LENGTH = 0;

LV_FOUND_IND = 0;
while (~LV_FOUND_IND)
   if feof(LV_FID)
      fclose(LV_FID);
      error('MSTAR:PHOENIX:NOHDREND','Could not find header end for %s',PR_IMAGE_FILE);
   end
   
   LV_RECORD = fgets(LV_FID);
   
   LV_TEMP = sscanf(LV_RECORD,'NumberOfColumns=%d');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.NumberOfColumns = LV_TEMP;
      continue;
   end
   
   LV_TEMP = sscanf(LV_RECORD,'NumberOfRows=%d');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.NumberOfRows = LV_TEMP;
      continue;
   end
   
   LV_TEMP = sscanf(LV_RECORD,'TargetAz=%f');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.TargetAz = LV_TEMP;
      continue;
   end
   
   
   LV_TEMP = sscanf(LV_RECORD,'CrossRangeResolution=%f');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.CrossRangeRes = LV_TEMP;
      continue;
   end
   
   LV_TEMP = sscanf(LV_RECORD,'RangeResolution=%f');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.RangeRes = LV_TEMP;
      continue;
   end
   
      LV_TEMP = sscanf(LV_RECORD,'CrossRangePixelSpacing=%f');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.CrossRangePix = LV_TEMP;
      continue;
   end
   
   LV_TEMP = sscanf(LV_RECORD,'RangePixelSpacing=%f');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.RangePix = LV_TEMP;
      continue;
   end
   
   
   LV_TEMP = sscanf(LV_RECORD,'MeasuredDepression=%f');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.MeasuredDepression = LV_TEMP;
      continue;
   end
   
   LV_TEMP = sscanf(LV_RECORD,'MeasuredRange=%f');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.MeasuredRange = LV_TEMP;
      continue;
   end
   
   LV_TEMP = sscanf(LV_RECORD,'Polarization=%s');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.Polarization = LV_TEMP;
      continue;
   end
   
   LV_TEMP = sscanf(LV_RECORD,'TargetType=%s');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.TargetType = LV_TEMP;
      continue;
   end 

   LV_TEMP = sscanf(LV_RECORD,'TargetSerNum=%s');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.TargetSerNum = LV_TEMP;
      continue;
   end 

      LV_TEMP = sscanf(LV_RECORD,'Bandwidth=%s');
   if ~isempty(LV_TEMP)
      RT_IMAGE_STRUCT.Bandwidth = LV_TEMP;
      continue;
   end 
   LV_TEMP = sscanf(LV_RECORD,'PhoenixHeaderLength=%d');
   if ~isempty(LV_TEMP)
      LV_PHOENIX_LENGTH = LV_TEMP;
      continue;
   end
   
   LV_TEMP = sscanf(LV_RECORD,'PhoenixSigSize=%d');
   if ~isempty(LV_TEMP)
      LV_PHOENIX_SIG_SIZE = LV_TEMP;
      continue;
   end

   LV_TEMP = sscanf(LV_RECORD,'native_header_length=%d');
   if ~isempty(LV_TEMP)
      LV_NATIVE_LENGTH = LV_TEMP;
      continue;
   end
   
   if strncmp(LV_RECORD,'[EndofPhoenixHeader]',20) == 1
      LV_FOUND_IND = 1;
   end
end

% Close the file and be sure we got valid image dimensions...

fclose(LV_FID);

if ~PR_HEADER_ONLY

   if ~((RT_IMAGE_STRUCT.NumberOfColumns > 0)&&(RT_IMAGE_STRUCT.NumberOfRows > 0))
      error('MSTAR:PHOENIX:NODIMS','Could not find dimensions for %s',PR_IMAGE_FILE);
   end

%  Reopen the file in binary mode, seek to the image data, and read it in...

   LV_FID = fopen(PR_IMAGE_FILE,'r','ieee-be');
   if LV_FID < 0
      error('MSTAR:PHOENIX:FILEOPEN','Could not open file %s',PR_IMAGE_FILE);
   end

   fseek(LV_FID,LV_PHOENIX_LENGTH+LV_NATIVE_LENGTH,'bof');

%    NOTE: Clutter imagery is stored in 16 bit unsigned integer values, and targets are 32 bit floating point values
%    NOTE: Phoenix data is stored by rows, but Matlab reads arrays by columns...

%    Decide whether to read as uint16 or float32...

   LV_BYTES_PER_VALUE = (LV_PHOENIX_SIG_SIZE - LV_PHOENIX_LENGTH - LV_NATIVE_LENGTH)/(2*RT_IMAGE_STRUCT.NumberOfColumns*RT_IMAGE_STRUCT.NumberOfRows);

   if (LV_BYTES_PER_VALUE == 4)
      LV_MAGNITUDE_DATA = fread(LV_FID,[RT_IMAGE_STRUCT.NumberOfColumns,RT_IMAGE_STRUCT.NumberOfRows],'float32').';
      LV_PHASE_DATA = fread(LV_FID,[RT_IMAGE_STRUCT.NumberOfColumns,RT_IMAGE_STRUCT.NumberOfRows],'float32').';
   elseif (LV_BYTES_PER_VALUE == 2)
%     NOTE: 2*pi/4096 seems to be the right multiplicative factor, because it maps phase into [0,2pi). We use this same factor for both magnitude and phase.
      LV_MAGNITUDE_DATA = fread(LV_FID,[RT_IMAGE_STRUCT.NumberOfColumns,RT_IMAGE_STRUCT.NumberOfRows],'uint16').' * 2*pi/4096;
      LV_PHASE_DATA = fread(LV_FID,[RT_IMAGE_STRUCT.NumberOfColumns,RT_IMAGE_STRUCT.NumberOfRows],'uint16').' * 2*pi/4096;
   else
      fclose(LV_FID);
      error('MSTAR:PHOENIX:BITDEPTH','Do not recognize bit depth for %s',PR_IMAGE_FILE);
   end

   fclose(LV_FID);

   RT_IMAGE_STRUCT.ImageData = LV_MAGNITUDE_DATA.*(cos(LV_PHASE_DATA) + sqrt(-1)*sin(LV_PHASE_DATA));
end