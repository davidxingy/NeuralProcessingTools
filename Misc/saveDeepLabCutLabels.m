function saveDeepLabCutLabels(deeplabcutConfigFile, simiFile, videoFileFullPath, h5File)

% read in the deeplabcut config file
config = yaml.ReadYaml(deeplabcutConfigFile);

% get marker names
markerNames = config.bodyparts;
nMarkers = length(markerNames);

% 

% next read data from simi file
[header, fileMarkerNames, fileMarkerIDs, myMarkerInds, data] = readPFile(markerNames,simiFile);

% double check that the marker names match
assert(strcmpi(join(markerNames),join(fileMarkerNames)));


% write labeled data file
fid = fopen(['CollectedData_' config.scorer '.csv'],'w');

delim = ',';

% first, the header lines
% scorer:
text = ['scorer', repmat({config.scorer},1,nMarkers*2)];
fprintf(fid, join(string(text),delim));
fprintf(fid, '\n');

% markers
namesString = [markerNames; markerNames];
namesString = namesString(:)';
text = ['bodyparts', namesString];
fprintf(fid, join(string(text),delim));
fprintf(fid, '\n');

% dimension
namesString = repmat({'x'; 'y'}, 1, nMarkers);
namesString = namesString(:)';
text = ['coords', namesString];
fprintf(fid, join(string(text),delim));
fprintf(fid, '\n');

% now for the actual data, go through each frame and save as an image, and
% write the tracked location to the file
framesWithData = find(sum(isnan(data),2)~=size(data,2));
vid = VideoReader(videoFileFullPath);
fps = vid.FrameRate;
xPixels = vid.Width;
yPixels = vid.Height;
vidName = vid.Name(1:end-4);
nDigits = ceil(log10(abs(size(data,1))));

for iFrame = 1:200
    
    %get the frame and save as image
    vid.CurrentTime = framesWithData(iFrame)/fps;
    frame = readFrame(vid);
    imageName = ['img', sprintf(['%.' num2str(nDigits) 'i'], framesWithData(iFrame)), '.png'];
    imwrite(frame,imageName);
    
    %now save the coordinates in the training data file
    markerValues = data(framesWithData(iFrame),:);
    %get it in pixels, and also python is index by 0
    markerValues = markerValues .* repmat([xPixels yPixels], 1, nMarkers) -1;
    %change into the string, and remove any nans
    text = string(cellfun(@(x) num2str(x), num2cell(markerValues), 'un', 0));
    text(strcmpi(text,'nan'))= '';
    text = [['labeled-data\\' vidName '\\' imageName] text];
    fprintf(fid, join(string(text),delim));
    fprintf(fid, '\n');
    
    %add to h5 data
    h5Data.index{iFrame} = ['labeled-data\' vidName '\' imageName];
    h5Data.values_block_0(iFrame,:) = markerValues;
end

% close files
fclose(fid)

% now write the data to the h5 file that deeplabcut uses
writeDeepLabCutH5File(h5File, h5Data, nMarkers)



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Read data from simi file
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [header,fileMarkerNames,fileMarkerIDs,myMarkerInds,data]=readPFile(myMarkerNames,filename)

header=string([]);
markerNames=string([]);
markerIDs=[];
data=[];

%open file
fID=fopen(filename);

if fID==-1
    warndlg('Unable to open file')
    return
end

%read header
header=string([]);
while true
    header(end+1)=fgetl(fID);
    
    if isempty(header{end}) || strcmp(header{end},'-1')
        break
    end
end

%get markers name and ids (for some reason split with '\t'
%doesn't work, have to use char(9))
fileMarkerIDs=split(string(fgetl(fID)));
fileMarkerIDs=fileMarkerIDs(1:2:end);
fileMarkerNames=split(string(fgetl(fID)),char(9));
fileMarkerNames=fileMarkerNames(1:2:end);


%load the data
fclose(fID);
data=importdata(filename,'\t',length(header)+2);
data=data.data;

%not all file markers have corresponding data
fileMarkerNamesWithData=fileMarkerNames(1:size(data,2)/2);

%make sure the file's markers names match my marker names, and
%get the index in the loaded data
myMarkerInds=zeros(1,length(myMarkerNames));
for iMarker=1:length(myMarkerNames)
    ind=find(myMarkerNames(iMarker)==fileMarkerNamesWithData,1);
    if isempty(ind)
        myMarkerInds(iMarker)=NaN;
    else
        myMarkerInds(iMarker)=ind;
    end
end



function writeDeepLabCutH5File(filepath, data, nMarkers)
% read in the compound data type that's stored in the HDF5 table
% fileattrib(filepath,'+w');
fid = H5F.open(filepath,'H5F_ACC_RDWR','H5P_DEFAULT');
dset_id = H5D.open(fid,'/df_with_missing/table');
compound_id = H5D.get_type(dset_id);

nMembers = H5T.get_nmembers(compound_id); 
%should be two (one as the string to the image file, the second as an array
%holding the marker pixel values)
assert(nMembers == 2, 'The compound datatype in the H5 file has more than two members!')

%make sure the first member is a string
strType = H5T.copy ('H5T_C_S1');
assert(H5T.get_class(strType) == H5T.get_member_class(compound_id,0),...
    'The compound datatype in the H5 file doesn''t have a string as the first member!');

% get how long the string is
strLength = H5T.get_member_offset(compound_id,1);

% next check to make sure that the second member is an array type
base_type_id = H5T.copy('H5T_NATIVE_DOUBLE');
arrayType = H5T.array_create(base_type_id, 12);
assert(H5T.get_class(arrayType) == H5T.get_member_class(compound_id,1),...
    'The compound datatype in the H5 file doesn''t have an array as the second member!');

% check that the number of elements in the array is the same as twice the
% number of tracked joints
assert(nMarkers == (H5T.get_size(compound_id)-strLength)/8/2, ...
    'The size of the array in the compound dataset doesn''t match the number of tracked markers!');

% ok, now write the actual data
H5D.set_extent(dset_id,length(data.index));
space_id = H5D.get_space(dset_id);
spaceMat_id = H5S.create_simple(1, 1, 1);
H5S.select_hyperslab(spaceMat_id ,'H5S_SELECT_SET',0,[],[],[]);
for iRow = 1:length(data.index)
    rowData.index = data.index{iRow};
    rowData.values_block_0 = data.values_block_0(iRow,:);
    H5S.select_hyperslab(space_id,'H5S_SELECT_SET',iRow-1,[],[],[]);
    try
        H5D.write(dset_id, compound_id, spaceMat_id, space_id, 'H5P_DEFAULT', rowData);
    catch
        H5D.close(dset_id);
        H5F.close(fid);
    end
end

% release file
H5D.close(dset_id);
H5F.close(fid);










% 

