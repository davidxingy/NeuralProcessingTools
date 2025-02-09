function makeKilosortChanMap(baseDir,spikeGLXMetaFile)

meta = ReadMeta(spikeGLXMetaFile, baseDir);

% total channel count
nChans = str2double(meta.nSavedChans);

spikeGLXCoords = textscan(meta.snsGeomMap, '(%d:%d:%d:%d', ...
    'EndOfLine', ')', 'HeaderLines', 1 );

assert(all(cellfun(@length,spikeGLXCoords) == repmat(384,1,4)),'# of chans in metafile does not equal 384!');

chanMap = (1:nChans)';
chanMap0Ind = chanMap-1;
connected = logical(ones(nChans,1)');
connected(end) = 0;

name = 'CustomNP2ChanMap';

shankSpacing = 250;

kcoords = [double(spikeGLXCoords{1})+1; 1];
xcoords = [double(spikeGLXCoords{2})+(kcoords(1:end-1)-1)*shankSpacing; 0];
ycoords = [double(spikeGLXCoords{3}); 0];

save(fullfile(baseDir,'kilosortChanMap.mat'),'chanMap','chanMap0Ind','connected','kcoords','xcoords','ycoords');