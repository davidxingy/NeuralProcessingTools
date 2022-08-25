function replaceWordsInFilenames(filesDir,wordsToReplace,replacingWords, caseSensitive)
% replaceWordsInFilenames(filesDir,wordsToReplace,replacingWords, [caseSensitive])
% 
% Renames all files within a specified directory so that any files
% containing specific words will have that word replaced with another word.
% The words can be case sensitive or not. Before doing the renaming, the
% function will list out the new filenames and confirm with the user that
% they want to carry out the renaming. Only renames files, not folders
% 
% Inputs:
% filesDir -        Directory containing the files
% 
% wordsToReplace -  Nx1 cell or string array containing the words to
%                   replace in the filenames
% 
% replacingWords -  Nx1 cell or string array the same size as
%                   wordsToReplace that contains the words that will
%                   replace the corresponding word in wordsToReplace.
% 
% caseSensitive -   Optional. Specifies whether the words being replaced
%                   should be case sensitive or not. Default true
% 
% David Xing 3/6/2019

% manage inputs
narginchk(3,4)
if nargin==3
    caseSensitive = true;
end

% if wordsToReplace or replacingWords is a single string, put it into a
% cell still
if ~iscell(wordsToReplace)
    wordsToReplace={wordsToReplace};
end
if ~iscell(replacingWords)
    replacingWords={replacingWords};
end

% check if directory exists
if ~isfolder(filesDir)
    error('Inputted path is not a folder!')
end

% number of words being replaced and words replacing must be the same
if length(wordsToReplace)~=length(replacingWords)
    error('Number of words to replace and number of replacing words must be the same!')
end

% get list of all filenames
filenames = dir(filesDir);
filenames([filenames.isdir])=[];
filenames = string({filenames.name})';

% make new filenames with the replaced words
newFilenames = filenames;

% go through each word
for iWord = 1:length(wordsToReplace)
    if caseSensitive
        newFilenames = strrep(newFilenames,wordsToReplace{iWord},replacingWords{iWord});
    else
        newFilenames = ...
            regexprep(newFilenames,wordsToReplace{iWord},replacingWords{iWord},'ignorecase');
    end
    
end

% notify user of the changes
changedFilesInds=[];
for iFile = 1:length(filenames)
    if ~strcmp(filenames{iFile},newFilenames{iFile})
        fprintf('Changing %s \t -> %s \n',filenames{iFile},newFilenames{iFile});
        changedFilesInds(end+1) = iFile;
    end
end

% ask for confirmation to continue
userInput = input('Continue? (y/n): ','s');
while ~strcmpi(userInput,'y') && ~strcmpi(userInput,'n')
    userInput = input('Please type ''y'' or ''n'': ','s');
end

if strcmpi(userInput, 'n')
    fprintf('Renaming cancelled\n')
    return
end

filenames = filenames(changedFilesInds);
newFilenames = newFilenames(changedFilesInds);

% do the actual renaming of the files
nRenamed = 0;
for iFile=1:length(filenames)
    [success, errorMsg] = movefile(fullfile(filesDir,filenames{iFile}),...
        fullfile(filesDir,newFilenames{iFile}),'f');
    
    if ~success
        warning('Unable to rename %s! \nReason: \n \t %s', filenames{iFile}, errorMsg)
    else
        nRenamed = nRenamed+1;
    end
end
fprintf('Renamed %u files\n',nRenamed);


% 
