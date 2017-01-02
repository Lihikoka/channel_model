% Copyright © 2016 NYU

runningFolder = pwd;
wordToFind = 'NYU_Sim';
indStop = strfind(runningFolder,wordToFind)+length(wordToFind)-1;
masterFolder = runningFolder(1:indStop(end));
cd(masterFolder)
addpath(genpath(masterFolder));

