clc;
clear variables;
close all;

dataFolder = fullfile(toolboxdir('nnet'),'nndemos','nndatasets','DigitDataset');

imds = imageDatastore(dataFolder, ...
    'IncludeSubfolders',true, ....
    'LabelSource','foldernames');

imds = shuffle(imds);

imds = shuffle(imds);

[imdsTrain,imdsVal,imdsTest] = splitEachLabel(imds,0.7,0.15,0.15,'randomized');

imdsTrain = transform(imdsTrain,@(x) rescale(x));
imdsVal = transform(imdsVal,@(x) rescale(x));
imdsTest = transform(imdsTest,@(x) rescale(x));

imdsInputTrain = transform(imdsTrain,@upsampLowRes);
imdsInputVal= transform(imdsVal,@upsampLowRes);
imdsInputTest = transform(imdsTest,@upsampLowRes);

dsTrain = combine(imdsInputTrain,imdsTrain);
dsVal = combine(imdsInputVal,imdsVal);
dsTest = combine(imdsInputTest,imdsTest);

layers = unetLayers([28,28,1],2,'encoderDepth',2);
deepNetworkDesigner(layers);