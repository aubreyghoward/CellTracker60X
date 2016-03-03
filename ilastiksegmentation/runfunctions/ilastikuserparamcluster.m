function ilastikuserparamcluster

global userparam
%
% file paths
% segfiledir: directory path of ilastik 2d segmentation probability density maps
% rawfiledir: directory path of the nuclear channel raw images (the same
% that was fed to FISH)
% mrnafilepath : directory path of spatzcells output.
% outfiledir:  output file path, wherever you want to save sampleoutput.mat
% imageoutputdir: if you want to save segmented image files, declare this
% path too. Else, comment it out.


userparam.segfiledir = '/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/160212FISH_images/spatzcellsimages_FISHnew/ilastiksegmentationmasks';
userparam.rawfiledir = '/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/160212FISH_images/spatzcellsimages_FISHnew/images00'; 
userparam.mrnafilepath = '/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/160212FISH_images/spatzcellsimages_FISHnew/spatzcellspotcount';
userparam.outfiledir =  '/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/160212FISH_images/spatzcellsimages_FISHnew/peaksoutputdirec';
userparam.imageoutputdir = '/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/160212FISH_images/spatzcellsimages_FISHnew/sementedimages';
userparam.fluorpdir = 0; % if there is no fluorescent protein quantification, set this to 0; else, specify the path corresponding to images with proteinchannel.
%%
userparam.nsamples = 23;
userparam.nucchannel = 0;
userparam.pchannel = 2; % fluorescence channel for protein quantification
userparam.channels = [1 2 3]; %fluorescent channel for which mRNA's are counted and need to be assigned
userparam.negativecontrol = 1;
userparam.imviews = 0;
userparam.imsave = 1; % to save all the segmented images for all channels separately (set this to 0 if you have a lot of samples.)

%%
% parameters for function primary filter
%  % number of zslices in the image
userparam.logfilter = 10;
userparam.bthreshfilter = 0.25;
userparam.diskfilter = 3;
userparam.area1filter = 100;

%%
% parameters for function secondary filter
userparam.minstartobj = 4;
userparam.minsolidity = [0.9, 0.8];
userparam.area2filter = 300;

%%
% tracing objects
% zmatch: zslices across which an object needs to be tracked. 
% if you think that one object can be present in a maximum of 5 zslices, set that as the limit.
% matchdistance: minimum distance in pixels between the centroid of two
% objectsin two different zslices to be considered as the same object.

userparam.zmatch = 8; 
userparam.matchdistance = 15;

%%
% overlap filter
userparam.overlapthresh = 70; 
%%
% distance between centroid of cell and mrna. 
userparam.cmcenter = 70;

end