function peaks = runmaskone(paramfile, samplenum,position_num, objnum)
%%
% segfiledir: directory path of ilastik 2d segmentation probability density maps
% rawfiledir: directory path of the nuclear channel raw images (the same
% that was fed to FISH)
% nzslices = no. of z slices
% colonyno: image no (1 colony implies I image or 1 unique position imaged)
% objno: cell no. for which you want to check no. of mRNA's assigned (for
% eg: if your raw image has around 20 different cells, then objno could be
% anywhere between 1 and 20.

global userparam
eval(paramfile);

ff = readFISHdir(userparam.rawfiledir, userparam.nsamples);



if(ff.positions(samplenum)< position_num)
    error('Error. Invalid position number. Try again!');
else
    nzslices = ff.zslices{samplenum}(position_num);
    nuc_ch = userparam.nucchannel;
    [pnuc, inuc] = readmaskfilesnew(userparam.segfiledir, userparam.rawfiledir, samplenum, position_num-1,  nzslices-1, nuc_ch);
    
    %%
    
    pmasks = primaryfilter(pnuc,userparam.logfilter, userparam.bthreshfilter, userparam.diskfilter, userparam.area1filter);
    %%
    
    [zrange, smasks] = secondaryfilter(pmasks, userparam.minstartobj, userparam.minsolidity, userparam.diskfilter, userparam.area2filter);
    %%
    %%
    if (zrange)
        %%
        zmatch = nzslices-1;
        [PILsn,PILsSourcen, CC, masterCCn, stats, nucleilist, zrange] = traceobjectszdistinct(smasks, userparam.matchdistance, zrange, zmatch);
        
        [nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, userparam.overlapthresh, userparam.imviews);
        
        %%
        masklabel{position_num} = masklabel3d(CC,nucleilist, zrange, nzslices);
        
        framenum = sum(ff.positions(1:samplenum-1)) + position_num;
        [peaks1{position_num}, spotsinfo{position_num}] = mrnapercells(nucleilist, stats, userparam.mrnafilepath, framenum, zrange, userparam.channels, userparam.cmcenter, masklabel{position_num});
        %%
        
        %%
        peaks{position_num} = peakscelltrackerformat(peaks1{position_num});
        if (~exist('objnum'))
            objnum = [5];
        end
        nucleimrnacheck(masterCC, inuc, zrange, peaks{position_num}, framenum, objnum, userparam.channels, userparam.mrnafilepath, userparam.cmcenter, masklabel{position_num});
        
    else
        peaks{position_num} = [];
        error('Error. \No cells found! Try another sample.')
    end
    
end
%%

