function runmaskdir_diffz(paramfile)
%%
tic;
global userparam
eval(paramfile);
mkdir (userparam.outfiledir);


ff = readFISHdir(userparam.rawfiledir, userparam.nsamples);
nuc_ch = userparam.nucchannel;

% if(userparam.imsave)
%     for i = 1:numel(userparam.channels)
%         newfolder = strcat(userparam.imageoutputdir, filesep, sprintf('channel%01d', userparam.channels(i)));
%         mkdir(newfolder);
%     end
% end


for samplenum = 1:userparam.nsamples
    
    npositions = ff.positions(samplenum);
    
    clear peaks1 peaks spotsinfo masklabel;
  
    for pos = 1:npositions
        nzslices = ff.zslices{samplenum}(pos);
        
        
        [pnuc, inuc] = readmaskfilesnew(userparam.segfiledir, userparam.rawfiledir, samplenum, pos-1,  nzslices-1, nuc_ch);
        %%
        
        pmasks = primaryfilter(pnuc,userparam.logfilter, userparam.bthreshfilter, userparam.diskfilter, userparam.area1filter);
        %%
        
        [zrange, smasks] = secondaryfilter(pmasks, userparam.minstartobj, userparam.minsolidity, userparam.diskfilter, userparam.area2filter);
        
        if (zrange)
            %%
            zmatch = nzslices-1;
            [PILsn,PILsSourcen, CC, masterCCn, stats, nucleilist, zrange] = traceobjectszdistinct(smasks, userparam.matchdistance, zrange, zmatch);
            %%
            
            if(zrange)
                [nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, userparam.overlapthresh, userparam.imviews);
                
                %%
                masklabel{pos} = masklabel3d(CC,nucleilist, zrange, nzslices);
                
                framenum = sum(ff.positions(1:samplenum-1)) + pos;
                [peaks1{pos}, spotsinfo{pos}] = mrnapercells(nucleilist, stats, userparam.mrnafilepath, framenum, zrange, userparam.channels, userparam.cmcenter, masklabel{pos});
                %%
                peaks{pos} = peakscelltrackerformat(peaks1{pos});
%                 if(userparam.imsave)
%                     nucleimrnasave(masterCC, inuc, zrange, peaks{pos}, framenum, samplenum, pos-1, userparam.channels,  userparam.cmcenter, userparam.mrnafilepath,userparam.imageoutputdir);
%                 end
                
                
                if(userparam.fluorpdir)
                    if(samplenum ~= userparam.negativecontrol)
                        peaks{pos} = assignproteinvalues(nucleilist, nzslices, userparam.fluorpdir, samplenum, pos-1, userparam.pchannel, inuc, peaks{pos}, masklabel{pos});
                    end
                end
            end
        else
            peaks{pos} = [];
            masklabel{pos} = [];
            spotsinfo{pos} = 0;
        end
        
        %%
        
        
        
    end
    
    outfile = strcat(userparam.outfiledir, filesep, sprintf('sample%02dout.mat', samplenum));
    save(outfile, 'peaks', 'spotsinfo', 'masklabel');
end
toc;