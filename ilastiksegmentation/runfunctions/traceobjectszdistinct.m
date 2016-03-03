%%
function [PILsn, PILsSourcen, CC, masterCCn, stats, nuclein1, zrange] = traceobjectszdistinct(smasks, matchdistance, zrange, zmatch)

% tracking objects
%very simple matching

clear nuclein1 PILsn PILsSourceN overlayn
n = zmatch;

CC = {};
stats = {};
zstart = zrange(1);
zend = zrange(end);

%%

% recalculating zstart and zend

m = 1;
objectsz(1:zend) = 0;

for z = zstart:zend
    objects = bwconncomp(smasks(:,:,z));
    if(objects.NumObjects > 2);
        zvalue(1,m) = z;
        m = m+1;
    end
end

if(exist('zvalue'))
    zstart = zvalue(1,1);
    zend = zvalue(1,end);
    %%
    
    for z = zvalue
        
        CC{z} = bwconncomp(smasks(:,:,z));
        stats{z} = regionprops(CC{z}, 'Centroid', 'Area');
    end
    
    
    
    nSlices = zend-zstart+1;
    nuclein1 = nan([CC{zstart}.NumObjects nSlices]);
    dmax = matchdistance;
    
    % master PixelIdxList
    PILsn = CC{zstart}.PixelIdxList;
    
    % keep track of which slice added an object to the PILs
    PILsSourcen = ones([numel(PILsn) 1])*zstart;
    
    for i = 1:CC{zstart}.NumObjects
        nuclein1(i,1) = i;
    end
    %%
    % Irrespective of the interval check size,  the logic for matching cells of two frames
    % will be the same.
    
    if(zstart~= zend)
        for z = zvalue(2)
            clear CM1 CM2 D
            
            CM1 = cat(1,stats{zvalue(1)}.Centroid);
            CM2 = cat(1,stats{z}.Centroid);
            
            D = distmat(CM1, CM2);
            
            [row,col] = find(D < dmax);
            
            % which of the new slice are matched to the previous
            matched = sum(D < dmax, 1)';
            
            % if matched add to current track
            for j = 1:sum(matched)
                trackIdx = find(nuclein1(:,1) == row(j));
                nuclein1(trackIdx,z-zstart+1) = col(j);
                PILsn{trackIdx} = cat(1,PILsn{trackIdx}(:), CC{z}.PixelIdxList{col(j)}(:));
            end
            
            % if not matched make new track
            for j = 1:CC{z}.NumObjects
                if ~matched(j)
                    newtrack = nan([1 zend-zstart+1]);
                    newtrack(z-zstart+1) = j;
                    nuclein1 = cat(1,nuclein1, newtrack);
                    PILsn = [PILsn, CC{z}.PixelIdxList(j)];
                    PILsSourcen = cat(1, PILsSourcen, z);
                end
            end
        end
        %%
        if numel(zvalue)>2
            for z = zvalue(3:end)
                
                
                zcheck = zvalue < z;
                zcheckvalues = zvalue(zcheck);
                
                framechecksize = numel(zcheckvalues);
                
                clear CM dm row col mm c0 c00 c1;
                
                for i = 1:framechecksize
                    CM{i} = cat(1,stats{zcheckvalues(i)}.Centroid);
                end
                
                CM{framechecksize+1} = cat(1,stats{z}.Centroid);
                
                for i = 1:framechecksize
                    dm{i} = distmat(CM{i}, CM{framechecksize+1});
                end
                
                for i = 1:framechecksize
                    [row{i}, col{i}] = find(dm{i} < dmax);
                    mm{i} = sum(dm{i} < dmax, 1)';
                end
                
                
                framelimit = framechecksize;
                
                framecount = 1;
                
                
                %%
                while framecount <= framelimit
                    
                    if(framecount ==1)
                        c1 = find(mm{framecount}(:,1)==1); % matched with frame1
                        c0 = find(mm{framecount}(:,1)==0); % did not match with frame1
                    end
                    
                    %%
                    if(exist('c1', 'var'))
                        for i1 = c1'
                            clear olabel
                            olabel = find(col{framecount}(:,1) == i1);
                            zsearch = zcheckvalues(framecount)-zstart+1;
                            clear trackIdx
                            trackIdx = find(nuclein1(:,zsearch) == row{framecount}(olabel,1));
                            nuclein1(trackIdx, z-zstart+1) = i1;
                            PILsn{trackIdx} = cat(1,PILsn{trackIdx}(:), CC{z}.PixelIdxList{i1}(:));
                        end
                    end
                    clear c1;
                    %%
                    if(exist('c0', 'var') && (~isempty(c0)))
                        if(zcheckvalues(framecount) < zcheckvalues(end))
                            clear c1 c00;
                            m1 = 1;
                            m2 = 1;
                            for i1 = c0'
                                if (mm{framecount+1}(i1,1) == 1)
                                    c1(m1,1) = i1;
                                    m1 = m1+1;
                                else
                                    c00(m2,1) = i1;
                                    m2 = m2+1;
                                end
                            end
                            clear c0;
                            if(exist('c00', 'var'))
                                c0 = c00;
                            end
                        end
                    end
                    %% assign a new track for elements that did not not match with any of the previous frames.
                    if(zcheckvalues(framecount) == zcheckvalues(end))&& exist('c0', 'var')&& (~isempty(c0))
                        for i1 = c0'
                            clear trackIdx;
                            trackIdx = size(nuclein1, 1)+1;
                            nuclein1(trackIdx, 1:nSlices) = NaN;
                            nuclein1(trackIdx, z-zstart+1) = i1;
                            PILsn = [PILsn, CC{z}.PixelIdxList(i1)];
                            PILsSourcen = cat(1, PILsSourcen, z);
                        end
                        
                    end
                    
                    framecount = framecount+1;
                    
                end
                
            end
        end
    end
    
    %%
    masterCCn = CC{zstart};
    masterCCn.NumObjects = numel(PILsn);
    masterCCn.PixelIdxList = PILsn;
    zrange = [zstart:zend];
    
    
else
    masterCCn = 0;
    PILsn = 0;
    PILsSourcen = 0;
    CC = 0;
    masterCCn = 0;
    stats = 0;
    nuclein1=0;
    zrange =0;
    
    
end


end


