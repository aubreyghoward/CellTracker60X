function nucleimrnacheck(masterCC, inuc, zrange, peaks, framenum, objno, channels, mrnafilepath, cmcenter, mask3d)

zstart = zrange(1);
zend = zrange(end);

%% centers of the objects
nucleicenterlist = peaks(:,1:3);
%%
% Reference figure;

masterLntr = labelmatrix(masterCC);
rgbLntr = label2rgb(masterLntr, 'jet', 'k', 'shuffle');
overlayntr = zeros(size(rgbLntr));

zstart = zrange(1);
zend = zrange(end);
for c1 = 1:3
    overlayntr(:,:,c1) = 0.5*mat2gray(rgbLntr(:,:,c1)) +...
        mat2gray(sum(inuc(:,:,zstart:zend),3));
end
figure; imshow(overlayntr);

hold on;

for i = 1:size(nucleicenterlist,1)
    text(nucleicenterlist(i,1), nucleicenterlist(i,2), int2str(i), 'Color', 'w', 'FontSize', 12);
    
end

%%
% assign mrna's to respective cells
for ii = 1:numel(channels)
    
    mrnafile = strcat(mrnafilepath, '/', sprintf('ch%dallspots.mat', channels(ii)));
    load(mrnafile);
    
    spots = spotinfomat(spotinfomat(:,1) == framenum,:);
    
    spots(any(isnan(spots),2), :) = [];
    
    spotspos = [spots(:,3:4), spots(:,2)];
    
    cellmrna = cell(1, size(nucleicenterlist,1));
    
    for i = 1: size((spots),1)
        
        x0 = spotspos(i,1);
        y0 = spotspos(i,2);
        z0 = spotspos(i,3)-1;
        
        %%
        rpxl = floor(y0+0.5);
        cpxl = floor(x0+0.5);
        
        cellnum = mask3d(rpxl, cpxl, z0);
        if cellnum > 0
            celln = cellnum;
            
        else
            mydist = sqrt((nucleicenterlist(:,1) - x0).^2 + (nucleicenterlist(:,2) - y0).^2 + (nucleicenterlist(:,3) - z0).^2);
            [dist, cellnum] = min(mydist);
            
            if(dist<cmcenter)
                celln = cellnum;
            end
        end
        
        
        if (~ isempty(cellmrna{celln}))
            newrow = size(cellmrna{celln},1)+1;
            cellmrna{celln}(newrow,:)= spotspos(i,:);
        else
            cellmrna{celln}(1,:) = spotspos(i,:);
        end
        
    end
end




figure; imshow(overlayntr);
hold on;

for j = 1:size(nucleicenterlist,1)
    text(nucleicenterlist(j,1), nucleicenterlist(j,2), int2str(j), 'Color', 'w', 'FontSize', 12);
    
end



for jj = 1:numel(objno)
    hold on;
    if(~isempty(cellmrna{objno(jj)}))
        plot(cellmrna{objno(jj)}(:,1), cellmrna{objno(jj)}(:,2), 'g.');
    end
end

end
