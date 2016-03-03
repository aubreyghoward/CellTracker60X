function nucleimrnasave(masterCC, inuc, zrange, peaks, framenum, samplenum, pos, channels, cmcenter, mrnafilepath, imageoutputpath)

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


% for i = 1:size(nucleicenterlist,1)
%     text(nucleicenterlist(i,1), nucleicenterlist(i,2), int2str(i), 'Color', 'w', 'FontSize', 12);
%     
% end

%%
% assign mrna's to respective cells
for ii = 1:numel(channels)
    
    mrnafile = strcat(mrnafilepath, filesep, sprintf('ch%dallspots.mat', channels(ii)));
    load(mrnafile);
    
    spots = spotinfomat(spotinfomat(:,1) == framenum,:);
    
    spots(any(isnan(spots),2), :) = [];
    
    spotspos = [spots(:,3:4), spots(:,2)];
    
    cellmrna = cell(1, size(nucleicenterlist,1));
    
    for i = 1: size((spots),1)
        
        x0 = spotspos(i,1);
        y0 = spotspos(i,2);
        z0 = spotspos(i,3);
        
        mydist = sqrt((nucleicenterlist(:,1) - x0).^2 + (nucleicenterlist(:,2) - y0).^2 + (nucleicenterlist(:,3) - z0).^2);
        [dist, celln] = min(mydist);
        
        if (dist < cmcenter)
            if (~ isempty(cellmrna{celln}))
                newrow = size(cellmrna{celln},1)+1;
                cellmrna{celln}(newrow,:)= spotspos(i,:);
            else
                cellmrna{celln}(1,:) = spotspos(i,:);
            end
            
        end
    end
    
    
    
    
    figure ('visible', 'off'); 
    imshow(overlayntr);
    hold on;
    
    for j = 1:size(nucleicenterlist,1)
        text(nucleicenterlist(j,1), nucleicenterlist(j,2), int2str(j), 'Color', 'w', 'FontSize', 12);
        
    end
    
    
    objno = [1:numel(cellmrna)];
    for jj = 1:numel(objno)
        hold on;
        if(~isempty(cellmrna{objno(jj)}))
            plot(cellmrna{objno(jj)}(:,1), cellmrna{objno(jj)}(:,2), 'g.');
        end
    end
    
    
    filename = sprintf('fish%02d_f%02d.tif', samplenum, pos);
    fullfilepath = strcat(imageoutputpath, filesep, sprintf('channel%01d', channels(ii)), filesep, filename);
    saveas(gcf, fullfilepath);
    
    close all;
end
