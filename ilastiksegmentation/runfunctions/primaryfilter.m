   
function pmasks = primaryfilter(pnuc, logfilter, bthreshfilter, diskfilter, areafilter)


    pmasks = false(size(pnuc));
    logim2 = zeros(size(pnuc));
    s = logfilter;
    h  = fspecial('log',s);
   
   
    
    for z= 1:size(pnuc,3)
        im = pnuc(:,:,z);
        logim2(:,:,z) = imfilter(im, h);
        
    end
    
   
    t = graythresh(logim2);
    se = strel('disk',diskfilter);
    
    for z = 1:size(pnuc,3)
        
        logimsl = adapthisteq(logim2(:,:,z));
        im1 = logimsl>bthreshfilter*t;
        im2 = imfill(im1, 'holes');
        im3 = imerode(im2, se);
        tmp = bwareaopen(im3,areafilter);
        
        pmasks(:,:,z) = tmp;
    end
    
end