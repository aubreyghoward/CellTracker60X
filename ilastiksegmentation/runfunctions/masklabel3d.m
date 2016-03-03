%%
function masklabel = masklabel3d(CC,nucleilist, zrange, nzslices)

masklabel = zeros(1024, 1024, nzslices);

%%
% making a labelled 3d mask
for i = 1:size(nucleilist,1)
    object = nucleilist(i,:);
    zmatch = find(~isnan(object));
    zmatchn = zrange(zmatch);
    
    
    objectsmatch = object(zmatch);
    
    m = 1;
    for z = zmatchn
        pixelz = CC{z}.PixelIdxList{objectsmatch(m)};
        
        [pixelz_row, pixelz_col] = ind2sub([1024 1024], pixelz);
        
        for count = 1:size(pixelz_row,1)
            masklabel(pixelz_row(count), pixelz_col(count), z) = i;
        end
        m = m+1;
        clear pixelz pixelz_row pixelz_col;
    end
    
    
end




