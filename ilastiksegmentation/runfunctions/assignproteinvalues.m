function peaks = assignproteinvalues(nucleilist, nzslices, fluorpdir, samplenum, pos, pchannel, inuc, peaks, mask3d)

%%

fluorchannel = size(peaks,2) +1;

m1 = 1;
    for j = 1:nzslices
        file2read = sprintf('fish%01d_f%04d_z%04d_w%04d.tif', samplenum, pos, j-1, pchannel);
        filename = strcat(fluorpdir,filesep, file2read);
        fnuc(:,:,m1) = imread(filename);
        m1 = m1+1;
    end

%%
for i = 1:size(nucleilist,1)
    
    objectlabel = find(mask3d == i);
    
    % dapi and fluorescent channel values corresponding to the same pixels
    fnuc_object = fnuc(objectlabel);
    fnuc_all = sum(fnuc_object);
    peaks(i, fluorchannel:fluorchannel+1) = fnuc_all;
    
    inuc_object = inuc(objectlabel);
    inuc_all = sum(inuc_object);
    
    peaks(i, 5) = inuc_all;
        
end