function [pnuc, inuc] = readmaskfilesnew(segfiledir, rawfiledir, imnum, pos,  nzslices, nuc_ch)


% reading raw data and masks
counter  =1;

for j = 0:nzslices
    file2read = sprintf('fish%01d_f%04d_z%04d_w%04d.tif', imnum, pos, j, nuc_ch);
    filename = strcat(rawfiledir,filesep, file2read);
    
    inuc(:,:,counter) = imread(filename);
    
    mask2read = sprintf('fish%01d_f%04d_z%04d_w%04d_Probabilities.h5', imnum, pos, j, nuc_ch);
    maskname = strcat(segfiledir, filesep, mask2read);
    
    mask = h5read(maskname, '/exported_data');
    
    nmask1 = squeeze(mask(1,:,:));
    nmask2 = squeeze(mask(2,:,:));
    
    pnuc(:,:,counter) = nmask1;
    pnuc(:,:,counter) = pnuc(:,:,counter)';
    
    counter = counter+1;
end


end
