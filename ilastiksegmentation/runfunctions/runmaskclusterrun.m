function runmaskclusterrun(direc, paramfile)

global userparam;

eval(paramfile);

sn = userparam.nsamples;


for i = 1:sn
    
       
    scriptfile = [direc filesep 'jobscripts' filesep 'jobcl' int2str(i)  '.slurm'];
    fid = fopen(scriptfile,'w');
    fprintf(fid, '#!/bin/bash \n');
    fprintf(fid,['#SBATCH --job-name=job' int2str(i) '\n']); 
    fprintf(fid, '#SBATCH --ntasks=1 \n');
    fprintf(fid,'#SBATCH --partition=commons\n');
    %fprintf(fid,'#SBATCH --nodes=2 \n');
    
    if( i>1) % positive control
       fprintf(fid,'#SBATCH --time=01:00:00 \n');
    %elseif(i == 3)
       %fprintf(fid,'#SBATCH --time=01:00:00 \n');
    else
        fprintf(fid,'#SBATCH --time=00:05:00 \n');
    end
    
    fprintf(fid,'#SBATCH --mail-type ALL \n');
    fprintf(fid,'#SBATCH --mail-user sc65@rice.edu \n');  
    fprintf(fid, '#SBATCH --account=commons \n'); 
    fprintf(fid, 'echo "My job ran on:" \n');
    fprintf(fid, 'cat $SLURM_NODELIST \n');
    
    %fprintf(fid, ['matlab -nodisplay -r "addpath(genpath(''/home/sc65/CellTracker''))']);
    fprintf(fid, ['matlab -nodisplay -r "runmaskdir_diffzcluster(''' paramfile ''', ' int2str(i) '); quit"']);
    fclose(fid);
    system(['sbatch ' scriptfile]);
    
  
end