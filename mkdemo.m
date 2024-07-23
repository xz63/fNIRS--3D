if 1<0
    % save a sample data from our analsis
oxy=x.glmres{1,2,1,2,1,2}.t;
deoxy=x.glmres{1,2,1,2,2,2}.t;
xyz=x.SubjectGroupChannel{2}.polhemus.NFRI_result.OtherC;
f=FNIR23DX;
f.setXYZ(xyz);
addpath /media/BFL/BFL-Share/MATLAB/spm8

save testFNIR23DX oxy deoxy xyz
end
load testFNIR23DX
f=FNIR23DX;
f.setXYZ(xyz);
M=f.getValue(deoxy,'');  % '' is quite meaningless but the program is OK. 
f.save('deoxy_1',M);


