function FDK = potato_FDK(proj,param,dim)

proj_filtered = filtering(proj,param);
FDK = CTbackprojection(proj_filtered, param);
