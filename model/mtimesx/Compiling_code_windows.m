libdir = 'mingw64';
comp = computer;
mext = mexext;
lc = length(comp);
lm = length(mext);
cbits = comp(max(1:lc-1):lc);
mbits = mext(max(1:lm-1):lm);
if( isequal(cbits,'64') || isequal(mbits,'64') )
compdir = 'win64';
largearraydims = '-largeArrayDims';
else
compdir = 'win32';
largearraydims = '';
end
lib_blas = [matlabroot '\extern\lib\' compdir '\' libdir '\libmwblas.lib'];
d = dir('mtimesx.m');
mname = [d.folder '/' d.name];
cname = [mname(1:end-2) '.c'];
mex(cname,largearraydims,lib_blas,'COMPFLAGS="$COMPFLAGS /openmp"');