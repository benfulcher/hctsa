function varargout = sparseinv_mex(varargin)          % compile and run function

loc = fileparts(mfilename('fullpath'));
src = ['-outdir ',loc,' ',loc,'/sparseinv_mex.c ',loc,'/sparseinv.c'];

if ~isempty(strfind(computer,'64'))
  fprintf ('Compiling sparseinv (64-bit)\n') ;
  eval(['mex -largeArrayDims ',src])
else
  fprintf ('Compiling sparseinv (32-bit)\n') ;
  eval(['mex ',src])
end
varargout = cell(nargout,1); [varargout{:}] = sparseinv_mex (varargin{:});