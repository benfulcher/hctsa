% minFunc
printf('Compiling minFunc files (octave version)...\n');
## working around the lack of an -outdir option in octave's mex

function mexme(fn)
  cmd = sprintf("mkoctfile --mex --output ../compiled/%s.mex %s.c", fn, fn) ;
  [ status output ] = system(cmd) ;
  if status!=0
    error("Executing command %s\n", cmd);
  else
    delete(sprintf("%s.o", fn));
    printf("%s compiled\n", fn);
  endif
endfunction

mexme("mcholC");
mexme("lbfgsC");
mexme("lbfgsAddC");
mexme("lbfgsProdC");

printf("Done.\n")