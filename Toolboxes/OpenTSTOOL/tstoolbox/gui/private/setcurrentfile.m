function setcurrentfile(filename, lboxhandle, currfilehandle)

set(lboxhandle, 'UserData', filename );

[path,name,ext] = fileparts(filename);

set(currfilehandle, 'String', name);
