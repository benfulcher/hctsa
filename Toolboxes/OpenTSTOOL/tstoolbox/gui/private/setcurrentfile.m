function setcurrentfile(filename, lboxhandle, currfilehandle)

set(lboxhandle, 'UserData', filename );

[path,name,ext,ver] = fileparts(filename);

set(currfilehandle, 'String', name);
