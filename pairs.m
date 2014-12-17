cmd=sprintf('export NUPACKHOME=/Users/bst/Dropbox/SynBio/src/nupack3.0.1; $NUPACKHOME/bin/energy -T 37 -dangles %s -material %s %s',args.nupackdangles, args.nupackparams, tmpfile);
if i==1
  fprintf('cmd=%s\n',cmd);
end
[s,r]=system(cmd);
if (s ~= 0)
  error('Error running energy');
end
newlines=strfind(r,char(10));
try 
  e(i)=str2double(r(newlines(end-1)+1:newlines(end)-1));
catch
  error('Error running NUPACK energy:\n%s', r);
end

export NUPACKHOME=/Users/bst/Dropbox/SynBio/src/nupack3.0.1
$NUPACKHOME/bin/energy -T 37 -dangles %s -material %s %s',args.nupackdangles, args.nupackparams, tmpfile);
