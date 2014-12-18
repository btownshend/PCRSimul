% Find pair probability for given ordered complex
function c=pairs(seqs,ocperm,varargin)
defaults=struct('maxsize',2,'temp',55,'cutoff',.001,'verbose',false,'sodium',1,'mg',0);
args=processargs(defaults,varargin);

nuprefix='export NUPACKHOME=/Users/bst/Dropbox/SynBio/src/nupack3.0.4; $NUPACKHOME/bin';
tmpfile=tempname('/tmp');
fd=fopen([tmpfile,'.in'],'w');
fprintf(fd,'%d\n',length(seqs));
for i=1:length(seqs)
  fprintf(fd,'%s\n',seqs{i});
end
fprintf(fd,'%d ',ocperm);
fprintf(fd,'\n');
fclose(fd);
cmd=sprintf('%s/pairs -T %f -material dna -cutoff %f -sodium %g -magnesium %g -multi %s',nuprefix,args.temp,args.cutoff, args.sodium, args.mg, tmpfile);
if args.verbose
  fprintf('cmd=%s\n',cmd);
end
[s,r]=system(cmd);
if (s ~= 0)
  error('Error running pairs');
end
if args.verbose
  fprintf('%s',r);
end

c=struct('seqs',{seqs},'temp',args.temp,'sodium',args.sodium,'mg',args.mg,'perm',ocperm);

% Parse .ppairs file
fd=fopen([tmpfile,'.ppairs'],'r');
c.npairs=[];
while true
  line=fgetl(fd);
  if ~ischar(line)
    break;
  end
  if length(line)==0 || line(1)=='%'
    continue;
  end
  if isempty(c.npairs)
    c.npairs=sscanf(line,'%d ');
  else
    vals=sscanf(line,'%d %d %g');
    c.pairfrac(vals(1),vals(2))=vals(3);
    if vals(2)<=c.npairs
      c.pairfrac(vals(2),vals(1))=vals(3);
    end
  end
end
fclose(fd);

eval(['!rm ',tmpfile,'.*']);
