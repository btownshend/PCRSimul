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

c=defaults;
c.seqs=seqs;
c.perm=ocperm;

% Parse .ppairs file
fd=fopen([tmpfile,'.ppairs'],'r');
t=textscan(fd,'%[^\n]');
fclose(fd);
lines=t{1};
c.npairs=[];
i=1;
while i<=length(lines)
  line=lines{i};
  i=i+1;
  if length(line)==0 || line(1)=='%'
    continue;
  end
  c.npairs=sscanf(line,'%d ');
  break;
end
c.pairfrac=zeros(c.npairs,c.npairs+1);
while i<=length(lines)
  line=lines{i};
  i=i+1;
  if length(line)==0 || line(1)=='%'
    continue;
  end
  vals=sscanf(line,'%d %d %g');
  c.pairfrac(vals(1),vals(2))=vals(3);
end
% Reflect in diagonal
c.pairfrac(:,1:end-1)=c.pairfrac(:,1:end-1)+c.pairfrac(:,1:end-1)';
eval(['!rm ',tmpfile,'.*']);
