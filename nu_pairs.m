% Find pair probability for given ordered complex
function c=nu_pairs(seqs,ocperm,varargin)
defaults=struct('maxsize',2,'temp',55,'cutoff',.001,'verbose',false,'sodium',1,'mg',0,'material','dna1998','dangles','some');
args=processargs(defaults,varargin);

nuprefix='export NUPACKHOME=/Users/bst/Dropbox/SynBio/src/nupack3.0.6; $NUPACKHOME/bin';
tmpfile=tempname('/tmp');
fd=fopen([tmpfile,'.in'],'w');
fprintf(fd,'%d\n',length(seqs));
for i=1:length(seqs)
  fprintf(fd,'%s\n',seqs{i});
end
fprintf(fd,'%d ',ocperm);
fprintf(fd,'\n');
fclose(fd);
cmd=sprintf('%s/pairs -T %f -material %s -dangles %s -cutoff %f -sodium %g -magnesium %g -multi %s',nuprefix,args.temp,args.material,args.dangles,args.cutoff, args.sodium, args.mg, tmpfile);
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

c=args;
c.seqs=seqs;
c.perm=ocperm;

% Parse .ppairs file
fd=fopen([tmpfile,'.ppairs'],'r');
while true
  line=fgetl(fd);
  if ~ischar(line)
    error('EOF on input file before reading npairs');
  end
  if isempty(line) || line(1)=='%'
    continue;
  end
  c.npairs=sscanf(line,'%d ');
  break;
end
c.pairfrac=zeros(c.npairs,c.npairs+1);
t=textscan(fd,'%d %d %f');
c.pairfrac(sub2ind(size(c.pairfrac),t{1},t{2}))=t{3};
fclose(fd);
% Reflect in diagonal
c.pairfrac(:,1:end-1)=c.pairfrac(:,1:end-1)+c.pairfrac(:,1:end-1)';
eval(['!rm ',tmpfile,'.*']);
