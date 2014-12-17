% Find complexes and optionally concentrations from mix of oligos
function c=complexes(seqs,varargin)
defaults=struct('maxsize',2,'temp',55,'cutoff',.001,'concentrations',[],'verbose',false);
args=processargs(defaults,varargin);

nuprefix='export NUPACKHOME=/Users/bst/Dropbox/SynBio/src/nupack3.0.4; $NUPACKHOME/bin';
tmpfile=tempname('/tmp');
fd=fopen([tmpfile,'.in'],'w');
fprintf(fd,'%d\n',length(seqs));
for i=1:length(seqs)
  fprintf(fd,'%s\n',seqs{i});
end
fprintf(fd,'%d\n',args.maxsize);
fclose(fd);
cmd=sprintf('%s/complexes -T %f -material dna -pairs -cutoff %f %s',nuprefix,args.temp,args.cutoff, tmpfile);
if args.verbose
  fprintf('cmd=%s\n',cmd);
end
[s,r]=system(cmd);
if (s ~= 0)
  error('Error running complexes');
end
if args.verbose
  fprintf('%s',r);
end

c=struct('seqs',{seqs},'maxsize',args.maxsize,'temp',args.temp,'strands',[],'dG',[]);
fmt='%d ';
for i=1:length(seqs)
  fmt=[fmt,'%d '];
end
fmt=[fmt,'%g'];
fd=fopen([tmpfile,'.cx'],'r');
while true
  line=fgetl(fd);
  if ~ischar(line)
    break;
  end
  if line(1)=='%'
    continue;
  end
  vals=sscanf(line,fmt);
  if length(vals)~=length(seqs)+2
    error('Unable to parse line: %s\n', line);
  end
  id=vals(1);
  c.strands(id,:)=vals(2:end-1);
  c.dG(id)=vals(end);
end
fclose(fd);

if ~isempty(args.concentrations)
  % Run concentrations if given
  if length(args.concentrations)~=length(seqs)
    error('concentrations must be same length as seqs\n');
  end
  fd=fopen([tmpfile,'.con'],'w');
  for i=1:length(args.concentrations)
    fprintf(fd,'%g\n',args.concentrations(i));
  end
  fclose(fd);
  c.concentrations=args.concentrations;
  cmd=sprintf('%s/concentrations -sort 0 -cutoff %f -pairs %s',nuprefix, args.cutoff, tmpfile);
  if args.verbose
    fprintf('cmd=%s\n',cmd);
  end
  [s,r]=system(cmd);
  if (s ~= 0)
    error('Error running concentrations');
  end  
  if args.verbose
    fprintf('%s',r);
  end
  
  % Parse .eq file 
  fd=fopen([tmpfile,'.eq'],'r');
  fmt=[fmt,' %g'];
  c.complexconc=[];
  while true
    line=fgetl(fd);
    if ~ischar(line)
      break;
    end
    if line(1)=='%'
      continue;
    end
    vals=sscanf(line,fmt);
    if length(vals)~=length(seqs)+3
      error('Unable to parse line: %s\n', line);
    end
    id=vals(1);
    c.strands(id,:)=vals(2:end-2);
    c.dG(id)=vals(end-1);
    c.complexconc(id)=vals(end);
  end
  fclose(fd);

  % Parse .fpairs file
  fd=fopen([tmpfile,'.fpairs'],'r');
  c.npairs=[];
  while true
    line=fgetl(fd);
    if ~ischar(line)
      break;
    end
    if line(1)=='%'
      continue;
    end
    if isempty(c.npairs)
      c.npairs=sscanf(line,'%d ');
    else
      vals=sscanf(line,'%d %d %g');
      c.pairfrac(vals(1),vals(2))=vals(3);
    end
  end
end
  
%eval(['!rm ',tmpfile,'.*']);


