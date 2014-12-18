% Find complexes and optionally concentrations from mix of oligos
function c=complexes(seqs,varargin)
defaults=struct('maxsize',2,'temp',55,'cutoff',.001,'concentrations',[],'verbose',false,'ordered',true,'sodium',1.0,'mg',0.0);
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
if args.ordered
  orderedopt='-ordered';
else
  orderedopt='';
end
cmd=sprintf('%s/complexes -sodium %f -mg %f -T %f -material dna -pairs -cutoff %f %s %s',nuprefix,args.sodium, args.mg,args.temp,args.cutoff, orderedopt, tmpfile);
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

c=struct('seqs',{seqs},'maxsize',args.maxsize,'sodium',args.sodium,'mg',args.mg,'temp',args.temp);

% Load .cx file and store result in c.strands, c.dG
% Contains the composition and free energy of each complex. The first column is an integer complex identifier, 
% and the next |Ψ0| + 1 columns are L1 L2 . . . L|Ψ0| ∆G.
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
  c.complex(id)=struct('id',id,'strands',vals(2:end-1),'dG',vals(end));
end
fclose(fd);

if args.ordered
  % Load .ocx and .ocx-key files and store result in
  % prefix.ocx:
  % Generated if -ordered is selected. Contains the composition and free energy of each ordered complex. 
  % The first and second columns are integer complex and ordered complex identifiers, respectively, 
  % and the remaining |Ψ0| + 1 columns are L1 L2 . . . L|Ψ0| ∆G for the ordered complex.
  fd=fopen([tmpfile,'.ocx'],'r');
  c.ocomplex=[];
  while true
    line=fgetl(fd);
    if ~ischar(line)
      break;
    end
    if line(1)=='%'
      continue;
    end
    vals=sscanf(line,['%d ',fmt]);
    if length(vals)~=length(seqs)+3
      error('Unable to parse line: %s\n', line);
    end
    id=vals(1);
    oid=vals(2);
    c.ocomplex=[c.ocomplex,struct('id',id,'oid',oid,'strands',vals(3:end-1),'dG',vals(end),'perm',[])];
  end
  fclose(fd);

  % prefix.ocx-key:
  % Generated if -ordered is selected. Contains the distinct circular permutation of strands for each 
  % ordered complex. The first and second columns are integer complex and ordered complex identifiers, 
  % respectively, and the remaining L columns are integers from the range 1 to |Ψ0|. Note that the value 
  % of L may be different for each complex.
  fd=fopen([tmpfile,'.ocx-key'],'r');
  while true
    line=fgetl(fd);
    if ~ischar(line)
      break;
    end
    if line(1)=='%'
      continue;
    end
    vals=sscanf(line,'%f '); 
    if length(vals)<2
      error('Unable to parse line: %s\n', line);
    end
    id=vals(1);
    oid=vals(2);
    sel=[c.ocomplex.id]==id & [c.ocomplex.oid]==oid;
    if sel~=1
      error('Unable to locate ordered complex %d,%d\n', id, oid);
    end
    c.ocomplex(sel).perm=vals(3:end);
  end
  fclose(fd);
end


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
  cmd=sprintf('%s/concentrations %s -sort 0 -cutoff %f -pairs %s',nuprefix, orderedopt, args.cutoff, tmpfile);
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
  while true
    line=fgetl(fd);
    if ~ischar(line)
      break;
    end
    if line(1)=='%'
      continue;
    end
    vals=sscanf(line,'%f ');
    if length(vals)<4
      error('Unable to parse line: %s\n', line);
    end
    id=vals(1);
    if args.ordered
      oid=vals(2);
      sel=[c.ocomplex.id]==id & [c.ocomplex.oid]==oid;
      if sum(sel)~=1
        error('Unable to locate ordered complex %d.%d in c.ocomplex\n',id,oid);
      end
      c.ocomplex(sel).eqconc=vals(end);
    else
      sel=[c.complex.id]==id;
      if sum(sel)~=1
        error('Unable to locate complex %d in c.complex\n',id);
      end
      c.complex(sel).eqconc=vals(end);
    end
  end
  fclose(fd);
  if args.ordered
    % Update c.complex(:).eqconc using c.ocomplex(:).eqconc
    for i=1:length(c.complex)
      c.complex(i).eqconc=sum(c.ocomplex([c.ocomplex.id]==c.complex(i).id).eqconc);
    end
  end
  
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
  fclose(fd);
end
  
eval(['!rm ',tmpfile,'.*']);
