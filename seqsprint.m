% seqsprint - pretty print a list of sequences
function seqsprint(seqs,concentrations,varargin)
  defaults=struct('verbose',false,'labels',containers.Map());
  args=processargs(defaults,varargin);
  %ntconc=0;
  for i=1:length(seqs)
    fprintf('%2d %s %20s %s\n',i, concfmt(concentrations(i)),label(seqs{i},args.labels,1),seqs{i});
    %  ntconc=ntconc+concentrations(i)*length(seqs{i});
  end
  %fprintf('Total nucleotide concentration = %s\n', concfmt(ntconc));
end

% Build a label for a sequence
function s=label(seq,l,depth)
  if isempty(seq)
    s='';
    return;
  end
  
  if l.isKey(seq)
    if depth==1
      s=l(seq);
    else
      s=['<',l(seq),'>'];
    end
    return;
  end

  % See if we can find a prefix label
  lend=0;
  for k=length(seq)-1:-1:6
    if l.isKey(seq(1:k))
      s=['<',l(seq(1:k)),'>',label(seq(k+1:end),l,depth+1)];
      return;
    end
  end

  % label suffix?
  for m=lend+1:length(seq)-6
    if l.isKey(seq(m:end))
      s=[label(seq(1:m),l,depth+1),'<',l(seq(m:end)),'>'];
      return;
    end
  end

  % No matches - trim by 1 on each end and recheck
  if length(seq)<6
    s='*';
    return;
  end

  s=label(seq(2:end-1),l,depth+1);
  if s(1)~='*'
    s=['*',s];
  end
  if s(end)~='*'
    s=[s,'*'];
  end
end

