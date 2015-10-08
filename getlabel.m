% Build a label for a sequence
function s=getlabel(seq,l,depth)
  maxliteral=6;
  s='xxx';
  if nargin<3
    depth=1;
  end
  %fprintf('getlabel(%s,l,%d)\n', seq,depth);
  if isempty(seq)
    s='';
    return;
  end
  
  if l.isKey(seq)
    if depth==1
      s=['### ',l(seq),' ###'];
    else
      s=['<',l(seq),'>'];
    end
    return;
  end

  % See if we can find a prefix getlabel
  for k=length(seq)-1:-1:maxliteral
    if l.isKey(seq(1:k))
      s=['<',l(seq(1:k)),'>',getlabel(seq(k+1:end),l,depth+1)];
      return;
    end
  end

  % Check for N* prefix
  nlen=find([seq,'X']~='N',1)-1;
  if nlen>1
    %fprintf('Prefix nlen=%d\n', nlen);
    s=sprintf('<N%d>%s',nlen,getlabel(seq(nlen+1:end),l,depth+1));
    return;
  end
  
  % getlabel suffix?
  for m=1:length(seq)-maxliteral
    if l.isKey(seq(m:end))
      s=[getlabel(seq(1:m-1),l,depth+1),'<',l(seq(m:end)),'>'];
      return;
    end
  end

  % Check for N* suffix
  nlen=find([seq(end:-1:1),'X']~='N',1)-1;
  if nlen>1
    %fprintf('Suffix nlen=%d\n', nlen);
    s=sprintf('%s<N%d>',getlabel(seq(1:end-nlen),l,depth+1),nlen);
    return;
  end
  
  % No matches - trim by 1 on each end and recheck
  if length(seq)<maxliteral
    s=seq;
    return;
  end

  s=getlabel(seq(2:end-1),l,depth+1);
  s=[seq(1),s,seq(end)];
  return;
  nlead=find(s=='<',1)-1;
  if isempty(nlead)
    nlead=length(s);
  end
  if nlead<maxliteral && s(1)~='*'
    s=[seq(1),s];
  else
    s=['*',s(nlead:end)];
  end
  ntail=find(s=='>',1,'last')-1;
  if isempty(ntail)
    ntail=length(s);
  end
  if ntail<maxliteral && s(end)~='*'
    s=[s,seq(end)];
  else
    s=[s(1:end-ntail),'*'];
  end
  if strcmp(s,'**')
    s='*';
  end
end

