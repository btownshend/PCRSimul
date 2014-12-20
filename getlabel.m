% Build a getlabel for a sequence
function s=getlabel(seq,l,depth)
  if nargin<3
    depth=1;
  end
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
  for k=length(seq)-1:-1:6
    if l.isKey(seq(1:k))
      s=['<',l(seq(1:k)),'>',getlabel(seq(k+1:end),l,depth+1)];
      return;
    end
  end

  % getlabel suffix?
  for m=1:length(seq)-6
    if l.isKey(seq(m:end))
      s=[getlabel(seq(1:m-1),l,depth+1),'<',l(seq(m:end)),'>'];
      return;
    end
  end

  % No matches - trim by 1 on each end and recheck
  if length(seq)<6
    s='*';
    return;
  end

  s=getlabel(seq(2:end-1),l,depth+1);
  if s(1)~='*'
    s=['*',s];
  end
  if s(end)~='*'
    s=[s,'*'];
  end
end

