% Build a label for a sequence
function s=getlabel(seq,labels)
  maxliteral=6;
  keys=labels.keys;
  % Build a table of label locations (lnum,startpos,length)
  locs=[];
  cover=[];
  for i=1:length(keys)
    p=strfind(seq,keys{i});
    for j=1:length(p)
      locs(end+1,:)=[i,p(j),length(keys{i})];
      cover(end+1,:)=zeros(1,length(seq));
      cover(end,p(j)+(1:length(keys{i}))-1)=1;
    end
  end
  [~,ord]=sort(locs(:,3),'descend');
  locs=locs(ord,:);
  cover=cover(ord,:);
  for i=1:size(locs,1)
    if sum(cover(i,:))>maxliteral
      % Keep this one, not needed for any others
      for j=i+1:size(locs,1)
        cover(j,:)=cover(j,:)&~cover(i,:);
      end
    else
      cover(i,:)=0;   % Not used
    end
  end
  sel=sum(cover')>maxliteral;
  locs=locs(sel,:); cover=cover(sel,:);
  i=1;
  s='';
  while i<=length(seq)
    sel=find(locs(:,2)<=i & locs(:,2)+locs(:,3)>i,1);   % Valid ones
    if isempty(sel)
      s(end+1)=seq(i);
      i=i+1;
      continue;
    end
    l=locs(sel,:);
    %fprintf('At i=%d, sel=%d, l=[%d,%d,%d]\n',i,sel,l);
    if l(2)<i
      s=[s,sprintf('%d',i-l(2))];  % Overlap
    end
    s=[s,'<',labels(keys{l(1)}),'>'];
    i=l(2)+l(3);
  end
end

