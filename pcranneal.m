% Perform PCR simulation with given sequences and concentrations
% negative concentrations indicate 3' end is blocked
% time: time to anneal in seconds
% ka: DNA association constant (/M/s)
function pcranneal(seqs,concentrations,varargin)
defaults=struct('temp',55,'maxsize',2,'ncycles',1,'verbose',false,'minconc',1e-12,'cutoff',1e-6,'mindisplayconc',1e-15,'labels',containers.Map(),'time',30,'ka',1e6);
args=processargs(defaults,varargin);
% Remove any blanks
for i=1:length(seqs)
  seqs{i}=strrep(seqs{i},' ','');
end

% Add generic ones
args.labels('CTTTTCCGTATATCTCGCCAG')='A_Primer';
args.labels('CGGAAATTTCAAAGGTGCTTC')='B_Primer';
args.labels('AATTTAATACGACTCACTATAGGGAAACAAACAAAGCTGTCACCGGA')='T7_W_Primer';
args.labels('TTTTTATTTTTCTTTTTGCTGTTTCGTCC')='X_Primer';
% Add RC
k=args.labels.keys();
for i=1:length(k)
  args.labels(rc(k{i}))=[args.labels(k{i}),'-RC'];
end
  
fprintf('Running simulation at T=%.0fC, Anneal time=%.0f sec, ka=%.1g /M/s\n', args.temp, args.time, args.ka);
fprintf('Initial concentrations:\n');
seqsprint(seqs,abs(concentrations),'labels',args.labels);
dstrack=[];   % Track number of ds bonds
trackseqs=seqs;
trackconc(:,1)=abs(concentrations);

for cycle=1:args.ncycles
  fprintf('\n******* Cycle %d ********\n',cycle);
  % Find all the complexes
  c=complexes(seqs,'temp',args.temp,'maxsize',args.maxsize,'cutoff',args.cutoff,'verbose',args.verbose,'concentrations',abs(concentrations));

  % Convert pair probs to concentrations of pairs
  pairconcs=c.pairfrac;
  endpos=0;
  initconds=zeros(size(pairconcs));
  for i=1:length(seqs);
    pos=endpos+(1:length(seqs{i}));
    pairconcs(pos,:)=pairconcs(pos,:)*abs(concentrations(i));
    initconds(pos,end)=abs(concentrations(i));
    endpos=endpos+length(seqs{i});
  end
  kd=args.ka./pairconcs(:,1:end-1);
  for i=1:size(kd,1)
    kd(i,:)=kd(i,:)*pairconcs(i,end);
    kd(:,i)=kd(:,i)*pairconcs(i,end);
  end
  kd(~isfinite(kd))=0;
  
  ka=(kd>0)*args.ka;   % Same fixed ka for anything that can form
  
  % Solve ODE
  endpts=cumsum(cellfun(@(z) length(z), seqs));
  trackpts=zeros(size(pairconcs));
  trackpts(endpts,end)=1;
  track=find(trackpts(:));
  fprintf('Tracking %s\n',sprintf('%d ',track));
  for i=1:length(track)
    fprintf('Seq%d 3''end @ %d: initial unpaired conc=%s, equilibrium fraction unpaired=%g, equilibrium conc=%s\n', i, endpts(i), concfmt(concentrations(i)),c.pairfrac(endpts(i),end),concfmt(pairconcs(track(i))));
  end
  % Set up sparse jacobian pattern
  % Each entry of jpattern(i,j) is non-zero iff dy(i)/dt depends on y(j)
  jpat=speye(length(initconds(:)));
  jtmp=ones(1,length(initconds(:)));
  for i=1:length(initconds(:))
    dc1=rxequation(0,reshape(jtmp,size(pairconcs)),kd,ka);
    jtmp(i)=0;
    dc2=rxequation(0,reshape(jtmp,size(pairconcs)),kd,ka);
    jtmp(i)=1;
    jpat(i,dc1~=dc2)=1;
  end
  
  options=odeset('Stats','on','RelTol',0.05,'AbsTol',1e-10,'NonNegative',1:length(initconds(:)),'NormControl','off','JPattern',jpat,'OutputFcn',@odeplot,'OutputSel',track);
  tic
  [t,y]=ode23s(@(t,y) rxequation(t,reshape(y,size(pairconcs)),kd,ka),[0,args.time],initconds(:),options);
  toc
  for i=1:length(track)
    fprintf('Seq%d 3''end @ %d: final unpaired conc=%s\n', i, endpts(i), concfmt(y(end,track(i))));
  end
  kineticpairconcs=reshape(y(end,:),size(pairconcs));
  
  % Gather seq start positions
  start=1;
  startpos=[];
  for i=1:length(seqs)
    startpos(i)=start;
    start=start+length(seqs{i});
  end
  startpos(end+1)=start;

  % Initialize for this cycle
  endpos=0;
  newseqs={};newconc=[];

  % Compute number of double stranded bonds
  ssconc=sum(kineticpairconcs(:,end));
  dsconc=sum(sum(kineticpairconcs(:,1:end-1)));
  fprintf('Unpaired nucleotides: %s, paired: %s, total: %s\n', concfmt(ssconc), concfmt(dsconc), concfmt(ssconc+dsconc));

  % Check each 3' end
  for i=1:length(seqs)
    if concentrations(i)<0
      % Blocked 3' end
      newseqs{end+1}=seqs{i};
      newconc(end+1)=concentrations(i);
      continue;
    end

    endpos=endpos+length(seqs{i});
    pc=kineticpairconcs(endpos,:);
    sel=find(pc>0);
    unpaired=1;
    for j=1:length(sel)
      pos=sel(j);
      if pos==size(c.pairfrac,2)
        % unpaired
        newseqs{end+1}=seqs{i};
        newconc(end+1)=pc(pos);
        fprintf('Strand %d.%d (%c) unannealed with frac=%g -> %s(eq),%s(kin)\n', i, length(seqs{i}), seqs{i}(end), c.pairfrac(endpos,pos),concfmt(pairconcs(endpos,pos)),concfmt(newconc(end)));
        continue;
      end
      strand=find(pos<startpos,1)-1;
      strandpos=pos-startpos(strand)+1;
      newseqs{end+1}=[seqs{i},rc(seqs{strand}(1:strandpos-1))];
      newconc(end+1)=pc(pos);
      fprintf('Strand %d.%d (%c) anneals to strand %2d.%-3d (%c) with frac=%g -> %s(eq),%s(kin)\n', i, length(seqs{i}), seqs{i}(end), strand, strandpos, seqs{strand}(strandpos), c.pairfrac(endpos,pos),concfmt(pairconcs(endpos,pos)),concfmt(newconc(end)));
      if newconc(end)>args.mindisplayconc
        off1=0;alignpos=length(seqs{i})-1;off2=strandpos-1+length(seqs{i})-length(seqs{strand});
        if off2<0
          off1=-off2;
          alignpos=alignpos-off2;
          off2=0;
        end
        fprintf(' %s%s\n',blanks(off1),seqs{i});
        fprintf(' %s',blanks(off1));
        for k=1:length(seqs{i})
          if length(seqs{i})-k > length(seqs{strand})-strandpos
            fprintf(' ');
          elseif pos+length(seqs{i})-k>0 && c.pairfrac(endpos-length(seqs{i})+k,pos+length(seqs{i})-k)>c.pairfrac(endpos,pos)/2
            fprintf('|');
          else
            fprintf(' ');
          end
        end
        fprintf('\n');
        fprintf(' %s%s\n\n',blanks(off2),seqs{strand}(end:-1:1));
      end
    end
  end
  
  dstrack(cycle)=dsconc;
  setfig('pcranneal-track');
  plot(dstrack*1e6,'o-');
  xlabel('Cycle');
  ylabel('Conc(dsDNA) µM');
  pause(0.1);
  
  % Merge duplicates
  [seqs,ia,ic]=unique(newseqs,'sorted');
  concentrations=zeros(1,length(seqs));
  for i=1:length(ic)
    assert(concentrations(ic(i))>=0);   % Make sure we con't collide with blocked ones
    concentrations(ic(i))=concentrations(ic(i))+newconc(i);
  end
  [concentrations,ord]=sort(concentrations,'descend');
  seqs=seqs(ord);
  sel=concentrations>=args.mindisplayconc;
  seqsprint(seqs(sel),abs(concentrations(sel)),'labels',args.labels);
  % Remove any below minconc
  if any(abs(concentrations)<args.minconc)
    fprintf('Removing %d sequences with concentration < %g M\n', sum(abs(concentrations)<args.minconc),args.minconc);
    seqs=seqs(abs(concentrations)>=args.minconc);
    concentrations=concentrations(abs(concentrations)>=args.minconc);
  end
  [trackseqs,ib,ic]=union(trackseqs,seqs,'stable');
  for i=1:length(trackseqs)
    sel=strcmp(seqs,trackseqs{i});
    if ~isempty(sel)
      trackconc(i,cycle+1)=concentrations(sel);
    end
  end
  setfig('seqs');clf;
  semilogy((1:size(trackconc,2))-1,trackconc'*1e12);
  xlabel('Cycle');
  ylabel('Concentration (pM)');
  legend(trackseqs);
  pause(0.1);
end

function dcv=rxequation(t,c,kd,ka)
global oldt prevt;
if isempty(prevt)
  prevt=t;
end
u=c(:,end);  % Unpaired concentrations
p=c(:,1:end-1);
% Pairs  dCij/dt = ka*Ui*Uj-kdCij
dc=-kd.*p+ka.*(u*u');
% Unpaired to balance
du=-sum(dc);
dc(:,end+1)=du;
dcv=dc(:);
if (isempty(oldt) || t-oldt>0.01 || t<oldt) && t>prevt
  relchange=abs(dc./c)*(t-prevt);
  maxrelchange=max(relchange(c>1e-12));
  %  fprintf('t=%f, maxrelchange=%.2f%%, maxabschange=%s\n',t,maxrelchange*100,concfmt((t-prevt)*max(abs(dc(c>0)))));
  oldt=t;
  set(gca,'YScale','log');
  %  keyboard
end
prevt=t;