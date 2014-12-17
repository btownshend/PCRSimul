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
  
  % Solve ODE
  options=odeset('Stats','on');
  [t,y]=ode45(@(t,y) rxequation(t,reshape(y,size(pairconcs)),kd,args.ka),[0,args.time],initconds,options);
  keyboard
    
  
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
  dsconc=0;
  ssconc=0;

  % Check each 3' end
  for i=1:length(seqs)
    % Compute number of double stranded bonds
    unpaired=sum(c.pairfrac(endpos+(1:length(seqs{i})),end));
    paired=sum(1-c.pairfrac(endpos+(1:length(seqs{i})),end));
    ssconc=ssconc+abs(concentrations(i))*unpaired;
    dsconc=dsconc+abs(concentrations(i))*paired;

    if concentrations(i)<0
      % Blocked 3' end
      newseqs{end+1}=seqs{i};
      newconc(end+1)=concentrations(i);
      continue;
    end

    endpos=endpos+length(seqs{i});
    pfrac=c.pairfrac(endpos,:);
    sel=find(pfrac>0);
    unpaired=1;
    for j=1:length(sel)
      pos=sel(j);
      if pos==size(c.pairfrac,2)
        % unpaired
        newseqs{end+1}=seqs{i};
        newconc(end+1)=concentrations(i)*unpaired;
        continue;
      end
      strand=find(pos<startpos,1)-1;
      strandpos=pos-startpos(strand)+1;
      assocfrac=(1-exp(-args.time*args.ka*max(concentrations([i,strand]))));
      unpaired=unpaired-assocfrac*pfrac(pos);
      newseqs{end+1}=[seqs{i},rc(seqs{strand}(1:strandpos-1))];
      newconc(end+1)=concentrations(i)*pfrac(pos)*assocfrac;
      fprintf('Strand %d.%d (%c) anneals to strand %2d.%-3d (%c) with frac=%g*%g=%g -> %s\n', i, length(seqs{i}), seqs{i}(end), strand, strandpos, seqs{strand}(strandpos), pfrac(pos),assocfrac, pfrac(pos)*assocfrac, concfmt(newconc(end)));
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
          elseif pos+length(seqs{i})-k>0 && c.pairfrac(endpos-length(seqs{i})+k,pos+length(seqs{i})-k)>pfrac(pos)/2
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
  fprintf('Unpaired nucleotides: %s, paired: %s\n', concfmt(ssconc), concfmt(dsconc));
  
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
end

keyboard

function dcv=rxequation(t,c,kd,ka)
fprintf('t=%f\n',t);
u=c(:,end);  % Unpaired concentrations
c=c(:,1:end-1);
% Pairs  dCij/dt = ka*Ui*Uj-kdCij
dc=-kd.*c;
for i=1:size(c,1)
  for j=1:size(c,2)
    dc(i,j)=dc(i,j)+ka*u(i)*u(j);
  end
end
% Unpaired to balance
for i=1:length(u)
  du(i)=-sum(dc(i,:));
end
dc(:,end+1)=du;
dcv=dc(:);

