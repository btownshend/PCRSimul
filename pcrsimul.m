% Perform PCR simulation with given sequences and concentrations
% Arguments:
%  seqs - cell array of sequences (oligos)
%  concentrations - array of initial concentration of each oligo
%		    negative concentrations indicate 3' end is blocked
% Optional arguments: 
%  time: time to anneal in seconds (default:30)
%  temp: annealing temperature (default: 55)
%  sodium: sodium concentration in M (default: 0.050)
%  mg: magnesium concentration in M (default: 0.002)
%  ka: DNA association constant (/M/s) (default 1e6)
%  ncycles: number of cycles to run (default 1)
%  minconc: minimum concentration in M of sequences to keep at the end of each cycle (default: 1e-6)
%  mindisplayconc: minimum concentration to display sequences, intermediates, etc (default: 1e-12)
%  labels: container.Map() that maps sequences to user-friendly names (default: a few preset labels)
%	   reverse complements are automatically added as well
%  verbose: true to make operation more verbose
% 
% Returns:
%  pcr: struct containing setup, and ordered complexes in each cycle
function pcr=pcrsimul(seqs,concentrations,varargin)
defaults=struct('temp',55,'maxsize',2,'ncycles',1,'verbose',false,'minconc',1e-12,'cutoff',1e-6,'mindisplayconc',1e-12,'labels',containers.Map(),'time',30,'ka',1e6,'sodium',0.050,'mg',0.002);
args=processargs(defaults,varargin);

% Remove any blanks, make sequence upper case
for i=1:length(seqs)
  seqs{i}=upper(strrep(seqs{i},' ',''));
end

% Add common labels
args.labels('CTTTTCCGTATATCTCGCCAG')='A_Primer';
args.labels('CGGAAATTTCAAAGGTGCTTC')='B_Primer';
args.labels('AATTTAATACGACTCACTATAGGGAAACAAACAAAGCTGTCACCGGA')='T7_W_Primer';
args.labels('TTTTTATTTTTCTTTTTGCTGTTTCGTCC')='X_Primer';

% Add RC of all labels
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

pcr=defaults;	% Return value

for cycle=1:args.ncycles
  tic;
  fprintf('\n******* Cycle %d ********\n',cycle);
  % Find all the complexes
  c=complexes(seqs,'temp',args.temp,'maxsize',args.maxsize,'cutoff',args.cutoff,'verbose',args.verbose,'concentrations',abs(concentrations),'sodium',args.sodium,'mg',args.mg);
  fprintf('Found %d possible ordered complexes,',length(c.ocomplex));
  for i=1:length(c.ocomplex)
    if length(c.ocomplex(i).perm)==1
      % Max is the initial individual component concentration
      c.ocomplex(i).maxconc=abs(concentrations(c.ocomplex(i).perm));
    elseif length(c.ocomplex(i).perm)==2
      % Compute maximum possible concentration of the complexes in isolation with an irreversible reaction
      [t,y]=ode45(@(t,y) [-args.ka*y(1)*y(2); -args.ka*y(1)*y(2); args.ka*y(1)*y(2)],[0,args.time],[abs(concentrations(c.ocomplex(i).perm)),0],odeset('AbsTol',args.minconc/10));
      c.ocomplex(i).maxconc=y(end,3);
    else
      error('Unsupported: complex with %d components\n', length(c.ocomplex(i).perm));
    end
  end
  % remove any complexes with maximum possible concentration < args.minconc
  c.allocomplex=c.ocomplex;
  c.ocomplex=c.ocomplex([c.ocomplex.maxconc]>=args.minconc);
  fprintf('reduced to %d with maxconc>=%s\n', length(c.ocomplex),concfmt(args.minconc));
  c=solvedynamics(c,args.time,'temp',args.temp,'ka',args.ka,'verbose',args.verbose);
  
  % Initialize for this cycle
  newseqs={};newconc=[];
  dsconc=0;ssconc=0;
  
  % Loop over all the ordered complexes
  for i=1:length(c.ocomplex)
    oc=c.ocomplex(i);
    if oc.conc==0
      continue;
    end

    p=pairs(seqs,oc.perm,'temp',args.temp,'cutoff',args.cutoff,'verbose',args.verbose,'sodium',args.sodium,'mg',args.mg);
    c.ocomplex(i).pairs=p;	% Keep for reference
    
    % Compute number of double stranded bonds
    ssconc1=sum(p.pairfrac(:,end))*oc.conc;
    dsconc1=sum(sum(p.pairfrac(:,1:end-1)))*oc.conc;
    if oc.conc>=args.mindisplayconc
      fprintf('Complex %d [%s]: %s (%s@eq), unpaired nucleotides: %s, paired: %s, total: %s\n', i, sprintf('%d ',oc.perm),concfmt(oc.conc), concfmt(oc.eqconc), concfmt(ssconc1), concfmt(dsconc1), concfmt(ssconc1+dsconc1));
    end
    dsconc=dsconc+dsconc1;
    ssconc=ssconc+ssconc1;

    % Map positions to strand, strandpos
    endpos=0;
    strand=[]; strandpos=[];
    for j=1:length(oc.perm)
      seq=seqs{oc.perm(j)};
      strand((1:length(seq))+endpos)=oc.perm(j);
      strandpos((1:length(seq))+endpos)=1:length(seq);
      endpos=endpos+length(seq);
    end
      
    % Check each 3' end
    endpos=0;
    for j=1:length(oc.perm)
      seq=seqs{oc.perm(j)};
      endpos=endpos+length(seq);
      
      if concentrations(oc.perm(j))<0
        % Blocked 3' end
        newseqs{end+1}=seqs{i};
        newconc(end+1)=concentrations(oc.perm(j));
        fprintf('Strand %d.%d blocked extension\n',oc.perm(j),length(seq));
        continue;
      end

      % Handle unpaired by keeping this component
      newseqs{end+1}=seq;
      newconc(end+1)=p.pairfrac(endpos,end)*oc.conc;
      if newconc(end)>=args.mindisplayconc
        fprintf('Strand %d.%d (%c) unpaired with frac=%g -> %s\n', ...
                oc.perm(j), length(seq), seq(end), p.pairfrac(endpos,end),concfmt(newconc(end)));
      end
      
      % Go through all the possible pairings
      for k=1:length(p.pairfrac(endpos,:))-1
        if p.pairfrac(endpos,k)==0
          continue;
        end
        pconc=p.pairfrac(endpos,k)*oc.conc;   % Concentration of this pair
        if pconc>=args.mindisplayconc
          fprintf('Strand %d.%d (%c) anneals to strand %2d.%-3d (%c) with frac=%g -> %s\n', ...
                  oc.perm(j), length(seq), seq(end), strand(k), strandpos(k), seqs{strand(k)}(strandpos(k)), p.pairfrac(endpos,k),concfmt(pconc));
        end
        newseqs{end+1}=[seq,rc(seqs{strand(k)}(1:strandpos(k)-1))];
        newconc(end+1)=pconc;
      
        if newconc(end)>args.mindisplayconc
          off1=0;
          alignpos=length(seq)-1;
          off2=strandpos(k)-1+length(seq)-length(seqs{strand(k)});
          if off2<0
            off1=-off2;
            alignpos=alignpos-off2;
            off2=0;
          end
          fprintf(' %s5''-%s-3''\n',blanks(off1),seq);
          fprintf(' %s   ',blanks(off1));
          for m=1:length(seq)
            if length(seq)-m > length(seqs{strand(k)})-strandpos(k)
              fprintf(' ');
            elseif k+length(seq)-m>0 && p.pairfrac(endpos-length(seq)+m,k+length(seq)-m)>p.pairfrac(endpos,k)/2
              fprintf('|');
            else
              fprintf(' ');
            end
          end
          fprintf('\n');
          fprintf(' %s3''-%s-5''\n\n',blanks(off2),seqs{strand(k)}(end:-1:1));
        end
      end
    end
  end
  
  pcr.cycle(cycle)=c;

  dstrack(cycle)=dsconc;
  setfig('pcrsimul-track');
  plot(dstrack*1e6,'o-');
  xlabel('Cycle');
  ylabel('Conc(dsDNA) \mu M');
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
    fprintf('Keeping %d sequences with concentration >= %s\n', sum(abs(concentrations)>=args.minconc),concfmt(args.minconc));
    seqs=seqs(abs(concentrations)>=args.minconc);
    concentrations=concentrations(abs(concentrations)>=args.minconc);
  end
  [trackseqs,ib,ic]=union(trackseqs,seqs,'stable');
  for i=1:length(trackseqs)
    sel=find(strcmp(seqs,trackseqs{i}));
    if ~isempty(sel)
      trackconc(i,cycle+1)=concentrations(sel);
    end
  end
  setfig('seqs');clf;
  semilogy((1:size(trackconc,2))-1,trackconc');
  xlabel('Cycle');
  ylabel('Concentration (M)');
  leg={};
  for i=1:length(trackseqs)
    leg{i}=sprintf('%s %s',trackseqs{i},concfmt(trackconc(i,end)));
  end
  
  legend(leg,'Location','EastOutside');
  pause(0.1);
  fprintf('Cycle took %.1f seconds\n', toc);
end
pcr.trackseqs=trackseqs;
pcr.trackconc=trackconc;
pcr.dstrack=dstrack;

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