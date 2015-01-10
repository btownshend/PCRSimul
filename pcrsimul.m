% Perform PCR simulation with given sequences and concentrations
% Arguments:
%  seqs - cell array of sequences (oligos)
%	'*' at end indicate 3' end is blocked
%  concentrations - array of initial concentration of each oligo
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
classdef PCRSimul < handle
  properties
    args;
    labels;
    cycle;   % Structure array of cycle data -- cycle(1) is input, cycle(2) is result of 1 cycle, etc.
    seqids;
    maxseqid;
  end

  methods
    function obj=PCRSimul(seqs,concentrations,varargin)
      defaults=struct('temp',55,'maxsize',2,'verbose',false,'minconc',1e-12,'cutoff',1e-6,'mindisplayconc',1e-12,'labels',containers.Map(),'time',30,'ka',1e6,'sodium',0.050,'mg',0.002,'maxstrands',50);
      obj.args=processargs(defaults,varargin);

      % Remove any blanks, make sequence upper case
      for i=1:length(seqs)
        seqs{i}=upper(strrep(seqs{i},' ',''));
      end

      % Add common labels
      obj.labels=obj.args.labels;
      obj.labels('AATTTAATACGACTCACTATAGGG')='T7';
      obj.labels('CTTTTCCGTATATCTCGCCAG')='A';
      obj.labels('CGGAAATTTCAAAGGTGCTTC')='B';
      obj.labels('AAACAAACAAA')='W';
      obj.labels('AAAAAGAAAAATAAAAA')='X';
      obj.labels('GCTGTCACCGGA')='s31';
      obj.labels('TCCGGTCTGATGAGTCC')='s12';
      obj.labels('GGACGAAACAGC')='s23';

      % Add RC of all labels
      k=obj.labels.keys();
      for i=1:length(k)
        obj.labels(rc(k{i}))=[obj.labels(k{i}),'-RC'];
      end

      % Initialize tracking
      obj.cycle=struct('cyclenum',0,'seqs',{seqs},'concentrations',concentrations,'dsconc',nan,'c',[]);

      % Initialize sequence IDs
      obj.seqids=containers.Map();
      obj.maxseqid=0;
      for i=1:length(seqs)
        obj.getid(seqs{i});
      end
    end
    
    function run(obj,ncycles)  
      fprintf('Running simulation at T=%.0fC, Anneal time=%.0f sec, ka=%.1g /M/s\n', obj.args.temp, obj.args.time, obj.args.ka);

      for cycle=length(obj.cycle)+(0:ncycles-1)
        tic;
        fprintf('\n******* Cycle %d ********\n',cycle);
        [seqs,concentrations,c,ds]=obj.onecycle(obj.cycle(end).seqs,obj.cycle(end).concentrations);
        fprintf('Cycle took %.1f seconds\n', toc);
        obj.cycle(end).dsconc=ds;
        obj.cycle(end).c=c;
        obj.cycle(end+1)=struct('cyclenum',cycle,'seqs',{seqs},'concentrations',concentrations,'dsconc',nan,'c',[]);
        obj.printseqs();
        obj.plotseqs();
        pause(0.1);
      end
    end

    function id=getid(obj,seq) 
      if ~obj.seqids.isKey(seq)
        obj.maxseqid=obj.maxseqid+1;
        obj.seqids(seq)=obj.maxseqid;
      end
      id=obj.seqids(seq);
    end
      
    function printseqs(obj,cycle,mindisplayconc)
      if nargin<2 || isempty(cycle)
        cycle=length(obj.cycle);
      end
      if nargin<3 || isempty(mindisplayconc)
        mindisplayconc=obj.args.mindisplayconc;
      end
      cy=obj.cycle(cycle);
      sel=abs(cy.concentrations)>=mindisplayconc;
      sel=sel|obj.labels.isKey(cy.seqs);
      fprintf('Sequences present after cycle %d ', cy.cyclenum);
      if sum(~sel)>0
        fprintf('(showing %d/%d that are >=%s)\n', sum(sel), length(sel), concfmt(mindisplayconc));
      else
        fprintf('\n');
      end
      seqs=cy.seqs(sel);
      concs=cy.concentrations(sel);
      for i=1:length(seqs)
        fprintf('%3d %s %-25s %s\n',obj.getid(seqs{i}), concfmt(concs(i)),getlabel(seqs{i},obj.args.labels,1),seqs{i});
      end
    end
    
    function plotds(obj)
      setfig('pcrsimul-dstrack');
      plot([obj.cycle.cyclenum],[obj.cycle.dsconc]*1e6,'o-');
      xlabel('Cycle');
      ylabel('Conc(dsDNA) \mu M');
    end


    function [seqs,concentrations,c,dsconc]=onecycle(obj,seqs,concentrations)
      % Initialize for this cycle
      newseqs={};newconc=[];
      dsconc=0;ssconc=0;
      
      % Reduce number of strands going into complexes
      if length(concentrations)>obj.args.maxstrands
        [sconc]=sort(concentrations,'descend');
        sel=concentrations>=max(sconc(obj.args.maxstrands),obj.args.minconc);
      else
        sel=concentrations>=obj.args.minconc;
      end
      sel=sel|obj.labels.isKey(seqs);   % Keep any labelled ones also
      if sum(~sel)>0
        fprintf('Finding complexes for %d/%d strands with concentration >= %s\n', sum(sel), length(seqs), concfmt(sconc(obj.args.maxstrands)));
        newseqs={newseqs{:},seqs{~sel}};
        newconc=[newconc,concentrations(~sel)];
        seqs=seqs(sel);
        concentrations=concentrations(sel);
      end

      % Find all the complexes
      fprintf('Running complexes on %d strands...',length(seqs));
      c=complexes(cellfun(@(z) strrep(z,'*',''), seqs, 'UniformOutput',false),'temp',obj.args.temp,'maxsize',obj.args.maxsize,'cutoff',obj.args.cutoff,'verbose',obj.args.verbose,'concentrations',abs(concentrations),'sodium',obj.args.sodium,'mg',obj.args.mg);
      fprintf('done\n');
      
      fprintf('Found %d possible ordered complexes,',length(c.ocomplex));

      if false
        % Estimate upper bound on possible concentration of each complex over time window
        for i=1:length(c.ocomplex)
          if length(c.ocomplex(i).perm)==1
            % Max is the initial individual component concentration
            c.ocomplex(i).maxconc=abs(concentrations(c.ocomplex(i).perm));
          elseif length(c.ocomplex(i).perm)==2
            % Compute maximum possible concentration of the complexes in isolation with an irreversible reaction
            [t,y]=ode45(@(t,y) [-obj.args.ka*y(1)*y(2); -obj.args.ka*y(1)*y(2); obj.args.ka*y(1)*y(2)],[0,obj.args.time],[abs(concentrations(c.ocomplex(i).perm)),0],odeset('AbsTol',obj.args.minconc/10));
            c.ocomplex(i).maxconc=y(end,3);
          else
            error('Unsupported: complex with %d components\n', length(c.ocomplex(i).perm));
          end
        end
        % remove any complexes with maximum possible concentration < obj.args.minconc to speed up solvedynamics()
        c.allocomplex=c.ocomplex;
        c.ocomplex=c.ocomplex([c.ocomplex.maxconc]>=obj.args.minconc);
        fprintf('reduced to %d with maxconc>=%s\n', length(c.ocomplex),concfmt(obj.args.minconc));
      end

      % compute concentration after a given time window (instead of using equilibrium)
      c=solvedynamics(c,obj.args.time,'temp',obj.args.temp,'ka',obj.args.ka,'verbose',obj.args.verbose);
      
      % Loop over all the ordered complexes
      for i=1:length(c.ocomplex)
        oc=c.ocomplex(i);
        if oc.conc<obj.args.minconc
          continue;
        end

        % Compute base pair probabilities for this complex
        p=pairs(cellfun(@(z) strrep(z,'*',''), seqs, 'UniformOutput',false),oc.perm,'temp',obj.args.temp,'cutoff',obj.args.cutoff,'verbose',obj.args.verbose,'sodium',obj.args.sodium,'mg',obj.args.mg);
        c.ocomplex(i).pairs=p;	% Keep for reference
        
        % Compute number of double stranded bonds
        ssconc1=sum(p.pairfrac(:,end))*oc.conc;
        dsconc1=sum(sum(p.pairfrac(:,1:end-1)))*oc.conc;
        if oc.conc>=obj.args.mindisplayconc
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
          endpos=endpos+length(strrep(seq,'*',''));
          
          if seq(end)=='*'
            % Blocked 3' end
            newseqs{end+1}=seq;
            newconc(end+1)=oc.conc;
            fprintf('Strand %d.%d blocked extension\n',oc.perm(j),length(seq));
            continue;
          end

          % Handle unpaired by keeping this component
          newseqs{end+1}=seq;
          newconc(end+1)=p.pairfrac(endpos,end)*oc.conc;
          if newconc(end)>=obj.args.mindisplayconc
            fprintf('Strand %d.%d (%c) unpaired with frac=%g -> %s (#%d)\n', ...
                    oc.perm(j), length(seq), seq(end), p.pairfrac(endpos,end),concfmt(newconc(end)),obj.getid(newseqs{end}));
          end
          
          % Go through all the possible pairings
          for k=1:length(p.pairfrac(endpos,:))-1
            if p.pairfrac(endpos,k)==0
              continue;
            end
            seqk=strrep(seqs{strand(k)},'*','');
            pconc=p.pairfrac(endpos,k)*oc.conc;   % Concentration of this pair
            newseqs{end+1}=[seq,rc(seqk(1:strandpos(k)-1))];
            newconc(end+1)=pconc;
            if pconc>=obj.args.mindisplayconc
              fprintf('Strand %d.%d (%c) anneals to strand %2d.%-3d (%c) with frac=%g -> %s (#%d)\n', ...
                      oc.perm(j), length(seq), seq(end), strand(k), strandpos(k), seqk(strandpos(k)), p.pairfrac(endpos,k),concfmt(pconc),obj.getid(newseqs{end}));
            end
            
            if newconc(end)>obj.args.mindisplayconc
              off1=0;
              alignpos=length(seq)-1;
              off2=strandpos(k)-1+length(seq)-length(seqk);
              if off2<0
                off1=-off2;
                alignpos=alignpos-off2;
                off2=0;
              end
              fprintf(' %s5''-%s-3''\n',blanks(off1),seq);
              fprintf(' %s   ',blanks(off1));
              for m=1:length(seq)
                if length(seq)-m > length(seqk)-strandpos(k)
                  fprintf(' ');
                elseif k+length(seq)-m>0 && p.pairfrac(endpos-length(seq)+m,k+length(seq)-m)>p.pairfrac(endpos,k)/2
                  fprintf('|');
                else
                  fprintf(' ');
                end
              end
              fprintf('\n');
              fprintf(' %s3''-%s-5''\n\n',blanks(off2),seqk(end:-1:1));
            end
          end
        end
      end
      
      [seqs,ia,ic]=unique(newseqs,'sorted');
      concentrations=zeros(1,length(seqs));
      for i=1:length(ic)
        if newconc(i)<0
          error('Error: newconc(%d)<0 (%g)\n', obj.getid(newseqs{i}), newconc(i));
        end
        concentrations(ic(i))=concentrations(ic(i))+newconc(i);
      end
      [concentrations,ord]=sort(concentrations,'descend');
      seqs=seqs(ord);
    end

    function plotseqs(obj,maxplot)
    % Plot track of sequences
      if nargin<2
        maxplot=20;   % Maximum number to plot
      end
      trackseqs={};
      trackconc=[];
      for cy=1:length(obj.cycle)
        % Merge in tracking
        [trackseqs,ib,ic]=union(trackseqs,obj.cycle(cy).seqs,'stable');
        for i=1:length(trackseqs)
          sel=find(strcmp(obj.cycle(cy).seqs,trackseqs{i}));
          if ~isempty(sel)
            trackconc(i,cy)=obj.cycle(cy).concentrations(sel);
          end
        end
      end
      
      setfig('seqs');clf;
      peakconc=max(trackconc,[],2);   % Peak of each sequence
      minplotconc=0;
      if length(peakconc)>maxplot
        minplotconc=prctile(peakconc,(1-maxplot/length(peakconc))*100);
      end
      sel=peakconc>=minplotconc;
      semilogy([obj.cycle.cyclenum],trackconc(sel,:)','.-');
      xlabel('Cycle');
      ylabel('Concentration (M)');
      leg={};
      for i=1:length(trackseqs)
        if sel(i)
          leg{end+1}=sprintf('%d %s %s',obj.getid(trackseqs{i}),getlabel(trackseqs{i},obj.labels),concfmt(trackconc(i,end)));
        end
      end
      
      h=legend(leg,'Location','EastOutside');
      set(h,'Interpreter','none');
    end
  end
end
