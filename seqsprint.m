% seqsprint - pretty print a list of sequences
function seqsprint(seqs,concentrations,varargin)
  defaults=struct('verbose',false,'labels',containers.Map());
  args=processargs(defaults,varargin);
  %ntconc=0;
  for i=1:length(seqs)
    fprintf('%2d %s %-25s %s\n',i, concfmt(concentrations(i)),getlabel(seqs{i},args.labels,1),seqs{i});
    %  ntconc=ntconc+concentrations(i)*length(seqs{i});
  end
  %fprintf('Total nucleotide concentration = %s\n', concfmt(ntconc));
end

