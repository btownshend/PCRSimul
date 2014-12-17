% seqsprint - pretty print a list of sequences
function seqsprint(seqs,concentrations,varargin)
defaults=struct('verbose',false,'labels',containers.Map());
args=processargs(defaults,varargin);
%ntconc=0;
for i=1:length(seqs)
  if args.labels.isKey(seqs{i})
    lbl=args.labels(seqs{i});
  else
    lbl='';
  end
  fprintf('%2d %s %20s %s\n',i, concfmt(concentrations(i)),lbl,seqs{i});
  %  ntconc=ntconc+concentrations(i)*length(seqs{i});
end
%fprintf('Total nucleotide concentration = %s\n', concfmt(ntconc));

