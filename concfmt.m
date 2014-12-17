function s=concfmt(c)
if c>=1
  units='M';
elseif c>=1e-3
  units='mM'; c=c*1000;
elseif c>=1e-6
  units='uM'; c=c*1e6;
elseif c>=1e-9
  units='nM'; c=c*1e9;
elseif c>=1e-12
  units='pM'; c=c*1e12;
elseif c>=1e-15
  units='fM'; c=c*1e15;
else
  s=sprintf('%6.2g M',c);
  return;
end
s=sprintf('%5.1f %-2.2s',c,units);

    