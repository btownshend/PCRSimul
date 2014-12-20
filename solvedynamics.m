% Solve dynamics of concentration changes to find concentration of each complex after time t
% Args: c - structure returned by complexes()
% Returns: d - same as c but with concentrations updated to dynamic values
function d=solvedynamics(c,t,varargin)
defaults=struct('temp',55,'ka',1e6,'mindisplayconc',1e-12,'verbose',false,'reltol',0.01,'abstol',1e-12);
args=processargs(defaults,varargin);
ncomplex=length(c.ocomplex);
fprintf('Solving dynamics for %d complexes over %f seconds\n', ncomplex, t);
% Setup up differential equation
% C(t)~ concentration of each complex at time t
% dC(k)/dt=sum(ka*C(i)*C(j)) - kd*C(m)
%    where complexes i,j can be combined to make complex k
%    or complex m can be broken down to from complex k
% Using fixed ka, find kd such that equilibrium stored in c.ocomplexconc will be reached
%   Or can use c.dG to find kd

% A2{k}(i,j) is the rate of production of k due to bimolecular reactions when scaled by c(i)*c(j);  since it is symmetric, total is A2{k}(i,j)+A2{k}(j,i) (i.e. each are 1/2)
% D2(i,j) is the rate of consumption of i and j in formation of a bimolecular complex
% A1(i,k) is the rate of production of i, when scaled by c(k)
% D1(k) is the rate of consumption of k, when scaled by c(k)
A2={};
for k=1:ncomplex
  A2{k}=spalloc(ncomplex,ncomplex,1);  % A{k}(i,j) is rate of association of i+j -> k; 
end
A1=zeros(ncomplex,ncomplex);  % A1(i,k) is rate of disocciation of k -> i
for i=1:ncomplex
  cstr{i}=sprintf('%d,',c.ocomplex(i).perm);
end

for i=1:ncomplex
  for j=1:ncomplex
    fusion=[cstr{i},cstr{j}];
    k=find(strcmp(fusion,cstr));
    if ~isempty(k)
      % Have a match
      if args.verbose
        fprintf('[%s]+[%s]->[%s]\n', sprintf('%d ',c.ocomplex(i).perm), sprintf('%d ',c.ocomplex(j).perm), sprintf('%d ',c.ocomplex(k).perm));
      end
      A2{k}(i,j)=A2{k}(i,j)+args.ka/2;
      A2{k}(j,i)=A2{k}(j,i)+args.ka/2;
      % Dissociation rate to reach equilibrium
      kd=args.ka*c.ocomplex(i).eqconc*c.ocomplex(j).eqconc/c.ocomplex(k).eqconc;
      A1(i,k)=A1(i,k)+kd;
      A1(j,k)=A1(j,k)+kd;
    end
  end
end
D1=sum(A1)'/2;
D2=zeros(ncomplex,ncomplex);
for k=1:length(A2)
  D2=D2+A2{k}*2;
end

initconds=c.concentrations;   % First complexes are strands by themselves
initconds(end+1:ncomplex)=0;
if args.verbose
  statsval='on';
else
  statsval='off';
end
if false
  jpat=speye(length(initconds));
  jtmp=ones(1,length(initconds));
  dc1=dyneqn(jtmp,A2,D);
  for i=1:length(initconds)
    jtmp(i)=0;
    dc2=dyneqn(jtmp,A2,D);
    jtmp(i)=1;
    jpat(i,dc1~=dc2)=1;
  end
end

options=odeset('Stats',statsval,'RelTol',args.reltol,'AbsTol',args.abstol,'NonNegative',1:length(initconds),'NormControl','off'); % ,'JPattern',jpat); % ,'OutputFcn',@odeplot);
[t,y]=ode15s(@(t,y) dyneqn(y,A2,D2,A1,D1),[0,t],initconds,options);
d=c;
d.time=t;
d.cconc=y;
for i=1:length(c.ocomplex)
  d.ocomplex(i).conc=d.cconc(end,i);
end

setfig('timecourse');clf;
sel=find(max(y,[],1)>args.mindisplayconc);
loglog(t,y(:,sel));
xlabel('Time (s)');
ylabel('Conc');
leg={};
for i=1:length(sel)
  leg{i}=sprintf('C%d [%s]',sel(i),sprintf('%d ',c.ocomplex(sel(i)).perm));
end
legend(leg,'Location','EastOutside');


function dC=dyneqn(c,A2,D2,A1,D1)
dC=zeros(size(c));
cc=c*c';
for i=1:length(c)
  dC(i)=c'*(A2{i}*c);   % Production from component parts
end
dC=dC ...
   -c.*(D2*c) ...    % Consumed in forming new complex
   -D1.*c ...	     % Degradation of complex to components
   +A1*c;	     % Formation due to breakdown of a complex
