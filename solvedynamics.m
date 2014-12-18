% Solve dynamics of concentration changes to find concentration of each complex after time t
% Args: c - structure returned by complexes()
% Returns: d - same as c but with concentrations updated to dynamic values
function d=solvedynamics(c,t,varargin)
defaults=struct('temp',55,'ka',1e6,'mindisplayconc',1e-12,'verbose',false,'reltol',0.01,'abstol',1e-12)
args=processargs(defaults,varargin);
ncomplex=length(c.complex);
% Setup up differential equation
% C(t)~ concentration of each complex at time t
% dC(k)/dt=sum(ka*C(i)*C(j)) - kd*C(m)
%    where complexes i,j can be combined to make complex k
%    or complex m can be broken down to from complex k
% Using fixed ka, find kd such that equilibrium stored in c.ocomplexconc will be reached
%   Or can use c.dG to find kd
A=zeros(ncomplex,ncomplex,ncomplex);  % A(i,j,k) is rate of association of i+j -> k; 
D=zeros(ncomplex,ncomplex,ncomplex);  % D(i,j,k) is rate of disocciation of k -> i+j
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
      A(i,j,k)=A(i,j,k)+args.ka/2;
      % Dissociation rate to reach equilibrium
      kd=args.ka*c.ocomplex(i).eqconc*c.ocomplex(j).eqconc/c.ocomplex(k).eqconc;
      D(i,j,k)=D(i,j,k)+kd/2;
      break;
    end
  end
end

initconds=c.concentrations;   % First complexes are strands by themselves
initconds(end+1:ncomplex)=0;
if args.verbose
  statsval='on';
else
  statsval='off';
end
options=odeset('Stats',statsval,'RelTol',args.reltol,'AbsTol',args.abstol,'NonNegative',1:length(initconds(:)),'NormControl','off'); % ,'OutputFcn',@odeplot);
[t,y]=ode23t(@(t,y) dyneqn(t,y,A,D),[0,t],initconds(:),options);
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


function dC=dyneqn(t,c,A,D)
dC=zeros(size(c));
for i=1:length(c)
  dC(i)=sum(sum((c*c').*A(:,:,i)))-sum(c(i)*c'*squeeze(A(i,:,:)))-sum(c(i)*c'*squeeze(A(:,i,:)))-sum(sum(D(:,:,i)))*c(i)+sum(squeeze(D(i,:,:))*c)+sum(squeeze(D(:,i,:))*c);
end
%total=[ 1 1 2 2 2]*c;
%fprintf(' C=[%s]\ndC=[%s]\n',sprintf('%.3g ',c),sprintf('%.3g ',dC));%,total);
                                                                              %if abs(total-4.1e-7)>1e-8
                                                                              %  keyboard
                                                                              %end
%set(gca,'YScale','log');
