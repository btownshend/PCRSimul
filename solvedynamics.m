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
    end
  end
end
Dsum=squeeze(sum(sum(D,1)));
sD2=squeeze(sum(D,2));
As={};
sA3=zeros(size(A,1),size(A,2));
for k=1:size(A,3)
  As{k}=sparse(A(:,:,k));
  sA3=sA3+As{k};
end
clear A;

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
  dc1=dyneqn(jtmp,As,D);
  for i=1:length(initconds)
    jtmp(i)=0;
    dc2=dyneqn(jtmp,As,D);
    jtmp(i)=1;
    jpat(i,dc1~=dc2)=1;
  end
end

options=odeset('Stats',statsval,'RelTol',args.reltol,'AbsTol',args.abstol,'NonNegative',1:length(initconds),'NormControl','off'); % ,'JPattern',jpat); % ,'OutputFcn',@odeplot);
[t,y]=ode15s(@(t,y) dyneqn(y,As,sA3,D,sD2,Dsum),[0,t],initconds,options);
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


function dC=dyneqn(c,As,sA3,D,sD2,Dsum)
dC=zeros(size(c));
cc=c*c';
for i=1:length(c)
  %  d2=sum(sum(cc.*A(:,:,i)))-c(i)*sum(c'*squeeze(A(i,:,:)))*2-sum(sum(D(:,:,i)))*c(i)+sum(squeeze(D(i,:,:))*c)*2;
  %  d2=sum(sum(cc.*A(:,:,i)))-c(i)*dot(c,sum(A(i,:,:),3))*2-sum(sum(D(:,:,i)))*c(i)+dot(squeeze(sum(D(i,:,:),2)),c)*2;
  %  d2=sum(sum(cc.*A(:,:,i)))-c(i)*sA3c(i)*2-sum(sum(D(:,:,i)))*c(i)+sD2c(i)*2;
  %  dC(i)=sum(sum(cc.*A(:,:,i)));   % Production from component parts
  dC(i)=c'*(As{i}*c);   % Production from component parts
  %    if abs(d2-dC(i))>1e-15
  %      keyboard
  %    end
end
sA3c=sA3*c;
sD2c=sD2*c;
dC=dC ...
   -c.*sA3c*2 ...    % Consumed in forming new complex
   -Dsum.*c ...	     % Degradation of complex to components
   +sD2c*2;	     % Formation due to breakdown of a complex
%total=[ 1 1 2 2 2]*c;
%fprintf(' C=[%s]\ndC=[%s]\n',sprintf('%.3g ',c),sprintf('%.3g ',dC));%,total);
                                                                              %if abs(total-4.1e-7)>1e-8
                                                                              %  keyboard
                                                                              %end
%set(gca,'YScale','log');
