function [xo,fo,gap]=genetic(f,x0,l,u,Np,Nb,Pc,Pm,eta,kmax,model)
% Genetic Algorithm to minimize f(x) s.t. l<=x<=u
N=length(x0);
if nargin<10, kmax=100; end % Number of iterations (generations)
if nargin<9|eta>1|eta<=0, eta=1; end % Learning rate(0<eta<1)
if nargin<8, Pm=0.01; end % Probability of mutation
if nargin<7, Pc=0.5; end % Probability of crossover
if nargin<6, Nb=8*ones(1,N); end % # of genes(bits) for each variable
if nargin<5, Np=10; end % Population size(number of chromosomes)
% Initialize the population pool
NNb=sum(Nb);
xo=x0(:)'; l=l(:)'; u=u(:)';
% fo=feval(f,xo)
% 
[gmo]=forward_gravity(xo,model); %Fungsi Forward


fx1o=gmo;
Nf=length(fx1o);

fo=1/Nf*(sum(sqrt((fx1o-f).^2))); %Misfit 
X(1,:)=xo;
for n=2:Np
  X(n,:)=l+rand(size(x0)).*(u-l); % Eq.(7.1.26)
end

P=gen_encode(X,Nb,l,u); % Eq.(7.1.27)
gap = [];
for k=1:kmax
  X=gen_decode(P,Nb,l,u); % Eq.(7.1.28)
  for n=1:Np, 
%       
      [gm]=forward_gravity(X(n,:),model); %Fungsi Forward
%       
      fx1=gm; 
      Nf=length(fx1);
%       
      fX(n)=1/Nf*(sum(sqrt((fx1-f).^2))); %Misfit 
      gap(k)=fX(n);
    
  end
 
  [fxb,nb]=min(fX); % Selection of the fittest
  if fxb<fo, 
      fo=fxb; xo=X(nb,:); 

  end
  fX1=max(fX)-fX; %make the nonnegative fitness vector by Eq.(7.1.29)
  fXm=fX1(nb);
  if fXm<eps, return; end %terminate if all the chromosomes are equal
  % Reproduction of next generation
  for n=1:Np
     X(n,:)=X(n,:)+eta*(fXm-fX1(n))/fXm*(X(nb,:)-X(n,:)); %Eq.(7.1.30)
  end
  P=gen_encode(X,Nb,l,u);
  % Mating/Crossover
  is=shuffle([1:Np]);
  for n=1:2:Np-1
    if rand<Pc, P(is(n:n+1),:)=crossover(P(is(n:n+1),:),Nb); end
  end
% Mutation
  P=mutation(P,Nb,Pm);
end