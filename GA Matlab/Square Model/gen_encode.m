function P=gen_encode(X,Nb,l,u)
% fprintf('ini: %d \n',X);
% Encode a population(X) of state into an array(P) of binary strings
Np=size(X,1); % Population size
N=length(Nb); % Dimension of the variable(state)
for n=1:Np
   b2=0;
   for m=1:N
      b1=b2+1; b2=b2+Nb(m);
      Xnm=(2^Nb(m)-1)*(X(n,m)-l(m))/(u(m)-l(m)); %Eq.(7.1.27)
      P(n,b1:b2)=dec2bin(Xnm,Nb(m)); %encoding to binary strings
   end 
end