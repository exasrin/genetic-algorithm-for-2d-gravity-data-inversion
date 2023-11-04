function P=mutation(P,Nb,Pm) % Mutation 
Nbb=length(Nb);
for n=1:size(P,1)
   b2=0;
   for m=1:Nbb
      if rand<Pm
        b1=b2+1; bi=b1+mod(floor(rand*Nb(m)),Nb(m)); b2=b2+Nb(m);
        P(n,bi)=num2str(1-eval(P(n,bi)));
      end 
   end
end