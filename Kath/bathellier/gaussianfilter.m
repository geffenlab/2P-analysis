function [S NormFactor]=gaussianfilter(S,dt,sig,Norm)
N=round(3*sig/dt);
t=(-N:N)*dt;
G=exp(-t.^2/(2*sig^2)); 
S0=S;
% Normalization
switch(Norm)
    case(1)
        NormFactor=[sum(G)/sqrt(sum(G.^2)) sum(G)/sqrt(sum(diff(G).^2))];
        G=G/sum(G);
    case(2)
        G=G/sqrt(sum(G.^2));
        NormFactor=sqrt(sum(G.^2))/sqrt(sum(diff(G).^2));
end

S=cat(2,S(:,N:-1:1,:,:),S,S(:,end:-1:end-N+1,:,:));
S=filter(G,1,S,[],2);
S=S(:,2*N+1:end,:,:);

