function [ BestFitness,BestPosiiton ] = MRM( X0,Posiiton,N,dim,fobj )
%UNTITLED2 此处显示有关此函数的摘要
%   X0待优化个体位置,Posiiton种群个体位置,N种群大小,dim问题维度,fobj适应度函数
Alfa=20;
Beta=-0.5;
NRM=200;
M=0.1*N;
Sigma=zeros(1,dim);
Fitness=zeros(1,N);
dis=zeros(M,dim);
Epsilon1=1.0e-150;
Epsilon2=1.0e-4;
Nl=15;
Xk=X0;
X=X0;
k=0;
k2=0;
for i=1:N
    Fitness(1,i)=fobj(Posiiton(i,:));
end
[sorted_fitness,sorted_indexes]=sort(Fitness);
for i=1:M
    for j=1:dim
        dis(i,j)=Posiiton(sorted_indexes(i),j)-X0(j);
    end
end
Sigma=0.1*sum(dis)/M;
Df=eye(dim,dim);
while((k2<2*Nl)&&(min(Sigma)>=1.0*10^(-(Alfa+NRM))))
    X=Xk;
    k1=0;
    Z=X;
    while(k1<Nl)
        Df=Df';
        for j=1:dim
            Y=X+Sigma(j)*Df(j,:);
            if fobj(Y)<fobj(X)
                X=Y;
                Sigma(j)=Alfa*Sigma(j);
            else
                Sigma(j)=Beta*Sigma(j);
            end
        end
        Df=Df';
        if (abs(fobj(Z)-fobj(X))/abs(fobj(X)+Epsilon1))<Epsilon2
            k2=k2+1;
        else
            k2=0;
        end
        k1=k1+1;
    end
    if ((fobj(X)<fobj(Xk))&&(min(Sigma)>=1e-28))
        k=k+1;
        Xk=X;
        P=zeros(dim,dim);
        for j=1:dim
            if Sigma(j)==0
                P(:,j)=Df(:,j);
            else
                for m=j:dim
                    P(:,j)=P(:,j)+Sigma(m)*Df(:,m);
                end
            end
        end
        [v]=GS(P);
        Df=v;
    end
end
BestPosiiton=X;
BestFitness=fobj(X);
end

