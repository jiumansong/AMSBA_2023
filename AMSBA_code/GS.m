function [v]=GS(A)

v(:,1)=A(:,1)/norm(A(:,1));%归一化
[Ahang,Alie]=size(A); %矩阵的行和列
for k=2:Alie %求解第j列正交向量
    res(:,1)=A(:,k);
    for i=1:k-1%减去待求向量在以求向量上的投影
        res=res-v(:,i)'*res*v(:,i);
    end
    v(:,k)=res/norm(res);%归一化
end
end