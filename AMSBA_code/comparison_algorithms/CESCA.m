% 混沌  SCA   Lin A, Wu Q, Heidari A A, Xu Y, Chen H, Geng W, li Y, Li C. Predicting Intentions of Students for Master Programs Using a Chaos-Induced Sine Cosine-Based Fuzzy K-Nearest Neighbor Classifier. IEEE Access, 2019, 7: 67235-67248.

function [Destination_position,Convergence_curve]=TentSCA(N,Max_iteration,lb,ub,dim,fobj)
tic

%Tent map
 ChaosVec(1)=0.7;
 for i=1:Max_iteration
     if ChaosVec(i)<0.7
         ChaosVec(i+1)=ChaosVec(i)/0.7;
     end
     if ChaosVec(i)>=0.7
         ChaosVec(i+1)=(10/3)*(1-ChaosVec(i));
     end
%      G(i)=(x(i))*Value;
 end

% ChaosVec(1)=0.7;
% for i=2:Max_iteration+1
%     ChaosVec(i)=4*ChaosVec(i-1)*(1-ChaosVec(i-1));
% end

%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);

Destination_position=zeros(1,dim);
Destination_fitness=inf;
t=0;
Convergence_curve=[];
Objective_values = zeros(1,size(X,1));
fes=0;
% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    Objective_values(1,i)=fobj(X(i,:));
    fes=fes+1;
    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
    
    All_objective_values(1,i)=Objective_values(1,i);
end
t=t+1;
Convergence_curve(t)=Destination_fitness;
%Main loop
 % start from the second iteration since the first iteration was dedicated to calculating the fitness
while fes<=Max_iteration
    
    % Eq. (3.4)
    a = 2;
    %     Max_iteration = Max_iteration;
    r1=a-fes*((a)/Max_iteration); % r1 decreases linearly from a to 0
    
    % Update the position of solutions with respect to destination
    for i=1:size(X,1) % in i-th solution
        for j=1:size(X,2) % in j-th dimension
            
            % Update r2, r3, and r4 for Eq. (3.3)
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();
            
            % Eq. (3.3)
            if r4<0.5
                % Eq. (3.1)
                X(i,j)= X(i,j)+(r1*sin(r2)*abs(ChaosVec(t)*Destination_position(j)-X(i,j)));
            else
                % Eq. (3.2)
                X(i,j)= X(i,j)+(r1*cos(r2)*abs(ChaosVec(t)*Destination_position(j)-X(i,j)));
            end
            % %% 混沌局部搜索
            % setCan = (Max_iteration-t+1)/Max_iteration;
            % %         setCan = 1- power(abs((t-1)/t),1000);
            % x = rand();
            % while(~(x~=0.25 && x~=0.5 && x~=0.75 && x~=1))
                % x=rand();
            % end
            % ch(1) = x;
            % for ant=1:N
                % ch(ant+1)=4*ch(ant)*(1-ch(ant));
                % CH(ant) = lb+ch(ant)*(ub-lb);    %ub大
                % V = (1-setCan)*Destination_position+setCan*CH(ant);
                % % ObjValV=feval(objfun,V);         %计算函数值
                
                % %% 边界控制
                % Flag4ub=V>ub';
                % Flag4lb=V<lb';
                % V=(V.*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;
                % %% 边界控制结束
                
                % FitnessV=fobj(V);%计算适应度值
                % %    if (ObjValV<GlobalMin)
                % if (FitnessV<Destination_fitness)
                    % Destination_fitness = FitnessV;
                    % Destination_position = V;
                    % %             break;
                % end
            % end
            %% 混沌最优值结束
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%           Moth_pos_m_gaus=X(i,:)*(1+randn(1));
%         Moth_fitness_m_gaus=fobj(Moth_pos_m_gaus);
%         Moth_fitness_s=fobj(X(i,:));
% %         Moth_fitness_f=fobj(Moth_pos(i,:));
%         
%         Moth_fitness_comb=[Moth_fitness_m_gaus,Moth_fitness_s];
%         [~,m]=min(Moth_fitness_comb);
%         if m==1        
%             X(i,:)=Moth_pos_m_gaus;
%             else
%             X(i,:)=X(i,:);
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
    end
    
    
    for i=1:size(X,1)
        
        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate the objective values
        Objective_values(1,i)=fobj(X(i,:));
        fes=fes+1;
        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    end
    % Increase the iteration counter
    t=t+1;
    Convergence_curve(t)=Destination_fitness;
    
    
end
toc
end

function o=Levy(d)
beta=1.5;
%Eq. (3.10)
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);
% Eq. (3.9)
o=step;
end