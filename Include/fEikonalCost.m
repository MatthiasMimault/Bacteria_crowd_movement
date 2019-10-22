function [Phi] = fEikonalCost(X,Y,C,PhiBd,PhiAt)
% CAPS_EIKONAL
% Solve the eikonal equation with the fast sweeping method with a cost C

%% Parameters
iter_MAX = 5;  

%% Variable
X1 = min(min(X));
X2 = max(max(X));
Y1 = min(min(Y));
Y2 = max(max(Y));
dx = X(1,2)-X(1,1);
dy = Y(1,2)-Y(1,1);
dPhi = sqrt(dx^2+dy^2);
MaxP = sqrt((X2-X1).^2+(Y2-Y1).^2);

%% Boundary conditions
Phi = MaxP*(ones(size(X)) - PhiAt);
[Lx,Ly]=size(Phi);  

%While loop on the error of XXX
% while e>e_thr 
for iter = 1:iter_MAX    
    %First loop
    for i = 2:Lx-1                 
        for j = 2:Ly-1
            a = min(Phi(i,j-1),Phi(i,j+1));
            b = min(Phi(i-1,j),Phi(i+1,j));
            
            %Evaluation of the cost function   
            cost = C(i,j);

%             %Computation of Phibar
%             if abs(a-b)<dPhi 
%                 Phibar = (a+b+sqrt(2*dPhi^2-(a-b)^2))/2;
%             else
%                 Phibar = min(a,b)+dPhi;
%             end
            %Computation of Phibar
            if abs(a-b)<cost*dPhi 
                Phibar = (a+b+sqrt(2*(cost*dPhi)^2-(a-b)^2))/2;
            else
                Phibar = min(a,b)+cost*dPhi;
            end
            
            %Update
            Phi(i,j) = min(Phibar,Phi(i,j)).*(1-PhiBd(i,j)>0)+MaxP.*PhiBd(i,j); 
        end
    end
    
    
    
    %Second loop
    for i = Lx-1:-1:2               
        for j = 2:Ly-1
            a = min(Phi(i,j-1),Phi(i,j+1));
            b = min(Phi(i-1,j),Phi(i+1,j));
            
            %Evaluation of the cost function   
            cost = C(i,j);  

%             %Computation of Phibar
%             if abs(a-b)<dPhi 
%                 Phibar = (a+b+sqrt(2*dPhi^2-(a-b)^2))/2;
%             else
%                 Phibar = min(a,b)+dPhi;
%             end
            %Computation of Phibar
            if abs(a-b)<cost*dPhi 
                Phibar = (a+b+sqrt(2*(cost*dPhi)^2-(a-b)^2))/2;
            else
                Phibar = min(a,b)+cost*dPhi;
            end
            
            %Update
            Phi(i,j) = min(Phibar,Phi(i,j)).*(1-PhiBd(i,j)>0)+MaxP.*PhiBd(i,j);
        end
    end
    
    
    %Third loop
    for i = Lx-1:-1:2                
        for j = Ly-1:-1:2
            a = min(Phi(i,j-1),Phi(i,j+1));
            b = min(Phi(i-1,j),Phi(i+1,j)); 
            
            %Evaluation of the cost function   
            cost = C(i,j);  

%             %Computation of Phibar
%             if abs(a-b)<dPhi 
%                 Phibar = (a+b+sqrt(2*dPhi^2-(a-b)^2))/2;
%             else
%                 Phibar = min(a,b)+dPhi;
%             end
            %Computation of Phibar
            if abs(a-b)<cost*dPhi 
                Phibar = (a+b+sqrt(2*(cost*dPhi)^2-(a-b)^2))/2;
            else
                Phibar = min(a,b)+cost*dPhi;
            end
            
            %Update
            Phi(i,j) = min(Phibar,Phi(i,j)).*(1-PhiBd(i,j)>0)+MaxP.*PhiBd(i,j); 
        end
    end
    
    
    %Fourth loop
    for i = 2:Lx-1                
        for j = Ly-1:-1:2
            a = min(Phi(i,j-1),Phi(i,j+1));
            b = min(Phi(i-1,j),Phi(i+1,j)); 
            
            %Evaluation of the cost function   
            cost = C(i,j);  

%             %Computation of Phibar
%             if abs(a-b)<dPhi 
%                 Phibar = (a+b+sqrt(2*dPhi^2-(a-b)^2))/2;
%             else
%                 Phibar = min(a,b)+dPhi;
%             end
            %Computation of Phibar
            if abs(a-b)<cost*dPhi 
                Phibar = (a+b+sqrt(2*(cost*dPhi)^2-(a-b)^2))/2;
            else
                Phibar = min(a,b)+cost*dPhi;
            end
            
            %Update
            Phi(i,j) = min(Phibar,Phi(i,j)).*(1-PhiBd(i,j)>0)+MaxP.*PhiBd(i,j);
        end
    end
    
    % Boundary
%     Phi = Phi.*(1-PhiBd>0)+2*sqrt((X2-X1).^2+(Y2-Y1).^2).*PhiBd;
%     surf(X,Y,Phi, 'EdgeColor','none')
%     figure
%     surf(X,Y,Phi.*(1-PhiBd-PhiAt), 'EdgeColor','none')
%     figure
%     surf(X,Y,2*sqrt((X2-X1).^2+(Y2-Y1).^2).*PhiBd, 'EdgeColor','none')
% %     Phi = Phi.*(1-PhiBd-PhiAt)+2*sqrt((X2-X1).^2+(Y2-Y1).^2).*PhiBd;
% %     Phi = Phi.*(1-PhiBd-PhiAt)+max(max(Phi.*(1-PhiBd-PhiAt))).*PhiBd;
%     figure    
end 

