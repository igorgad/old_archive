clc
clear
M=7;
N=7;
LX=1;
LY=1;
dx=LX/M;
dy=LY/N;
aa=-1;

for i=1:M;
aa=aa+1;
x(i)=aa*dx;
end
%figure('Position',[20 20 1400 700])

aa=-1;

for j=1:N;
aa=aa+1;
y(j)=aa*dy;
end

for i=1:M;
for j=1:N;
%plot(i,j,'*r');
hold on
end
end

plot([8 8],[8 8]);
plot([0 0],[0 0]);

for i=1:M;
for j=1:N;
u(j,i)=y(j);
end
end

for i=1:M;
for j=1:N;
v(j,i)=(1/i)*(cos(x(i))-sin(x(j)));
w_angle(j,i) = atan(v(j,i)/u(j,i));
w_mag(j,i) = sqrt((v(j,i)^2)+(u(j,i)^2));
end
end

quiver(u,v,'b')
grid on
xlabel('x');
ylabel('y')
title('Vector Field')
x = [1 2 2 2 3 4 5 6 7];
y = [1 2 3 4 5 6 7 7 7];
plot(x,y, 'black')
%set(gca,'XLim',[0 M+2],'YLim',[0 N])

sigma = 0.5;
Vmin = 20;
Wmax = 15;
C = 50;
theta_a = 0;
v_angle = 0;
lambda = 0.5;

A = [0 45 90 135 180 -135 -90 -45];
A = degtorad(A);

i_obj = 7; j_obj = 7;

w_mag = w_mag*Wmax; %velocidade do vento em cada célula

%montar a matriz de probabilidades
Paux1 = zeros(49,49);
Paux2 = zeros(49,49);
Paux3 = zeros(49,49);
Paux4 = zeros(49,49);
Paux5 = zeros(49,49);
Paux6 = zeros(49,49);
Paux7 = zeros(49,49);
Paux8 = zeros(49,49);
Paux = zeros(49,49);
Ra1 = zeros(1,49);
Ra2 = zeros(1,49);
Ra3 = zeros(1,49);
Ra4 = zeros(1,49);
Ra5 = zeros(1,49);
Ra6 = zeros(1,49);
Ra7 = zeros(1,49);
Ra8 = zeros(1,49);


for i = 1:7
    for j = 1:7
        for k = 1:8
            Fx = Vmin*cos(A(1,k)) + w_mag(i,j)*cos(w_angle(i,j));
            Fy = Vmin*sin(A(1,k)) + w_mag(i,j)*sin(w_angle(i,j));
            F = sqrt((Fx^2)+(Fy^2));
            w = atan(Fy/Fx);
            mu = w;
            
            lim = [A(1)-(pi/8) A(1)+(pi/8)];
            gauss = normcdf(lim,mu,sigma);
            p1 = gauss(2)-gauss(1);
            
            lim = [A(2)-(pi/8) A(2)+(pi/8)];
            gauss = normcdf(lim,mu,sigma);
            p2 = gauss(2)-gauss(1);
            
            lim = [A(3)-(pi/8) A(3)+(pi/8)];
            gauss = normcdf(lim,mu,sigma);
            p3 = gauss(2)-gauss(1);
            
            lim = [A(4)-(pi/8) A(4)+(pi/8)];
            gauss = normcdf(lim,mu,sigma);
            p4 = gauss(2)-gauss(1);
            
            lim = [A(5)-(pi/8) A(5)+(pi/8)];
            gauss = normcdf(lim,mu,sigma);
            p5 = gauss(2)-gauss(1);
            
            lim = [A(6)-(pi/8) A(6)+(pi/8)];
            gauss = normcdf(lim,mu,sigma);
            p6 = gauss(2)-gauss(1);
            
            lim = [A(7)-(pi/8) A(7)+(pi/8)];
            gauss = normcdf(lim,mu,sigma);
            p7 = gauss(2)-gauss(1);
            
            lim = [A(8)-(pi/8) A(8)+(pi/8)];
            gauss = normcdf(lim,mu,sigma);
            p8 = gauss(2)-gauss(1);
            
            estado = (7*(i-1))+j;
            
            if (k == 1)
                if (mod(estado,7) ~= 0)
                    Paux1(estado,estado+1) = p1;
                end
                if (mod(estado,7) ~= 0) && (floor(estado/7) < 6)
                    Paux1(estado,estado+8) = p2;
                end
                if (floor(estado/7) < 6)
                    Paux1(estado,estado+7) = p3;
                end                
                if (floor(estado/7) < 6) && (mod(estado,7) ~= 1)
                    Paux1(estado,estado+6) = p4;
                end                
                if (mod(estado,7) ~= 1)
                    Paux1(estado,estado-1) = p5;
                    theta_t = atan((i_obj-i)/(j_obj-j));
                    Ra1(1,estado) = ((w_mag(i,j)*cos(w_angle(i,j)+theta_t))/Wmax)*C;
                end                
                if (mod(estado,7) ~= 1) && (floor(estado/7) > 1)
                    Paux1(estado,estado-8) = p6;
                end                
                if (floor(estado/7) > 1)
                    Paux1(estado,estado-7) = p7;
                end                
                if (floor(estado/7) > 1) && (mod(estado,7) ~= 0)
                    Paux1(estado,estado-6) = p8;
                end
            end
            
            if (k == 2)
                if (mod(estado,7) ~= 0)
                    Paux2(estado,estado+1) = p1;
                end                
                if (mod(estado,7) ~= 0) && (floor(estado/7) < 6)
                    Paux2(estado,estado+8) = p2;
                end                
                if (floor(estado/7) < 6)
                    Paux2(estado,estado+7) = p3;
                end                
                if (floor(estado/7) < 6) && (mod(estado,7) ~= 1)
                    Paux2(estado,estado+6) = p4;
                end                
                if (mod(estado,7) ~= 1)
                    Paux2(estado,estado-1) = p5;
                end                
                if (mod(estado,7) ~= 1) && (floor(estado/7) > 1)
                    Paux2(estado,estado-8) = p6;
                    theta_t = atan((i_obj-i)/(j_obj-j));
                    Ra2(1,estado) = ((w_mag(i,j)*cos(w_angle(i,j)+theta_t))/Wmax)*C;
                end                
                if (floor(estado/7) > 1)
                    Paux2(estado,estado-7) = p7;
                end                
                if (floor(estado/7) > 1) && (mod(estado,7) ~= 0)
                    Paux2(estado,estado-6) = p8;
                end
            end

            if (k == 3)
                if (mod(estado,7) ~= 0)
                    Paux3(estado,estado+1) = p1;
                end                
                if (mod(estado,7) ~= 0) && (floor(estado/7) < 6)
                    Paux3(estado,estado+8) = p2;
                end                
                if (floor(estado/7) < 6)
                    Paux3(estado,estado+7) = p3;
                end                
                if (floor(estado/7) < 6) && (mod(estado,7) ~= 1)
                    Paux3(estado,estado+6) = p4;
                end                
                if (mod(estado,7) ~= 1)
                    Paux3(estado,estado-1) = p5;
                end                
                if (mod(estado,7) ~= 1) && (floor(estado/7) > 1)
                    Paux3(estado,estado-8) = p6;
                end                
                if (floor(estado/7) > 1)
                    Paux3(estado,estado-7) = p7;
                    theta_t = atan((i_obj-i)/(j_obj-j));
                    Ra3(1,estado) = ((w_mag(i,j)*cos(w_angle(i,j)+theta_t))/Wmax)*C;
                end                
                if (floor(estado/7) > 1) && (mod(estado,7) ~= 0)
                    Paux3(estado,estado-6) = p8;
                end
            end

            if (k == 4)
                if (mod(estado,7) ~= 0)
                    Paux4(estado,estado+1) = p1;
                end                
                if (mod(estado,7) ~= 0) && (floor(estado/7) < 6)
                    Paux4(estado,estado+8) = p2;
                end                
                if (floor(estado/7) < 6)
                    Paux4(estado,estado+7) = p3;
                end                
                if (floor(estado/7) < 6) && (mod(estado,7) ~= 1)
                    Paux4(estado,estado+6) = p4;
                end                
                if (mod(estado,7) ~= 1)
                    Paux4(estado,estado-1) = p5;
                end                
                if (mod(estado,7) ~= 1) && (floor(estado/7) > 1)
                    Paux4(estado,estado-8) = p6;
                end                
                if (floor(estado/7) > 1)
                    Paux4(estado,estado-7) = p7;
                end                
                if (floor(estado/7) > 1) && (mod(estado,7) ~= 0)
                    Paux4(estado,estado-6) = p8;
                    theta_t = atan((i_obj-i)/(j_obj-j));
                    Ra4(1,estado) = ((w_mag(i,j)*cos(w_angle(i,j)+theta_t))/Wmax)*C;
                end
            end

            if (k == 5)
                if (mod(estado,7) ~= 0)
                    Paux5(estado,estado+1) = p1;
                    theta_t = atan((i_obj-i)/(j_obj-j));
                    Ra5(1,estado) = ((w_mag(i,j)*cos(w_angle(i,j)+theta_t))/Wmax)*C;
                end                
                if (mod(estado,7) ~= 0) && (floor(estado/7) < 6)
                    Paux5(estado,estado+8) = p2;
                end                
                if (floor(estado/7) < 6)
                    Paux5(estado,estado+7) = p3;
                end                
                if (floor(estado/7) < 6) && (mod(estado,7) ~= 1)
                    Paux5(estado,estado+6) = p4;
                end                
                if (mod(estado,7) ~= 1)
                    Paux5(estado,estado-1) = p5;
                end                
                if (mod(estado,7) ~= 1) && (floor(estado/7) > 1)
                    Paux5(estado,estado-8) = p6;
                end                
                if (floor(estado/7) > 1)
                    Paux5(estado,estado-7) = p7;
                end                
                if (floor(estado/7) > 1) && (mod(estado,7) ~= 0)
                    Paux5(estado,estado-6) = p8;
                end
            end

            if (k == 6)
                if (mod(estado,7) ~= 0)
                    Paux6(estado,estado+1) = p1;
                end                
                if (mod(estado,7) ~= 0) && (floor(estado/7) < 6)
                    Paux6(estado,estado+8) = p2;
                    theta_t = atan((i_obj-i)/(j_obj-j));
                    Ra6(1,estado) = ((w_mag(i,j)*cos(w_angle(i,j)+theta_t))/Wmax)*C;
                end                
                if (floor(estado/7) < 6)
                    Paux6(estado,estado+7) = p3;
                end                
                if (floor(estado/7) < 6) && (mod(estado,7) ~= 1)
                    Paux6(estado,estado+6) = p4;
                end                
                if (mod(estado,7) ~= 1)
                    Paux6(estado,estado-1) = p5;
                end                
                if (mod(estado,7) ~= 1) && (floor(estado/7) > 1)
                    Paux6(estado,estado-8) = p6;
                end                
                if (floor(estado/7) > 1)
                    Paux6(estado,estado-7) = p7;
                end                
                if (floor(estado/7) > 1) && (mod(estado,7) ~= 0)
                    Paux6(estado,estado-6) = p8;
                end
            end

            if (k == 7)
                if (mod(estado,7) ~= 0)
                    Paux7(estado,estado+1) = p1;
                end                
                if (mod(estado,7) ~= 0) && (floor(estado/7) < 6)
                    Paux7(estado,estado+8) = p2;
                end                
                if (floor(estado/7) < 6)
                    Paux7(estado,estado+7) = p3;
                    theta_t = atan((i_obj-i)/(j_obj-j));
                    Ra7(1,estado) = ((w_mag(i,j)*cos(w_angle(i,j)+theta_t))/Wmax)*C;
                end                
                if (floor(estado/7) < 6) && (mod(estado,7) ~= 1)
                    Paux7(estado,estado+6) = p4;
                end                
                if (mod(estado,7) ~= 1)
                    Paux7(estado,estado-1) = p5;
                end                
                if (mod(estado,7) ~= 1) && (floor(estado/7) > 1)
                    Paux7(estado,estado-8) = p6;
                end                
                if (floor(estado/7) > 1)
                    Paux7(estado,estado-7) = p7;
                end                
                if (floor(estado/7) > 1) && (mod(estado,7) ~= 0)
                    Paux7(estado,estado-6) = p8;
                end
            end

            if (k == 8)
                if (mod(estado,7) ~= 0)
                    Paux8(estado,estado+1) = p1;
                end                
                if (mod(estado,7) ~= 0) && (floor(estado/7) < 6)
                    Paux8(estado,estado+8) = p2;
                end                
                if (floor(estado/7) < 6)
                    Paux8(estado,estado+7) = p3;
                end                
                if (floor(estado/7) < 6) && (mod(estado,7) ~= 1)
                    Paux8(estado,estado+6) = p4;
                    theta_t = atan((i_obj-i)/(j_obj-j));
                    Ra8(1,estado) = ((w_mag(i,j)*cos(w_angle(i,j)+theta_t))/Wmax)*C;
                end                
                if (mod(estado,7) ~= 1)
                    Paux8(estado,estado-1) = p5;
                end                
                if (mod(estado,7) ~= 1) && (floor(estado/7) > 1)
                    Paux8(estado,estado-8) = p6;
                end                
                if (floor(estado/7) > 1)
                    Paux8(estado,estado-7) = p7;
                end                
                if (floor(estado/7) > 1) && (mod(estado,7) ~= 0)
                    Paux8(estado,estado-6) = p8;
                end
            end            
        end
    end
end

for i = 1:49
    sobra = 1 - sum(Paux1(i,:));
    [l,c] = max(Paux1(i,:));
    Paux1(i,c) = Paux1(i,c) + sobra;
end

for i = 1:49
    sobra = 1 - sum(Paux2(i,:));
    [l,c] = max(Paux2(i,:));
    Paux2(i,c) = Paux2(i,c) + sobra;
end

for i = 1:49
    sobra = 1 - sum(Paux3(i,:));
    [l,c] = max(Paux3(i,:));
    Paux3(i,c) = Paux3(i,c) + sobra;
end

for i = 1:49
    sobra = 1 - sum(Paux4(i,:));
    [l,c] = max(Paux4(i,:));
    Paux4(i,c) = Paux4(i,c) + sobra;
end

for i = 1:49
    sobra = 1 - sum(Paux5(i,:));
    [l,c] = max(Paux5(i,:));
    Paux5(i,c) = Paux5(i,c) + sobra;
end

for i = 1:49
    sobra = 1 - sum(Paux6(i,:));
    [l,c] = max(Paux6(i,:));
    Paux6(i,c) = Paux6(i,c) + sobra;
end

for i = 1:49
    sobra = 1 - sum(Paux7(i,:));
    [l,c] = max(Paux7(i,:));
    Paux7(i,c) = Paux7(i,c) + sobra;
end

for i = 1:49
    sobra = 1 - sum(Paux8(i,:));
    [l,c] = max(Paux8(i,:));
    Paux8(i,c) = Paux8(i,c) + sobra;
end


P(:,:,1) = Paux1;
P(:,:,2) = Paux2;
P(:,:,3) = Paux3;
P(:,:,4) = Paux4;
P(:,:,5) = Paux5;
P(:,:,6) = Paux6;
P(:,:,7) = Paux7;
P(:,:,8) = Paux8;

Ra1(1,49) = Ra1(1,48);
Ra2(1,49) = Ra2(1,48);
Ra3(1,49) = Ra3(1,48);
Ra4(1,49) = Ra4(1,48);
Ra5(1,49) = Ra5(1,48);
Ra6(1,49) = Ra6(1,48);
Ra7(1,49) = Ra7(1,48);
Ra8(1,49) = Ra8(1,48);

R(:,1) = Ra1';
R(:,2) = Ra2';
R(:,3) = Ra3';
R(:,4) = Ra4';
R(:,5) = Ra5';
R(:,6) = Ra6';
R(:,7) = Ra7';
R(:,8) = Ra8';

mdp_check(P,R)
discount = 0.01;
[V, policy] = mdp_policy_iteration(P, R, discount);
[policy2] = mdp_value_iteration(P,R,discount);
[V3, policy3] = mdp_LP(P, R, discount);
[~, V4, policy4] = mdp_Q_learning(P, R, discount);

% V = zeros(7,7);
% Vaux = zeros(1,8);
% Prob = zeros(1,8);
% Raaux = zeros(1,8);
% Q = zeros(1,8);
% pol = zeros(7,7);
% 
% for k = 1:10
%     for i = 1:7
%         for j = 1:7
%             estado = (7*(i-1))+j;
%             for z = 1:8
%                 if (mod(estado,7) ~= 0)
%                     Vaux(1,1) = V(i,j+1);
%                     Raaux(1,1) = R(estado,z);
%                     Paux = P(:,:,z);
%                     Prob(1,1) = Paux(i,j+1);
%                 end
%                 if (mod(estado,7) ~= 0) && (floor(estado/7) < 6)
%                     Vaux(1,2) = V(i+1,j+1);
%                     Raaux(1,2) = R(estado,z);
%                     Paux = P(:,:,z);
%                     Prob(1,2) = Paux(i+1,j+1);
%                 end
%                 if (floor(estado/7) < 6)
%                     Vaux(1,3) = V(i+1,j);
%                     Raaux(1,3) = R(estado,z);
%                     Paux = P(:,:,z);
%                     Prob(1,3) = Paux(i+1,j);
%                 end
%                 if (floor(estado/7) < 6) && (mod(estado,7) ~= 1)
%                     Vaux(1,4) = V(i+1,j-1);
%                     Raaux(1,4) = R(estado,z);
%                     Paux = P(:,:,z);
%                     Prob(1,4) = Paux(i+1,j-1);
%                 end
%                 if (mod(estado,7) ~= 1)
%                     Vaux(1,5) = V(i,j-1);
%                     Raaux(1,5) = R(estado,z);
%                     Paux = P(:,:,z);
%                     Prob(1,5) = Paux(i,j-1);
%                 end
%                 if (mod(estado,7) ~= 1) && (floor(estado/7) > 1)
%                     Vaux(1,6) = V(i-1,j-1);
%                     Raaux(1,6) = R(estado,z);
%                     Paux = P(:,:,z);
%                     Prob(1,6) = Paux(i-1,j-1);
%                 end
%                 if (floor(estado/7) > 1)
%                     Vaux(1,7) = V(i-1,j);
%                     Raaux(1,7) = R(estado,z);
%                     Paux = P(:,:,z);
%                     Prob(1,7) = Paux(i-1,j);
%                 end
policy(1) = 2;
policy(9) = 3;
policy(16) = 2;
policy(24) = 3;
policy(31) = 2;
policy(39) = 2;
policy(47) = 1;
policy(48) = 1;
policy(49) = 1;
%                 if (floor(estado/7) > 1) && (mod(estado,7) ~= 0)
%                     Vaux(1,8) = V(i-1,j+1);
%                     Raaux(1,8) = R(estado,z);
%                     Paux = P(:,:,z);
%                     Prob(1,8) = Paux(i-1,j+1);
%                 end
%                 soma = 0;
%                 for q = 1:8
%                     soma = soma + Prob(1,q)*Vaux(1,q);
%                 end
%                 soma = soma * lambda;
%                 Q(1,z) = R(estado,z) + soma;
%             end
%             [l, c] = max(Q);
%             pol(i,j) = A(c);
%             V(i,j) = Vaux(1,c);
%             v_angle = A(c);
%         end
%     end
% end