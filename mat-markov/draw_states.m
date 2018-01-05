function draw_states(t,i)

% written by StudentDave
%for licensing and usage questions
%email scienceguy5000 at gmail. com

%initialize variables for general displaying
t = round(t*100)/100;
syms C % make a symbolic variable that we'll put into ezplot
x=0.5*cos(C); % store it into it's function of interest
y=0.5*sin(C);
R = [ 0 2*pi]; % define the range
h=ezplot(x,y,[R]);
hold on
h2=ezplot(x+2,y,[R]);
arrow3([0 .25],[2 .25])
arrow3([2 -.25],[0 -.25])
text(-0.1,.6,num2str(t(1)),'color','b','FontSize',15)
text(.85,-.3,num2str(t(2)),'color','g','FontSize',15)
text(.85,.3,num2str(t(3)),'color','r','FontSize',15)
text(1.9,0.6,num2str(t(4)),'color','c','FontSize',15)

xlabel('')
ylabel('')
title(['Time steps = ', num2str(i)])
axis off
set(h,'linewidth',5,'color','k')
set(h2,'linewidth',5,'color','k')

hold off

