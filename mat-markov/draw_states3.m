function draw_states3(t,i)

% written by StudentDave
%for licensing and usage questions
%email scienceguy5000 at gmail. com


t=round(t*100)/100;
%initialize variables for general displaying
syms C % make a symbolic variable that we'll put into ezplot
x=0.5*cos(C); % store it into it's function of interest
y=0.5*sin(C);
R = [ 0 2*pi]; % define the range
h=ezplot(x,y,[R]);
hold on
h2=ezplot(x+2,y+2,[R]);
h3=ezplot(x+4,y,[R]);
arrow3([0 .25],[2 2.25])
arrow3([2 1.75],[0 -.25])
arrow3([2.2 2.25],[4 .25])
arrow3([4 -.25],[2.2 1.75])
arrow3([0 .15],[4 .15])
arrow3([4 -.35],[0 -.35])

text(0.5,1.3,num2str(t(4)),'color','k','FontSize',15)
text(1.2,0.8,num2str(t(2)),'color','k','FontSize',15)
text(2.5,0.8,num2str(t(6)),'color','k','FontSize',15)
text(3.2,1.3,num2str(t(8)),'color','k','FontSize',15)
text(1.8,.3,num2str(t(7)),'color','k','FontSize',15)
text(1.8,-.5,num2str(t(3)),'color','k','FontSize',15)

text(-.2,-.7,num2str(t(1)),'color','k','FontSize',15)
text(3.8,-.7,num2str(t(9)),'color','k','FontSize',15)
text(1.8,2.7,num2str(t(5)),'color','k','FontSize',15)

xlabel(['Time steps = ', num2str(i)])
ylabel('')
title('')
axis off
set(h,'linewidth',5,'color','r')
set(h2,'linewidth',5,'color','y')
set(h3,'linewidth',5,'color','b')

hold off

