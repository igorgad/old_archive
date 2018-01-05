function draw_states(t,i)

% written by StudentDave
%for licensing and usage questions
%email scienceguy5000 at gmail. com

%round off matrix values for easy of visualizing

t=round(t*100)/100;


%initialize variables for general displaying
syms C % make a symbolic variable that we'll put into ezplot
x=0.5*cos(C); % store it into it's function of interest
y=0.5*sin(C);
R = [ 0 2*pi]; % define the range
h=ezplot(x,y,[R]);
hold on
h2=ezplot(x+2,y,[R]);
h3=ezplot(x+4,y,[R]);
h4=ezplot(x+2,y-2,[R]);
arrow3([0.2 0],[1.8 0])
arrow3([1.8 0.3],[0.2 0.3],'--r' )
arrow3([3.8 0.3],[2.2 0.3],'--r' )
arrow3([2.2 0],[3.8 0])
arrow3([0 0],[1.8 -2.1])
arrow3([1.6 -2.1],[-0.2 0])
arrow3([2.1 0],[2.1 -2.1])
arrow3([1.9 -2.1],[1.9 0],'--r' )
arrow3([4.2 0],[2.4 -2.1])
arrow3([2.1 -2.1],[3.9 0],'--r' )


text(1.9,-2.7,num2str(t(1)),'FontSize',15)
text(.9,-0.9,num2str(t(2)),'FontSize',15)
text(2.1,-1,num2str(t(3)),'FontSize',15)
text(1.4,-1,num2str(t(9)),'FontSize',15,'color','r')
text(3.5,-1,num2str(t(4)),'FontSize',15)
text(2.7,-0.8,num2str(t(13)),'FontSize',15,'color','r')
text(.60,-1.4,num2str(t(5)),'FontSize',15)
text(.8,-0.1,num2str(t(10)),'FontSize',15)
text(2.8,-0.1,num2str(t(15)),'FontSize',15)
text(.8,0.4,num2str(t(7)),'FontSize',15','color','r')
text(2.8,-0.1,num2str(t(15)),'FontSize',15)
text(2.8,0.4,num2str(t(8)),'FontSize',15','color','r')

text(-0.15,0.6,num2str(t(6)),'FontSize',15','color','r')
text(1.85,0.6,num2str(t(11)),'FontSize',15','color','r')
text(3.85,0.6,num2str(t(16)),'FontSize',15','color','r')


xlabel('')
ylabel('')
title('')
axis off
set(h,'linewidth',5,'color','r')
set(h2,'linewidth',5,'color','y')
set(h3,'linewidth',5,'color','b')
set(h4,'linewidth',5,'color','k')

hold off

