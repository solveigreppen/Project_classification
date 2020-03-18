clear;

N = 100;

m1 = -1;
m2 = 1;

sigma1 = 0.25;
sigma2 = 0.49;
sigma3 = 1.0;

class11 = sigma1.*randn(N,1)+m1;
class12 = sigma1.*randn(N,1)+m2;

class21 = sigma2.*randn(N,1)+m1;
class22 = sigma2.*randn(N,1)+m2;

class31 = sigma3.*randn(N,1)+m1;
class32 = sigma3.*randn(N,1)+m2;


subplot(3,1,1)
histogram(class11)
hold on
histogram(class12);
title("variance = 0.25");

subplot(3,1,2)
histogram(class21)
hold on
histogram(class22);
title("variance = 0.49");

subplot(3,1,3)
histogram(class31)
hold on
histogram(class32);
title("variance = 1.00");

%NT = 5;
%newClass11 = zeros(NT);
%newClass12 = zeros(NT);
%{
for i = 0:NT-1
    newClass11(i) = rand(class11,1);
    newClass12(i) = rand(class12,1);
end
%}

%C1 = confusionmat(class11, class12);
%Cm = confusionchart(C1);

<<<<<<< Updated upstream
hei hallo hadet
print("Hello world");
=======
disp("hello world")
%hei
disp("hello world")

>>>>>>> Stashed changes
'Heisann'
'Hallo';