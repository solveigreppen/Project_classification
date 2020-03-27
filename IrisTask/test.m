%Class 1 = Iris Setosa
%Class 2 = Iris Versicolour
%Class 3 = Iris Virginica


clear;
Ntrain = 30;
Ntest = 20;
C = 3; % # classes
D = 4; % feature dimension
M = 1000; % # iterations (velger vi denne størrelsen selv?)
alpha = 0.05; %teste seg fram med denne
W0 = zeros(C,D); %skal vi fylle noe i denne?
w0 = zeros(C,1); %hva skal vi fylle i denne?
W = [W0 w0]; %forenkler W til senere uttrykk
nablaW_MSEs = zeros(1,M);
MSEs= zeros(1,M); %creating space for the matrix

%load the datas
x1all = load('class_1','-ascii');
x2all = load('class_2','-ascii');
x3all = load('class_3','-ascii');

% petal widths in cm
x1= [x1all(:,4)];
x2= [x2all(:,4)];
x3= [x3all(:,4)];

% train sets
x1_train = x1(1:Ntrain);
x2_train = x2(1:Ntrain);
x3_train = x3(1:Ntrain);
x_train = [x1_train, x2_train, x3_train];

% test sets
x1_test = x1(Ntrain+1:end);
x2_test = x2(Ntrain+1:end);
x3_test = x3(Ntrain+1:end);
x_test = [x1_test, x2_test, x3_test];

%training the classifier
for m = 1:M
    nablaW_MSE = 0;
    MSE=0; 
    
    for k = 1:Ntrain
        xk = x_train(k,:).'; %legge in x = [xk' 1]' ?
        zk = W.*xk; %forenkle til zk = Wx?
    
        gk = (1+exp(-zk)).^-1; %bruke innebygd sigmoid eller lage egen funksjon for å forenkle koden?
        
        %kopiert kode, bør endres
        tk = zeros(C,1);
        c = floor((k-1)/Ntrain * C) + 1;
        tk(c) = 1;
        
        MSE1 = (gk-tk).*(gk).*(1-gk);     
        nablaW_MSE = nablaW_MSE + MSE1.*xk;
        MSE = MSE + 0.5*(gk-tk)'*(gk-tk);
        
    end
    W = W - alpha.*nablaW_MSE;
    MSEs(m,1:5) = MSE;
    nablaW_MSEs(m) = norm(nablaW_MSE); %riktig å gjøre det sånn?
end
%}

%Plots 

plot(nablaW_MSEs); 
title('MSE gradient'); 
ylabel('Magnitude');
xlabel('iterations');