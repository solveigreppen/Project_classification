%Class 1 = Iris Setosa
%Class 2 = Iris Versicolour
%Class 3 = Iris Virginica


clear;
Ntrain = 30;
Ntest = 20;
C = 3; % # classes
D = 4; % feature dimension
M = 1000; % # iterations (velger vi denne st�rrelsen selv?)
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
x_train = [x1_train; x2_train; x3_train];

% test sets
x1_test = x1(Ntrain+1:end);
x2_test = x2(Ntrain+1:end);
x3_test = x3(Ntrain+1:end);
x_test = [x1_test; x2_test; x3_test];

%training the classifier
for m = 1:M
    nablaW_MSE = 0;
    MSE=0; 
    
    for k = 1:Ntrain
        xk = x_train(k,:).'; %legge in x = [xk' 1]' ?
        zk = W.*xk; %forenkle til zk = Wx?
    
        gk = (1+exp(-zk)).^-1; %bruke innebygd sigmoid eller lage egen funksjon for � forenkle koden?
        
        %kopiert kode, b�r endres
        tk = zeros(C,1);
        c = floor((k-1)/Ntrain * C) + 1;
        tk(c) = 1;
        
        MSE1 = (gk-tk).*(gk).*(1-gk);     
        nablaW_MSE = nablaW_MSE + MSE1.*xk;
        %MSE = MSE + 0.5*(gk-tk)'*(gk-tk);
        
    end
    W = W - alpha.*nablaW_MSE;
    %MSEs(m) = MSE;
    nablaW_MSEs(m) = norm(nablaW_MSE); %riktig � gj�re det s�nn?
end
%}
%N� har vi trent classifier v�r, n� m� vi klassifisere testene v�re ??
%oppgave 1 c) 
%confusion matrix and the error rate for both the training set and the test
%set
%W er den trente classsifier vi m� benytte for � sortere test dataen

%testing the classifier:
all_test_t = [repmat([1 0 0], (50-Ntrain), 1); repmat([0 1 0], (50-Ntrain), 1); repmat([0 0 1], (50-Ntrain), 1)]; 

x_vec = [x_test, ones(C*Ntest, 1)];  
disp(W);
disp(x_vec); 
z= W.*x_vec; 
g= sigmoid(z); 

testset_class = zeros(1,60); %totalt 60 test set
testset_est= zeros(1,60); 

for t=1:60
    for c=1:3
        if all_test_t(c,t) == max(all_test_t(:,t))
            testset_class(t)= c; 
        end
        if g(c,t) == max(g(:,t))
            testset_est(t) = c; 
        end
    end
end


% 
% conf_matrix= zeros(C,C); % trenger en tabell som er 3x3, en med sann klasse, en med plassering. 
% %fill inn confusion matrix: 
% 
% %for x=1:Ntest 
    


%Plots 
plot(nablaW_MSEs); 
title('MSE gradient'); 
ylabel('Magnitude');
xlabel('iterations');