%Class 1 = Iris Setosa
%Class 2 = Iris Versicolour
%Class 3 = Iris Virginica


clear;
Ntrain = 30;
Ntest = 20;
C = 3; % # classes
D = 4; % feature dimension
M = 1000; % # iterations (velger vi denne størrelsen selv?)
alpha = 0.008; %teste seg fram med denne
W0 = zeros(C,D); 
w0 = zeros(C,1); 
W = [W0 w0]; %forenkler W til senere uttrykk
nablaW_MSEs = zeros(1,M);
MSEs= zeros(1,M); %creating space for the matrix
g_all = zeros(Ntrain*C,C); %er denne nødvendig?

%load the datas
x1all = load('class_1','-ascii');
x2all = load('class_2','-ascii');
x3all = load('class_3','-ascii');

%train sets
x1_train = x1all(1:Ntrain,:);
x2_train = x2all(1:Ntrain,:);
x3_train = x3all(1:Ntrain,:);
x_train = [x1_train; x2_train; x3_train]; %90x4

% test sets
x1_test = x1all(Ntrain+1:end,:);
x2_test = x2all(Ntrain+1:end,:);
x3_test = x3all(Ntrain+1:end,:);
x_test = [x1_test; x2_test; x3_test]; %60x4


%training the classifier
for m = 1:M
    nablaW_MSE = 0;
    MSE=0; 
        
        for k = 1:Ntrain*C
            xk = x_train(k,:)';
            x = [xk' 1]';
            zk = W*x; %forenkle til zk = Wx?
           
            gk = sigmoid(zk); %bruke innebygd sigmoid eller lage egen funksjon for å forenkle koden?
            g_all(k,:) = gk';
            %kopiert kode, bør endres
            tk = zeros(C,1);
            c = floor((k-1)/(Ntrain*C)*C) + 1;
            tk(c) = 1;
       
            MSE1 = (gk-tk).*(gk).*(1-gk);     
            nablaW_MSE = nablaW_MSE + MSE1*x';
            MSE = MSE + 0.5*(gk-tk)'*(gk-tk);
        end
    %end
    W = W - alpha.*nablaW_MSE;
    MSEs(m) = MSE; %brukes til å tune alpha til riktig verdi (konvergering)
    nablaW_MSEs(m) = norm(nablaW_MSE); %riktig å gjøre det sånn?
end


trainset_class=zeros(1,90);
%fyller inn sanne verdier
trainset_class(1,1:Ntrain) = 1;
trainset_class(1,31:60) = 2;
trainset_class(1,61:90) = 3;

trainset_est=zeros(1,90);


 for x=1:Ntrain*C 
    for c=1:C
        %{
        if all_test_t(c,x) == max(all_test_t(:,x))
            trainset_class(x)=c;
        end
        %}
        if g_all(x,c) == max(g_all(x,:))
            trainset_est(x)= c;
        end
    end
 end

testset_class = zeros(1,60);
%fyller inn sanne verdier
testset_class(1,1:20) = 1;
testset_class(1,21:40) = 2;
testset_class(1,41:60) = 3;

testset_est = zeros(1,60);
g_all_test = zeros(C*Ntest,C);

for k=1:C*Ntest
    xk = x_test(k,:)';
    x = [xk' 1]';
    zk = W*x; 
           
    gk_test = sigmoid(zk); 
    g_all_test(k,:) = gk_test';
    
end

for x = 1:C*Ntest
    for c = 1:C
        if g_all_test(x,c) == max(g_all_test(x,:))
            testset_est(x)= c;
        end
    end
end

%confusion matrix for train set
conf_matrix= zeros(C);  
for t=1:Ntrain*C
    x=trainset_class(t); 
    y=trainset_est(t); 
    conf_matrix(x,y)= conf_matrix(x,y) +1;
end
disp('Confusion matrix, train set, four features');
disp(conf_matrix);

%confusion matrix for test set
conf_matrix_test= zeros(C);  
for t=1:Ntest*C
    x=testset_class(t); 
    y=testset_est(t); 
    conf_matrix_test(x,y)= conf_matrix_test(x,y) +1;
end
disp('Confusion matrix, test set, four features');
disp(conf_matrix_test);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finner error rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_error_train=0;
sum_error_test=0;
for i = 1:C
    for j = 1:C
        if j ~= i
            sum_error_train = sum_error_train +conf_matrix(j,i);
            sum_error_test = sum_error_test +conf_matrix_test(j,i);
        end
    end
end

error_rate_train = sum_error_train/Ntrain;
error_rate_test = sum_error_test/Ntest;


%Plots 
%{
plot(nablaW_MSEs); 
title('MSE gradient'); 
ylabel('Magnitude');
xlabel('iterations');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oppgave 1 d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
alpha_2 = 0.008; %teste seg fram med denne
W0_2 = zeros(C,D); 
w0_2 = zeros(C,1); 
W_2 = [W0_2 w0_2]; %forenkler W til senere uttrykk
nablaW_MSEs_2 = zeros(1,M);
MSEs_2= zeros(1,M); %creating space for the matrix
g_all_2 = zeros(Ntrain*C,C); %er denne nødvendig?

x4_train = x1all(Ntest+1:end,:);
x5_train = x2all(Ntest+1:end,:);
x6_train = x3all(Ntest+1:end,:);
x2_train = [x4_train; x5_train; x6_train]; %90x4

% test sets

x4_test = x1all(1:Ntest,:);
x5_test = x2all(1:Ntest,:);
x6_test = x3all(1:Ntest,:);
x_test_2 = [x4_test; x5_test; x6_test]; %60x4


%training the classifier
for m = 1:M
    nablaW_MSE_2 = 0;
    MSE_2=0; 
    
        for k = 1:Ntrain*C
            xk = x2_train(k,:)';
            x = [xk' 1]';
            zk = W_2*x; %forenkle til zk = Wx?
           
            gk = sigmoid(zk); %bruke innebygd sigmoid eller lage egen funksjon for å forenkle koden?
            g_all_2(k,:) = gk';
            %kopiert kode, bør endres
            tk = zeros(C,1);
            c = floor((k-1)/(Ntrain*C)*C) + 1;
            tk(c) = 1;
       
            MSE2 = (gk-tk).*(gk).*(1-gk);     
            nablaW_MSE_2 = nablaW_MSE_2 + MSE2*x';
            MSE_2 = MSE_2 + 0.5*(gk-tk)'*(gk-tk);
        end
    %end
    W_2 = W_2 - alpha_2.*nablaW_MSE_2;
    MSEs_2(m) = MSE_2; %brukes til å tune alpha til riktig verdi (konvergering)
    nablaW_MSEs_2(m) = norm(nablaW_MSE_2); %riktig å gjøre det sånn?
end

trainset_class_2=zeros(1,90);
%fyller inn sanne verdier
trainset_class_2(1,1:Ntrain) = 1;
trainset_class_2(1,31:60) = 2;
trainset_class_2(1,61:90) = 3;

trainset_est_2=zeros(1,90);


 for x=1:Ntrain*C 
    for c=1:C
        %{
        if all_test_t(c,x) == max(all_test_t(:,x))
            trainset_class(x)=c;
        end
        %}
        if g_all_2(x,c) == max(g_all_2(x,:))
            trainset_est_2(x)= c;
        end
    end
 end

testset_class_2 = zeros(1,60);
%fyller inn sanne verdier
testset_class_2(1,1:20) = 1;
testset_class_2(1,21:40) = 2;
testset_class_2(1,41:60) = 3;

testset_est_2 = zeros(1,60);
g_all_test_2 = zeros(C*Ntest,C);

for k=1:C*Ntest
    xk = x_test_2(k,:)';
    x = [xk' 1]';
    zk = W_2*x; 
           
    gk_test_2 = sigmoid(zk); 
    g_all_test_2(k,:) = gk_test_2';
    
end

for x = 1:C*Ntest
    for c = 1:C
        if g_all_test_2(x,c) == max(g_all_test_2(x,:))
            testset_est_2(x)= c;
        end
    end
end

%confusion matrix for train set
conf_matrix_2= zeros(C);  
for t=1:Ntrain*C
    x=trainset_class_2(t); 
    y=trainset_est_2(t); 
    conf_matrix_2(x,y)= conf_matrix_2(x,y) +1;
end
disp(conf_matrix_2);

%confusion matrix for test set
conf_matrix_test_2= zeros(C);  
for t=1:Ntest*C
    x=testset_class_2(t); 
    y=testset_est_2(t); 
    conf_matrix_test_2(x,y)= conf_matrix_test_2(x,y) +1;
end
disp(conf_matrix_test_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finner error rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_error_train_2=0;
sum_error_test_2=0;
for i = 1:C
    for j = 1:C
        if j ~= i
            sum_error_train_2 = sum_error_train_2 +conf_matrix_2(j,i);
            sum_error_test_2 = sum_error_test_2 +conf_matrix_test_2(j,i);
        end
    end
end

error_rate_train_2 = sum_error_train_2/Ntrain;
error_rate_test_2 = sum_error_test_2/Ntest;
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oppgave 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
x = F0s(find(talker_group_code==1));   
figure(1);   
subplot(2,2,1);   
hist(x,20);  % use 20 “bins”   
set(gca,'XLim',[50 500]);  % set x-axis limits between 50 & 500 Hz   
%xlim([50 500]);
title('adult males')
%}

%{
xh = x1_train(:,1);
figure(1);
subplot(3,4,1);
hist(xh,20)
set(gca,'XLim',[4 8]);  % set x-axis limits between 0-10  
title('Class 1, sepal length');

xh2 = x1_train(:,2);
figure(1);
subplot(3,4,2);
hist(xh2,20)
set(gca,'XLim',[1.5 4.5]);  % set x-axis limits between 0-10  
title('Class 1, sepal width');

xh3 = x1_train(:,3);
figure(1);
subplot(3,4,3);
hist(xh3,20)
set(gca,'XLim',[0.5 7.5]);  % set x-axis limits between 0-10  
title('Class 1, petal length');

xh4 = x1_train(:,4);
figure(1);
subplot(3,4,4);
hist(xh4,20)
set(gca,'XLim',[0 3]);  % set x-axis limits between 0-10  
title('Class 1, petal width');

xh5 = x2_train(:,1);
figure(1);
subplot(3,4,5);
hist(xh5,20)
set(gca,'XLim',[4 8]);  % set x-axis limits between 0-10  
title('Class 2, sepal length');

xh6 = x2_train(:,2);
figure(1);
subplot(3,4,6);
hist(xh6,20)
set(gca,'XLim',[1.5 4.5]);  % set x-axis limits between 0-10  
title('Class 2, sepal width');

xh7 = x2_train(:,3);
figure(1);
subplot(3,4,7);
hist(xh7,20)
set(gca,'XLim',[0.5 7.5]);  % set x-axis limits between 0-10  
title('Class 2, petal length');

xh8 = x2_train(:,4);
figure(1);
subplot(3,4,8);
hist(xh8,20)
set(gca,'XLim',[0 3]);  % set x-axis limits between 0-10  
title('Class 2, petal width');

xh9 = x3_train(:,1);
figure(1);
subplot(3,4,9);
hist(xh9,20)
set(gca,'XLim',[4 8]);  % set x-axis limits between 0-10  
title('Class 3, sepal length');

xh10 = x3_train(:,2);
figure(1);
subplot(3,4,10);
hist(xh10,20)
set(gca,'XLim',[1.5 4.5]);  % set x-axis limits between 0-10  
title('Class 3, sepal width');

xh11 = x3_train(:,3);
figure(1);
subplot(3,4,11);
hist(xh11,20)
set(gca,'XLim',[0.5 7.5]);  % set x-axis limits between 0-10  
title('Class 3, petal length');

xh12 = x3_train(:,4);
figure(1);
subplot(3,4,12);
hist(xh12,20)
set(gca,'XLim',[0 3]);  % set x-axis limits between 0-10  
title('Class 3, petal width');
%}

%the sepal width-feature has the biggest overlap between the classes - we
%remove this one

W0_3 = zeros(C);
W_3 = [W0_3 w0];
nablaW_MSEs_3 = zeros(1,M);
MSEs_3= zeros(1,M); %creating space for the matrix
g_all_3 = zeros(Ntrain*C,C); %er denne nødvendig?

%train sets. Removes the second feature
x1_train_3 = x1_train;
x1_train_3(:,2) = [];

x2_train_3 = x2_train;
x2_train_3(:,2) = [];

x3_train_3 = x3_train;
x3_train_3(:,2) = [];

x_train_3 = [x1_train_3; x2_train_3; x3_train_3]; %90x4

% test sets
x1_test_3 = x1_test;
x1_test_3(:,2) = [];

x2_test_3 = x2_test;
x2_test_3(:,2) = [];

x3_test_3 = x3_test;
x3_test_3(:,2) = [];

x_test_3 = [x1_test_3; x2_test_3; x3_test_3]; %90x4

%training the classifier
for m = 1:M
    nablaW_MSE_3 = 0;
    MSE_3=0; 
        
        for k = 1:Ntrain*C
            xk = x_train_3(k,:)';
            x = [xk' 1]';
            zk = W_3*x; %forenkle til zk = Wx?
           
            gk = sigmoid(zk); %bruke innebygd sigmoid eller lage egen funksjon for å forenkle koden?
            g_all_3(k,:) = gk';
            %kopiert kode, bør endres
            tk = zeros(C,1);
            c = floor((k-1)/(Ntrain*C)*C) + 1;
            tk(c) = 1;
       
            MSE1 = (gk-tk).*(gk).*(1-gk);     
            nablaW_MSE_3 = nablaW_MSE_3 + MSE1*x';
            MSE_3 = MSE + 0.5*(gk-tk)'*(gk-tk);
        end
    %end
    W_3 = W_3 - alpha.*nablaW_MSE_3;
    MSEs_3(m) = MSE; %brukes til å tune alpha til riktig verdi (konvergering)
    nablaW_MSEs_3(m) = norm(nablaW_MSE); %riktig å gjøre det sånn?
end

trainset_class_3=zeros(1,90);
%fyller inn sanne verdier
trainset_class_3(1,1:Ntrain) = 1;
trainset_class_3(1,31:60) = 2;
trainset_class_3(1,61:90) = 3;

trainset_est_3=zeros(1,90);


 for x=1:Ntrain*C 
    for c=1:C
        %{
        if all_test_t(c,x) == max(all_test_t(:,x))
            trainset_class(x)=c;
        end
        %}
        if g_all_3(x,c) == max(g_all_3(x,:))
            trainset_est_3(x)= c;
        end
    end
 end

testset_class_3 = zeros(1,60);
%fyller inn sanne verdier
testset_class_3(1,1:20) = 1;
testset_class_3(1,21:40) = 2;
testset_class_3(1,41:60) = 3;

testset_est_3 = zeros(1,60);
g_all_test_3 = zeros(C*Ntest,C);

for k=1:C*Ntest
    xk = x_test_3(k,:)';
    x = [xk' 1]';
    zk = W_3*x; 
           
    gk_test_3 = sigmoid(zk); 
    g_all_test_3(k,:) = gk_test_3';
    
end

for x = 1:C*Ntest
    for c = 1:C
        if g_all_test_3(x,c) == max(g_all_test_3(x,:))
            testset_est_3(x)= c;
        end
    end
end

%confusion matrix for train set
conf_matrix_3= zeros(C);  
for t=1:Ntrain*C
    x=trainset_class_3(t); 
    y=trainset_est_3(t); 
    conf_matrix_3(x,y)= conf_matrix_3(x,y) +1;
end
disp('Confusion matrix, train set, three features');
disp(conf_matrix_3);

%confusion matrix for test set
conf_matrix_test_3= zeros(C);  
for t=1:Ntest*C
    x=testset_class_3(t); 
    y=testset_est_3(t); 
    conf_matrix_test_3(x,y)= conf_matrix_test_3(x,y) +1;
end
disp('Confusion matrix, test set, three features');
disp(conf_matrix_test_3);
