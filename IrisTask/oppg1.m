%Class 1 = Iris Setosa
%Class 2 = Iris Versicolour
%Class 3 = Iris Virginica


clear;
Ntrain = 30;
Ntest = 20;
C = 3; % # classes
D = 4; % feature dimension
M = 1000; % # iterations (velger vi denne st�rrelsen selv?)
alpha = 0.008; %teste seg fram med denne
W0 = zeros(C,D); 
w0 = zeros(C,1); 
W = [W0 w0]; %forenkler W til senere uttrykk
nablaW_MSEs = zeros(1,M);
MSEs= zeros(1,M); %creating space for the matrix
g_all = zeros(Ntrain*C,C); %er denne n�dvendig?

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
% %function [W, MSEs, nablaW_MSEs] = train_classifier(M, N, C, x_vec, Weight)
[W,MSEs, g_all]=train_classifier(M, Ntrain, C, x_train, W, alpha);  

trainset_class = fill_in_truevalues(C,Ntrain);
trainset_est=zeros(1,90);


 for x=1:Ntrain*C 
    for c=1:C
        if g_all(x,c) == max(g_all(x,:))
            trainset_est(x)= c;
        end
    end
 end

testset_class = fill_in_truevalues(C,Ntest);
testset_est = zeros(1,60);
g_all_test = zeros(C*Ntest,C);

for k=1:C*Ntest
    [gk_test, x] = g_value(k, x_test, W); 
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
conf_matrix_train = compute_confusion(C, Ntrain, trainset_class, trainset_est);
disp('Confusion matrix, train set, four features');
disp(conf_matrix_train);

%confusion matrix for test set
conf_matrix_test = compute_confusion(C, Ntest, testset_class, testset_est);
disp('Confusion matrix, test set, four features');
disp(conf_matrix_test);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finner error rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_train = compute_error(C,Ntrain,conf_matrix_train);
disp('Error rate for train set:');
disp(error_train);
error_test = compute_error(C,Ntest,conf_matrix_test);
disp('Error rate for test set:');
disp(error_test);


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



alpha_2 = 0.008; %teste seg fram med denne
W0_2 = zeros(C,D); 
w0_2 = zeros(C,1); 
W_2 = [W0_2 w0_2]; %forenkler W til senere uttrykk
nablaW_MSEs_2 = zeros(1,M);
MSEs_2= zeros(1,M); %creating space for the matrix
g_all_2 = zeros(Ntrain*C,C); %er denne n�dvendig?

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
[W_2,MSEs_2, g_all_2]=train_classifier(M, Ntrain, C, x2_train, W_2, alpha_2); 

trainset_class_2 = fill_in_truevalues(C,Ntrain);
trainset_est_2=zeros(1,90);

%classifying the inputs
for x=1:Ntrain*C 
    for c=1:C
        if g_all_2(x,c) == max(g_all_2(x,:))
            trainset_est_2(x)= c;
        end
    end
end

testset_class_2 = fill_in_truevalues(C,Ntest);
testset_est_2 = zeros(1,60);
g_all_test_2 = zeros(C*Ntest,C);

for k=1:C*Ntest
    gk_test_2= g_value(k, x_test_2, W_2); 
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
conf_matrix_train_2= compute_confusion(C, Ntrain, trainset_class_2, trainset_est_2); 
disp('Confusion matrix, train set, four features, 30 last samples');
disp(conf_matrix_train_2);

%confusion matrix for test set
conf_matrix_test_2 = compute_confusion(C, Ntest, testset_class_2, testset_est_2); 
disp('Confusion matrix, test set, four features, 20 first samples');
disp(conf_matrix_test_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finner error rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_train_2 = compute_error(C,Ntrain,conf_matrix_train_2);
disp('Error rate for train set:');
disp(error_train_2);
error_test_2 = compute_error(C,Ntest,conf_matrix_test_2);
disp('Error rate for test set:');
disp(error_test_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oppgave 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
x = F0s(find(talker_group_code==1));   
figure(1);   
subplot(2,2,1);   
hist(x,20);  % use 20 �bins�   
set(gca,'XLim',[50 500]);  % set x-axis limits between 50 & 500 Hz   
%xlim([50 500]);
title('adult males')
%}

%%{
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
%%}

%the sepal width-feature has the biggest overlap between the classes - we
%remove this one

W0_3 = zeros(C);
W_3 = [W0_3 w0];
nablaW_MSEs_3 = zeros(1,M);
MSEs_3= zeros(1,M); %creating space for the matrix
g_all_3 = zeros(Ntrain*C,C); %er denne n�dvendig?

%train sets. Removes the second feature
x_train_3 = x_train; 
x_train_3(:,2)=[]; 


% test sets
x_test_3=x_test; 
x_test_3(:,2)=[]; 


% %training the classifier
% function [W, MSEs, g_all] = train_classifier(M, N, C, x_vec, W, alpha) 
[W_3, MSEs_3, g_all_3]=train_classifier(M, Ntrain, C, x_train_3, W_3, alpha); 

trainset_class_3=fill_in_truevalues(C,Ntrain); 


trainset_est_3=zeros(1,90);


 for x=1:Ntrain*C 
    for c=1:C
        if g_all_3(x,c) == max(g_all_3(x,:))
            trainset_est_3(x)= c;
        end
    end
 end
testset_class_3=fill_in_truevalues(C, Ntest); 

testset_est_3 = zeros(1,60);
g_all_test_3 = zeros(C*Ntest,C);

for k=1:C*Ntest
    gk_test_3= g_value(k, x_test_3, W_3); 
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
conf_matrix_train_3= compute_confusion(C, Ntrain, trainset_class_3, trainset_est_3); 
disp('Confusion matrix, train set, three features');
disp(conf_matrix_train_3);

%confusion matrix for test set
conf_matrix_test_3 = compute_confusion(C, Ntest, testset_class_3, testset_est_3); 
disp('Confusion matrix, test set, three features');
disp(conf_matrix_test_3);

%error rates
error_train_3 = compute_error(C,Ntrain,conf_matrix_train_3);
disp('Error rate for train set:');
disp(error_train_3);
error_test_3 = compute_error(C,Ntest,conf_matrix_test_3);
disp('Error rate for test set:');
disp(error_test_3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Oppgave 2b) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W0_4 = zeros(C,2);
W_4 = [W0_4 w0];
nablaW_MSEs_4 = zeros(1,M);
MSEs_4= zeros(1,M); %creating space for the matrix
g_all_4 = zeros(Ntrain*C,C);

%remove a feautre from trainset and testset
x_train_4 = x_train_3; 
x_train_4(:,3)=[]; 


x_test_4=x_test_3; 
x_test_4(:,3)=[]; 

%train classifier
[W_4, MSEs_4, g_all_4]=train_classifier(M, Ntrain, C, x_train_4, W_4, alpha); 

% fill inn true values
trainset_class_4=fill_in_truevalues(C,Ntrain); 

%estimation
trainset_est_4=zeros(1,90);


 for x=1:Ntrain*C 
    for c=1:C
        if g_all_4(x,c) == max(g_all_4(x,:))
            trainset_est_4(x)= c;
        end
    end
 end
 %fill in true values testset
 testset_class_4=fill_in_truevalues(C, Ntest); 

 %testset estimation
testset_est_4 = zeros(1,60);
g_all_test_4 = zeros(C*Ntest,C);

for k=1:C*Ntest
    %gk_test_4= g_value(k, x_test_4, W_4); 
    xk = x_test(k,:)';
    x = [xk' 1]'; %sp�rre solveig om denne, hvorfor transposer man to ganger
    zk = W*x; 
           
    gk_test_4 = sigmoid(zk);
    g_all_test_4(k,:) = gk_test_4';
    
end

for x = 1:C*Ntest
    for c = 1:C
        if g_all_test_4(x,c) == max(g_all_test_4(x,:))
            testset_est_4(x)= c;
        end
    end
end


conf_matrix_train_4= compute_confusion(C, Ntrain, trainset_class_4, trainset_est_4); 
disp('Confusion matrix, train set, two features');
disp(conf_matrix_train_4);

conf_matrix_test_4 = compute_confusion(C, Ntest, testset_class_4, testset_est_4); 
disp('Confusion matrix, test set, two features');
disp(conf_matrix_test_4);

%error rates
error_train_4 = compute_error(C,Ntrain,conf_matrix_train_4);
disp('Error rate for train set:');
disp(error_train_4);
error_test_4 = compute_error(C,Ntest,conf_matrix_test_4);
disp('Error rate for test set:');
disp(error_test_4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  OPPGAVE 2b - 1 feature
W0_5 = zeros(C,1);
W_5 = [W0_5 w0];
nablaW_MSEs_5 = zeros(1,M);
MSEs_5= zeros(1,M); %creating space for the matrix
g_all_5 = zeros(Ntrain*C,C);

%remove a feautre from trainset and testset
x_train_5 = x_train_4; 
x_train_5(:,2)=[]; 


x_test_5=x_test_4; 
x_test_5(:,2)=[]; 

%train classifier
[W_5, MSEs_5, g_all_5]=train_classifier(M, Ntrain, C, x_train_5, W_5, alpha); 

% fill inn true values
trainset_class_5=fill_in_truevalues(C,Ntrain); 

%estimation
trainset_est_5=zeros(1,90);


 for x=1:Ntrain*C 
    for c=1:C
        if g_all_5(x,c) == max(g_all_5(x,:))
            trainset_est_5(x)= c;
        end
    end
 end
 %fill in true values testset
 testset_class_5=fill_in_truevalues(C, Ntest); 

 %testset estimation
testset_est_5 = zeros(1,60);
g_all_test_5 = zeros(C*Ntest,C);

for k=1:C*Ntest
    gk_test_5= g_value(k, x_test_5, W_5); 
%     xk = x_test(k,:)';
%     x = [xk' 1]'; %sp�rre solveig om denne, hvorfor transposer man to ganger
%     zk = W_5*x; 
%            
%     gk_test_5 = sigmoid(zk);
    g_all_test_5(k,:) = gk_test_5';
    
end

for x = 1:C*Ntest
    for c = 1:C
        if g_all_test_5(x,c) == max(g_all_test_5(x,:))
            testset_est_5(x)= c;
        end
    end
end


conf_matrix_train_5= compute_confusion(C, Ntrain, trainset_class_5, trainset_est_5); 
disp('Confusion matrix, train set, one features');
disp(conf_matrix_train_5);


conf_matrix_test_5 = compute_confusion(C, Ntest, testset_class_5, testset_est_5); 
disp('Confusion matrix, test set, one features');
disp(conf_matrix_test_5);

%error rates
error_train_5 = compute_error(C,Ntrain,conf_matrix_train_5);
disp('Error rate for train set:');
disp(error_train_5);
error_test_5 = compute_error(C,Ntest,conf_matrix_test_5);
disp('Error rate for test set:');
disp(error_test_5);

%funskjon for � trene classifieren 
function [W_0, MSEs_0, g_all_0] = train_classifier(M, N, C, x_vec, W_0, alpha) 
for m=1:M
    nablaW_MSE=0; 
    MSE=0;
    for k=1:N*C
        [gk,x]=g_value(k, x_vec, W_0); 
        g_all_0(k,:) = gk';
        tk=zeros(C,1);
        c = floor((k-1)/(N*C)*C) + 1;
        tk(c)=1; 
        MSE1 = (gk-tk).*(gk).*(1-gk);     
        nablaW_MSE = nablaW_MSE + MSE1*x';
        MSE = MSE + 0.5*(gk-tk)'*(gk-tk);
    end
    W_0 = W_0 - alpha.*nablaW_MSE;
    MSEs_0(m) = MSE; %brukes til � tune alpha til riktig verdi (konvergering)
    %nablaW_MSEs(m) = norm(nablaW_MSE); %riktig � gj�re det s�nn?
end
end
        
%plasserers i en loop der i er den itererende tingen. 
function [g,x] = g_value(i, x_vec, Weight) %x_vec er vectoren som inneholder alle klassedatene
 xk=x_vec(i,:)'; 
 x=[xk' 1]'; 
 zk=Weight*x; 
 g=sigmoid(zk);
end

function conf_mat = compute_confusion(C,N, set_class, set_est)  
    conf_mat=zeros(C);
    for t=1:N*C
        x=set_class(t); 
        y=set_est(t); 
        conf_mat(x,y)=conf_mat(x,y)+1; 
    end
end

function truevalues = fill_in_truevalues(C,N)
truevalues=zeros(1,C*N);
truevalues(1, 1:N)=1; 
truevalues(1, (N+1):(N*2))=2; 
truevalues(1, (N*2+1):(N*C))=3; 
end

function error_rate = compute_error(C,N,conf_mat)
sum_error=0;
for i = 1:C
    for j = 1:C
        if j ~= i           
            sum_error = sum_error +conf_mat(j,i);
        end
    end
end
error_rate = sum_error/(C*N);
end

%%heihei
