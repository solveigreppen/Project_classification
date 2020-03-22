%Class 1 = Iris Setosa
%Class 2 = Iris Versicolour
%Class 3 = Iris Virginica

clear;
Ntrain = 30;
Ntest = 20;
C = 3; % # classes
D = 4; % # columns in textfile
M = 1000; % # iterations (velger vi denne størrelsen selv?)
alpha = 0.05; %teste seg fram med denne
W = zeros(C,D);
w0 = zeros(C,1); %hva skal vi fylle i denne?
%W = [W0 w0]; skal denne legges inn? I så fall hvorfor?

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
    
    for k = 1:Ntrain
        xk = x_train(k,:).'; %ønskelig å transponere denne? Ja, virker sånn pga matrix dimension
        zk = W.*xk + w0;
        gk = (1+exp(-zk)).^-1; %bruke innebygd sigmoid eller lage egen funksjon?
    end
    %{
    nabla_gk_MSE = gk-Ntrain;
    nabla_zk_gk = gk.*(1-gk);
    nabla_W_zk = xk.';

    nablaW_MSE = nablaW_MSE + nabla_gk_MSE.* nabla_zk_gk .* nabla_W_zk;
    W = W - aplha.*nablaW_MSE;
    %}
end


