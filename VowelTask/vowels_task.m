clear;
N = 139;
Ntrain = 70;
Ntest = 69;
Nclass = 12;
Nfeatures = 3;
[files,dur,F0s,F1s,F2s,F3s,F4s,F120,F220,F320,F150,F250,F350,F180,F280,F380] = textread('vowdata_nohead.dat','%s%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f'); 
nfiles = 1668;
Prob_w = 1/Nclass;
% character 1:     m=man, w=woman, b=boy, g=girl 
% characters 2-3:  talker number 
% characters 4-5:  vowel (ae=�had�, ah=�hod�, aw=�hawed�, eh=�head�,  
%    er=�heard�, ei=�haid�, ih=�hid�, iy=�heed�, oa=/o/ as in �boat�,  
%    oo=�hood�, uh=�hud�, uw=�who�d�) 

%data processing
vowel = str2mat('ae','ah','aw','eh','er','ei','ih','iy','oa','oo','uh','uw'); 
talker_group = char('m','w','b','g');  

vowel_code = zeros(1,nfiles);
talker_number = zeros(1,nfiles);

files=char(files); % convert cell array to character matrix [nfiles,nchar]=size(filenames); 
for i=1:nfiles  
    vowel_code(i) = strmatch(files(i,4:5),vowel);  
    talker_group_code(i) = strmatch(files(i,1),talker_group);  
    talker_number(i) = str2num(files(i,2:3)); 
end; 




%Task 1 a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample mean for each klasse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mean_trainF1, mean_testF1] = find_mean(F1s,Nclass,N,Ntrain);
[mean_trainF2, mean_testF2] = find_mean(F2s,Nclass,N,Ntrain);
[mean_trainF3, mean_testF3] = find_mean(F3s,Nclass,N,Ntrain);


%12x3 matrix which consist of the means of f1, f2 og f3 for each
%class
means_train = [mean_trainF1 mean_trainF2 mean_trainF3];
means_test = [mean_testF1 mean_testF2 mean_testF3];

%make (12*70)x3 matrise for training
Fs = make_string(F1s,F2s,F3s,N,Ntrain,Nclass);
Ftot = [F1s F2s F3s];

% (12*69)x3 matrix for testing
test_vals = make_test_matrix(Ftot,N,Ntest,Nclass,Nfeatures);
test456 = cov(Fs(281:350,:));

%making histograms for each class
%{
for i =1:Nclass
%x = Ftot((i-1)*N+1:i*N,:);    %the total set
x = Fs((i-1)*Ntrain+1:i*Ntrain,:);  % the train set
%x = test_vals((i-1)*Ntest+1:i*Ntest,:);    % the test set
xn = i;
figure(1);   
subplot(3,4,i);   
hist(x,30);  % use 30 �bins�   
set(gca,'XLim',[1400 4700]);  % set x-axis limits between 50 & 500 Hz  
set(gca,'YLim',[0 25]);
%sgtitle('Histogram of F1 for each class','t');
title("Class " + i); 
end
%}

%histograms of F1 for each talker group
%{
%men
x = F1s(find(vowel_code==1));
figure(1);   
subplot(2,2,1);   
hist(x,20);  % use 20 �bins�   
set(gca,'XLim',[200 1400]);  % set x-axis limits between 50 & 500 Hz   
%xlim([50 500]);
title('adult males')

%women
x = F1s(find(talker_group_code==2));
figure(1);
subplot(2,2,2);
hist(x,20);
set (gca, 'Xlim',[200 1400]);
%xlim([50 500]);
title('adult females');

%boys
x = F1s(find(talker_group_code==3));
figure(1);
subplot(2,2,3);
hist(x,20);
set (gca, 'Xlim',[200 1400]);
%xlim([50 500]);
title('boys');

%girls
x = F1s(find(talker_group_code==4));
figure(1);
subplot(2,2,4);
hist(x,20);
set (gca, 'Xlim',[200 1400]);
%xlim([50 500]);
title('girls');
%}

% (12*3)x3 matrix consisting of the 12 covariance matrices
cov_matrices = zeros(Nclass*Nfeatures,Nfeatures);
for i = 1:Nclass
    cov_matrices((i-1)*3+1:(i*3),:) = find_cov(Fs, i, Ntrain);
end

%make covariance matrix from eq. (17) in the compendium
cov_mat2 = zeros(Nclass*Nfeatures,Nfeatures);
for c = 1:Nclass
   sigma2 = 0;
   for k = 1:Ntrain
       x = Fs((c-1)*Ntrain+k,:);
       mu = means_train(c,:);
       sigma2 = sigma2 + (x-mu)'*(x-mu);
   end
   cov_mat2((c-1)*Nfeatures+1:c*Nfeatures,:) = sigma2/Ntrain;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1b: designe gaussian classifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (69*12)x12 matrix which contains all discriminant functions for each
% sample, tested for each class
g_all_test = zeros(Ntest*Nclass, Nclass);
for i = 1:Ntest*Nclass
    for c = 1:Nclass
        x = test_vals(i,:);
        mu = means_train(c,:);
        cov_mat = cov_matrices((c-1)*Nfeatures+1:c*Nfeatures,:);
        g_all_test(i,c) = discriminant2(cov_mat, mu, x, Nfeatures, Prob_w);
    end
end

%testing the classifier
true_val_test = fill_in_truevalues(Nclass,Ntest);
testset = test_classifier(Nclass, Ntest,g_all_test,true_val_test);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make confusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conf_mat = compute_confusion(Nclass,Ntest, true_val_test, testset);
disp('Confusion matrix, full covariance pdfs');
disp(conf_mat);
%finds error rate
error_test = compute_error(Nclass,Ntest,conf_mat);
disp('error_rate, full covariance pdfs');
disp(error_test);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1c) diagonal covariance matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make (12*3)x3 matrix consisting of the diagonal cov matrix for each class
cov_diags = zeros(Nclass*Nfeatures,Nfeatures);
for i = 1:Nclass
    cov_diag = find_cov(Fs, i, Ntrain);
    cov_diags((i-1)*3+1:(i*3),:) = make_diag_cov(cov_diag, Nfeatures, Nfeatures);
end 

%fill in all discriminant functions in one matrix
g_all_diag = zeros(Ntest*Nclass, Nclass);
for i = 1:Ntest*Nclass
    for c = 1:Nclass
        x = test_vals(i,:);
        mu = means_train(c,:);
        cov_mat = cov_diags((c-1)*Nfeatures+1:c*Nfeatures,:);
        g_all_diag(i,c) = discriminant2(cov_mat, mu, x, Nfeatures, Prob_w);
    end
end

%testing the classifier
true_val_diag = fill_in_truevalues(Nclass,Ntest);
testset_diag = test_classifier(Nclass, Ntest,g_all_diag,true_val_diag);

%make confusion matrix
conf_mat_diag = compute_confusion(Nclass,Ntest, true_val_diag, testset_diag);
disp('Confusion matrix, diagonal covariance pdfs');
disp(conf_mat_diag);
%finds error rate
error_diag = compute_error(Nclass,Ntest,conf_mat_diag);
disp('Error rate, diagonal covariance pdfs');
disp(error_diag);


%%%%%%%%%%%%%%%%%%%%%
% Oppgave 2 test
%%%%%%%%%%%%%%%%%%%%%%

% Oppgave 2a) 
%make gmm distributions with diagonal matrices
%{
gmms2_diag = cell(Nclass,1);
for c = 1:Nclass
    data = Fs((c-1)*Ntrain+1:c*Ntrain,:);
    try
    gmms2_diag{c,1} = fitgmdist(data, 2,'RegularizationValue',0.1,'CovarianceType','diagonal'); %legger til regulariztionvalue for � unng� feilmelding
    catch exeption
        disp ('Noe er feil') 
        error=exeption.message
    end  
end

gmms3_diag = cell(Nclass,1);  
for c = 1:Nclass
    data = Fs((c-1)*Ntrain+1:c*Ntrain,:);
    try
    gmms3_diag{c,1} = fitgmdist(data, 3,'RegularizationValue',0.1,'CovarianceType','diagonal'); %legger til regulariztionvalue for � unng� feilmelding
    catch exeption
        disp ('Noe er feil') 
        error=exeption.message
    end
end
%}

%%%%%% 2b %%%%%%%%%%
% modify the classifier from 1b)

% make GMM models for M=2 and M = 3
%GMM with M = 2
gmms2 = cell(Nclass,1);
for c = 1:Nclass
    data = Fs((c-1)*Ntrain+1:c*Ntrain,:);
    %gmdistribution.fit(Fs,2);
    try
    gmms2{c,1} = fitgmdist(data, 2,'RegularizationValue',0.5); %legger til regularizationvalue for � unng� feilmelding
    catch exeption
        disp ('Noe er feil') 
        error=exeption.message
    end
    
end

%GMM with M=3
gmms3 = cell(Nclass,1);   
for c = 1:Nclass
    data = Fs((c-1)*Ntrain+1:c*Ntrain,:);
    try
    gmms3{c,1} = fitgmdist(data, 3,'RegularizationValue',0.5); %legger til regulariztionvalue for � unng� feilmelding
    catch exeption
        disp ('Noe er feil') 
        error=exeption.message
    end
end


%sort the data
cov_mats2 = zeros(Nclass*Nfeatures,2*Nfeatures);
cov_mats3 = zeros(Nclass*Nfeatures,3*Nfeatures);
means2 = zeros(Nclass,2*Nfeatures);
means3 = zeros(Nclass,3*Nfeatures);


%get means and cov matrices for each class, put in matrix
for i = 1:2
    for c=1:Nclass 
    sigma =gmms2{c,1}.Sigma(:,:,i);
    cov_mats2((c-1)*Nfeatures+1:c*Nfeatures,(i-1)*Nfeatures+1:i*Nfeatures)= sigma;
    means2(c,(i-1)*Nfeatures+1:i*Nfeatures) = gmms2{c,1}.mu(i,:);
    end
end


for i = 1:3
    for c=1:Nclass 
    sigma =gmms3{c,1}.Sigma(:,:,i);
    cov_mats3((c-1)*Nfeatures+1:c*Nfeatures,(i-1)*Nfeatures+1:i*Nfeatures)= sigma;
    means3(c,(i-1)*Nfeatures+1:i*Nfeatures) = gmms3{c,1}.mu(i,:);
    end
end

%matrix that will be filled with the discriminant functions
pdf1 = zeros(Ntest*Nclass,Nclass);

%for GMM with M=2
for k = 1:Ntest*Nclass
    x = test_vals(k,:);
    for c = 1:Nclass  
        pdf = 0;
        for i=1:2             
            cov = cov_mats2((c-1)*Nfeatures+1:c*Nfeatures,(i-1)*Nfeatures+1:i*Nfeatures);
            mu = means2(c,(i-1)*Nfeatures+1:i*Nfeatures);
            weight = gmms2{c, 1}.ComponentProportion(1,i);  
            frac = (sqrt((2*pi)^(Nfeatures)*det(cov)))^(-1);
            expo = exp((-0.5)*(x-mu)*cov^(-1)*(x-mu)');
            pdf = pdf + frac*expo*weight;
        end
    pdf1(k,c) = pdf;
    end
end

%for GMM with M=3
pdf3 = zeros(Ntest*Nclass,Nclass);
for k = 1:Ntest*Nclass
    x = test_vals(k,:);
    for c = 1:Nclass
        pdf = 0;
        for i=1:3   
            cov = cov_mats3((c-1)*Nfeatures+1:c*Nfeatures,(i-1)*Nfeatures+1:i*Nfeatures);
            mu = means3(c,(i-1)*Nfeatures+1:i*Nfeatures);
            weight = gmms3{c, 1}.ComponentProportion(1,i);   
            frac = (sqrt((2*pi)^(Nfeatures)*det(cov)))^(-1);
            expo = exp((-0.5)*(x-mu)*cov^(-1)*(x-mu)');
            pdf = pdf + frac*expo*weight;
        end
    pdf3(k,c) = pdf;
    end
end

% testing the classifiers
testset_gmm = test_classifier(Nclass, Ntest,pdf1,true_val_diag);
testset_gmm3 = test_classifier(Nclass, Ntest,pdf3,true_val_diag);

%make confusion matrix
%M=2
conf_mat_gmm = compute_confusion(Nclass,Ntest, true_val_diag, testset_gmm);
disp('conf mat, 2 mixtures');
disp(conf_mat_gmm);
%M=3
conf_mat_gmm3 = compute_confusion(Nclass,Ntest, true_val_diag, testset_gmm3);
disp('conf mat, 3 mixtures');
disp(conf_mat_gmm3);

%find error rate
%M=2
error_gmm = compute_error(Nclass,Ntest,conf_mat_gmm);
disp('error rate, 2 mixtures');
disp(error_gmm);
%M=3
error_gmm3 = compute_error(Nclass,Ntest,conf_mat_gmm3);
disp('error rate, 3 mixtures');
disp(error_gmm3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%make covariance matrix from data set
function cov_matrix = find_cov(string, class_num, N)
    x_string = string((class_num-1)*N+1:class_num*N,:);
    cov_matrix = cov(x_string);
end

%function 3.4 from compendium
function g_i = discriminant2(cov_matrix, mu, x,Nfeatures, prior)
    frac = (sqrt((2*pi)^(Nfeatures)*det(cov_matrix)))^(-1);
    expo = exp(-0.5*(x-mu)*cov_matrix^(-1)*(x-mu)');
    g_i = frac*expo*prior;
end

%Nclassx1 matrix consisting of the mean of one feature for each klasse
function [mean_train, mean_test] = find_mean(string,Nclass,Ntot,Ntrain)
    mean_train = zeros(Nclass,1);
    mean_test = zeros(Nclass,1);

    for i = 1:Nclass    
        y = string((i*Ntot-Ntot)+1:i*Ntot);
        mean_train(i,1) = mean(y(1:Ntrain));
        mean_test(i,1) = mean(y(Ntrain+1:Ntot));
    end
end

%matrise best�ende av inputs fra 3 ulike features for � gj�re behandlingen
%av data lettere
function matrix = make_string(string1,string2,string3,Ntot,Ntrain,Nclass)
    matrix_tot = [string1 string2 string3];
    matrix = zeros(Ntrain*Nclass,3);
    for i = 1:Nclass
        x_string = matrix_tot((i-1)*Ntot+1:(i-1)*Ntot+Ntrain,:);
        matrix((i*Ntrain-Ntrain)+1:i*Ntrain,:) = x_string;
    end
end

%make matrix for test set
function test_inputs = make_test_matrix(inputs,N,Ntest,C,Nfeatures)
    test_inputs = zeros(C*Ntest,Nfeatures);
    for c=1:C
        x_string = inputs((c-1)*N+71:c*N,:);
        test_inputs((c-1)*Ntest+1:c*Ntest,:) = x_string;
    end
end

%fyller inn sanne verdier til bruk til confusion matrix
function truevalues = fill_in_truevalues(C,N)
truevalues=zeros(C*N,1);
    for c = 1:C
        for i = 1:N
            truevalues(i+(c-1)*N,1) = c;
        end
    end
end

%lager confusion matrix
function conf_mat = compute_confusion(C,N, set_class, set_est)  
    conf_mat=zeros(C);
    for t=1:N*C
        x=set_class(t); 
        y=set_est(t); 
        conf_mat(x,y)=conf_mat(x,y)+1; 
    end
end

%finner error rate 
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

%tester klassifisereren, fyller inn resultatene fra discriminant functions
function testset = test_classifier(C,N,g_all,true_set)
testset = zeros(N*C,1);
for x = 1:N*C
    for c = 1:C
        if g_all(x,c) == max(g_all(x,:))
            testset(x,1)= c;
        end
    end
end
end

function cov_diag = make_diag_cov(cov_matrix,width,length)
    for i = 1:width
        for j = 1:length
            if j ~= i
                cov_matrix(i,j) = 0;
            end
        end
    end
    cov_diag = cov_matrix;
end