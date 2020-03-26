%Iris forsøk to 

D= 4; %dimensjon på inputvektorer, kan ses i filen til class1, class2 osv.
C= 3; %antall klasser
data_size_r= 50; %antall rader i klassefilene
data_size = 50 * C; %Lager en variabel med totalt antall rader
N_train = 30; %antall vi skal benytte i trening
N_test= 20;  % antall vi benytter testing
train_tot= N_train * C; %totale antallet vi skal benytte i trening sum av c1,c2,c3
test_tot = N_test * C; %totale antallet vi skal benytte i testing, sum av c1, c2, c3

M=30; %Antall iterasjoner vi skal benytte i MSE
alpha = 0.005; %step factor, vi må teste oss frem med denne, til vi er fornøyd

w0= 

%%Load data

x1 = load('class_1', '-ascii'); %Laster inn dataene fra de ulike klassene
x2 = load('class_2', '-ascii'); 
x3 = load ('class_3', '-ascii');

x1_train = x1(1:N_train, :);  % legger radene som skal benyttes i treningen som egen matrise
x2_train = x2(1:N_train, :); 
x3_train = x3(1:N_train, :);
x_train = [x1_train, x2_train, x3_train]; %egen matrise for all train data for alle klassene

x1_test = x1(N_train+1:end, :); % tar de 20 siste radene fra klassene som skal brukes i test
x2_test = x2(N_train+1:end, :);
x3_test = x3(N_train+1:end, :);
x_test = [x1_test, x2_test, x3_test]; %egen matrise for all test data



%%trening

%Vi ønsker å iterere M ganger for å finne MSE. 

for m= 1:M
    MSE =0; 
    grad=0; %opretter ny for hver gang
    for k= 1:Ntrain
        c = floor((k-1)/Ntrain * C) + 1;
        t = zeros(C, 1);
        t(c) = 1;  %lager en matrise over hvilke klasse testen vår tilhører
        %disp(t);
        
        xk=[x_train(k, :)'; 1];  %Verdiene i treningsmatrisen vår, tar ut en rad av gangen, itererer med k, altså 30 ganger
        disp(xk);
        z= W*xk +w0; 
        g=sigmoid(z); 
        grad= grad +((g-t)*g*(1-g))*xk; 
        MSE= MSE + 0.5*(g-t)^2; 
    end
end





