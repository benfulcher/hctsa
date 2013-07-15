echo on;
clc;

% load and plot data on geyser eruption durations and waiting times

load geyser;                          
                                      
clf; 
subplot(211);
plot(geyser(:,1),geyser(:,2),'o');  
xlabel('dur'); ylabel('wait time');
axis('square'); hold on;

% Hit any key to continue 

pause;
clc;

% divide up dataset

Xtrain=geyser(1:200,:); Xtest=geyser(201:295,:);

% train HMM with 3 states for 30 cycles of EM or until convergence

[Mu,Cov,P,Pi,LL]=hmm(Xtrain,200,3,30); 

% plot log likelihood (log L) per sample

subplot(212);
plot(LL/200);                          
ylabel('Log likelihood per sample'); 
xlabel('Iterations of EM'); 
hold on;

% calculate log L for test data

lik=hmm_cl(Xtest,95,3,Mu,Cov,P,Pi);  

plot(length(LL),lik/95,'go');

% examine state transition matrix and plot means

P
                   
subplot(211);
plot(Mu(:,1),Mu(:,2),'r*');
echo off;

hold off;