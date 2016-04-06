% SSA_matlab.m
%
% Audrey Lustig and Cedric Gaucherel, 2010-2011
%
% computes SSA according to:
%
% For more details see      Cyclostationnarity analysis of enso memory
%                           Student research report - 2011
%                           Audrey LUSTIG
%
%
% inputs:  - name: name of the analysed signal
%		 - M : Wawelet length


clear all;
close all;
clc;

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LECTURE DES DONNEES ET PARAMETRES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% lecture des données
name='SOI.txt';
[S,Time]=readFile(name);

% Noramlisation
X = S - mean(S); 
X = X/std(X);


% Paramètre d'entrée
M =60; % Paramètre à entrer par l'utilisateur

% Paramètre fixe
N=length(X);
iunit = 'months'; %(pour affichage)
fspec=1;
itime=[1:N];




%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGO SSA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Etape 1: reorganisation des données 
% Faire une matrice de trajectoire de taille [M * N-M+1], M etant 
% le paramètre choisi par l'utilisateur
% Exemple : X = 1 , 2 , 3 , 4, 5, 6 et M = 3
% D= [ 1 , 2 , 3 , 4
%      2 , 3 , 4 , 5
%      3 , 4 , 5 , 6 ]

if M > N; error(['Window length greater than series length! M should be less than N']); return; end;
D = NaN .* ones(M,N-M+1);
for i=1:M
  D(i,:)=X(1,i:N-M+i); 
end                



%%% Etape 2 : gestion des valeurs manquantes (cela devrait être bon)
% Now, do the covariance matrix calculation;
% First, find NaNs and calculate sample size N (not unbiased following Ghil et al. 2002)
nnan = ~isnan(D); 
% renvoie une matrice de boolean avec des 1 pour les valeurs non manquantes 
% et des zeros pour les valeurs manquantes
xnan = ind2sub(size(nnan),nnan); 
% transforme la matrice de boolean en matrice de scalaire
xsize = ((xnan*xnan')); 
% Renvoie un matrice de facteurs de noramlisation utilisé pour ne pas biaiser l'estimateur 
% de la matrice de variance covariance C (following Ghil et al. 2002)) 
% Attention il s'agit bien d'un produit matriciel ici.

% Set NaNs in trajectory matrix to zero if there are any
D(find(isnan(D))) = 0;




%%% Etape 3: Calcul de la matrice de variance covariance (C)
% now we can calculate the covariance matrix, again _not_ unbiased 
% following Ghil et al. 2002)
C = (D*D').*(1./xsize); 
% attention il s'agit ici d'un prduit terme à terme et non matriciel 
% il existe des alternative pour cet étape que l'on pourra envisager si le
% temps nous le permet




%%% Etape 4 : decomposition en valeur singulière de la matrice C
% Find Eigenvalues and Eigenvectors using Singular Value Decomposition of
% the covariance matrix
[U,S,V] = svd(C);  
% produces a diagonal matrix S, of the same 
% dimension as C and with nonnegative diagonal elements in
% decreasing order, and unitary matrices U and V so that
% C = U*S*V'.




%%% Etape 5: Classifification des valeurs manquantes en ordre décroissant
% get eigenvalues and fraction variance explained
evec = U; % vecteurs propres
eval = diag(S);  % valeurs propres
varexp = diag(S)/trace(S); % valeurs propres noramlisées
[eval_sort,ind]=sort(varexp,'descend' ); % Valeurs propres triées en ordre
                                         % décroissant


                                         
%%% Etape 6 : Calcul des composantes principales puis reconstruction du signal inital
% Make Principal Component time series, A 
A = U'*D;
% Now make the reconstructed components
RC = zeros(1,N);
for i=1:M
    RC(i,:) = conv(U(:,i),A(i,:));
end
% Now, normalize the RCs appropriately 
% (Ghil et al. 2002, with error corrected)
for i=1:N
    if i >= 1 && i <= M-1
        RC(:,i) = RC(:,i)*(1/i);
    elseif i <= N-M+1 && i >= M
        RC(:,i) = RC(:,i)./M;
    else
        RC(:,i) = RC(:,i)/(N-i+1);
    end
end


%% Graphiques                                         
% Graphe des valeurs singulières (sans l'intervalle de confiance ici)
figure(1);
clf;
plot(eval_sort,'x-');
set(gca,'XLim',[1 length(eval_sort)]);
xlabel('Valeurs singulières');
ylabel('Explication de la variance totale');
title('SSA - Singular Spectrum analysis')
grid;

% Graphe des vecteurs propres en pair
figure (2);
clf;
subplot(3,1,1);
plot(evec(:,1:2), '-');
legend('1', '2');
title('Eigen vectors 1 & 2')

subplot(3,1,2);
plot(evec(:,3:4), '-');
legend('3', '4');
title('Eigen vectors 3 & 4')

subplot(3,1,3);
plot(evec(:,5:6), '-');
legend('5', '6');
title('Eigen vectors 5 & 6')

% Graphique des éléments reconstruits
figure(3);
clf;
for i=1:7
  subplot(7,1,i);
  plot(RC(:,i),'k-');
end;

% Graphiques: Serie temporelle + éléments reconstruits
figure(4);
clf;
subplot(3,1,1)
mm = ones(1,N)*mean(X); % testing for randomness 
plot(Time,X);
hold on
pmean=plot(Time,mm,'g');
hold off
set(gcf,'Name','Southern Oscillation Index (SOI)');
set(gca,'XLim',[Time(1) Time(N)]);
set(gca,'YLim',[min(X) max(X)]);
Xlabel('Date (Year)'); Ylabel('Amplitude'); title('Southern Oscillation Index (SOI) : from NOAA/CPC Jan to Dec')
legend([pmean(1)],'mean');
grid;

subplot(3,1,2)
plot(Time,X,'b');
hold on
plot(Time,sum(RC(1:2,:),1),'r-','LineWidth',2);% seulement pour les deux 1ere VP
legend('Original','5 RCs');
hold off
set(gca,'XLim',[Time(1) Time(N)]);
Xlabel('Date (Year)'); Ylabel('Amplitude');


subplot(3,1,3)
plot(Time,X,'b');
hold on
plot(Time,sum(RC(1:7,:),1),'r-','LineWidth',2); % Avce 7 valeurs propres
legend('Original','7 RCs');
hold off
set(gca,'XLim',[Time(1) Time(N)]);
Xlabel('Date (Year)'); Ylabel('Amplitude');






%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERVALLES DE CONFIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Monte Carlo (MC) based significance test for white noise spectrums

nb_simul=10; % nombre de simulation MC
for i=1:nb_simul  
a=randn(1,1)*10
randn('seed',floor(sum(10*clock)/a)); %pour changer la seed
R=randn(1,N)*sqrt(nanstd(X)^2); % genère un bruit blanc de variance égale au signal étudié

  % Etape 1 : calcul de la matrice de trajectoire
  for j=1:M
      Dr(j,:)=R(:,j:N-M+j); 
  end
  
% Etape 2 : Gestion des valeurs manquantes
% il est inutile ici de traiter les valeurs manquantes car on crée des 
% bruits blancs sans valeurs manquantes

%  Etape 3 : Calcul de la matrice de variance covariance (C)
Cr = (Dr*Dr').*(1./(N-M+1)); 

% Etape 4 (differente de la precedente): calul des valeurs propres avec la
% même matrice U que celle du jeu de données initial
varexpr = diag(U'*Cr*U)/trace(U'*Cr*U)

[evalr,ind]=sort(varexpr,'descend' );
lambda_white(i,:) = evalr;
clear evalr Cr Dr varexpr R 
end


% Graphique des valeurs propres avec l'intervalle de confiance à 99% pour
% du bruit blanc.
figure(1);
clf;
plot(eval_sort,'x-');
set(gca,'XLim',[1 length(eval_sort)]);
hold on
plot(quantile(lambda_white,0.99),'r-');
hold off
xlabel('Rank');
ylabel('Variance');
title('SSA - Singular Spectrum analysis')
legend('Singular value spectrum', '99% confident level')
grid;








%%%%%%%%%%%%%%%%%%%

%  Monte Carlo (MC) based significance test for red noise spectrums
% Preliminaire:
% Now get simple first order autocorrelation and variance for Monte Carlo tests
temp=X(~isnan(X)); % eliminer les valeurs manquantes pour le calcul de l'autocorrelation
[ACF,lags,bounds] = autocorr(temp,1);
foac=ACF(2);
%disp(['First Order Autocorrelation calculated: ', num2str(foac)])
%disp(['Variance of original datset calculated: ', num2str(nanstd(temp)^2)])


nb_simul=10;
for i=1:nb_simul
a=randn(1,1)*10;
randn('seed',floor(sum(10*clock)/a)); % pour changer la seed

% Creation d'un bruit rouge
G(1,1) = randn(1,1);
for k=1:N
  G(1,k+1) = (foac*G(1,k)) + ((1-foac)*randn(1,1));
end

G = G - repmat(mean(G),1,1);
G = G/std(G);

for j=1:N
  G(1,j) = G(1,j) .* sqrt(nanstd(X)^2);
end
 
for j=1:M
  Dr(j,:)=G(:,j:N-M+j); 
end

Cr = (Dr*Dr').*(1./(N-M+1));
varexpr = diag(U'*Cr*U)'/trace(U'*Cr*U)';
[evalr,ind]=sort(varexpr,'descend');
lambda_red(i,:) = evalr ;
clear evalr Cr Dr evalr varexpr
end


% Graphique des valeurs propres avec l'intervalle de confiance à 99% pour
% du bruit blanc.
figure(1);
clf;
plot(eval_sort,'x-');
set(gca,'XLim',[1 length(eval_sort)]);
hold on
plot(quantile(lambda_red,0.99),'r-');
hold off
xlabel('Rank');
ylabel('Variance');
title('SSA - Singular Spectrum analysis')
legend('Singular value spectrum', '99% confident level')
grid;



%%
% A faire : une analyse de Fourier ou de la fréquence dominante des
% elements reconstruit choisi par l'utilisateur
