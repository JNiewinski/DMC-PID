clear all;
%%  wczytywanie danych 
dane18=dlmread('dane18.txt');
u=dane18(:,1);
y=dane18(:,2);
%% wybór regulatora 
dl=200; %dlugosc symulacji
%%%%%%%%%%%%%%%%%%%%%% wybor regulatora %%%%%%%%%%%%%%%%%%%%%%
% regulator='PID';
regulator='DMC';
% regulator='DMCzo'; 
%DMCzo - DMC z ograniczeniami wartosci sygnalu sterujacego do zad6
ogru=0;
ogrdu=0;
% wlaczenie ograniczen (0-wylaczone, 1-wlaczone)
% ogru - ograniczenia sygnalu sterujacego
% ogrdu - ograniczenia zmiany sygnalu sterujacego
zaklocenia=1; %tylko dla DMC
% wlaczenie zaklocen do zad5 
% zaklocenia=0 jesli brak
tau=1;
%wzmocnienie do zadania dodatkowego
% tau=1 jesli nie ma zmiany wzmocnienia
%% wspolczynniki DMC 
lambda=3000; % wspolczynnik kary
N=30; %horyzont predykcji
Nu=2; %horyzont sterowania
D=22; %horyzont dynamiki
umin=-0.4;                                                                   
umax=0.4;
dumin=-0.01;
dumax=0.01;
% dobor ograniczen
%% wspolczynniki PID 
K=0.5;    % wzmocnienie
Td=1.7;   % czas wyprzedzenia
Ti=11;   % czas zdwojenia
% ogólny wzór pid=k(1+1/Ti*s+Td*s)
% dla eksperymentu Zieglera-Nicholsa Td=0, Ti=inf
%%
yzad=zeros(dl,1); % wartosc zadana
skok=zeros(1,dl); % wektor do badania skoku jednostkowego
yzad(1:10)=0;
yzad(10:dl)=1;
skok(1)=0;
skok(1:dl)=1;
%%  model 
E=0;
ymod=zeros(dl:1);
t=6; %opoznienie
Y=y(2+t : end);
U=[u(2:end-t) u(1:end-t-1) y(1+t:end-1) y(t:end-2)];
wsp=U\Y;
b1=wsp(1);b2=wsp(2);a1=-wsp(3);a2=-wsp(4);
ymod(1:t+1)=0;
% symulacja modelu
for k=t+2:dl
    ymod(k)=b1 * u(k-t) + b2 * u(k-t-1)-a1*ymod(k-1)-a2*ymod(k-2);
    E=E+(ymod(k)-y(k))^2;
end
% transmitancja
H=tf([0,0,0,0,0,b1,b2],[1,a1,a2,0,0,0,0,0],1);
% stairs(ymod);
% hold on;
% stairs(y,':');
% title(strcat('E=',num2str(E),' opoznienie=', num2str(t)));
% legend('ymod','y');
% saveas(figure(1),strcat('model_opoznie=',num2str(t),'.png'));
%% obliczanie parametrów
if strcmp(regulator,'DMC') || strcmp(regulator,'DMCzo')
            
         for k=t+2:dl
            ymod(k)=b1*skok(k-t)+b2*skok(k-t-1)-a1*ymod(k-1)-a2*ymod(k-2);
            s(k)=ymod(k);
         end
    %wzmocnienie statyczne
    wzmocnie_statyczne=s(D);
    s=s(t:end);     % paramatery s
    dU=zeros(Nu,1);
    dUp=zeros(D-1,1);   % macierz poprzednich sterowan
    M=zeros(N,Nu);      % macierz M
    Mp=zeros(N,D-1);% maczierz Mp
        for i=1:N           
            for j=1:Nu
                if i>=j
                M(i,j)=s(i-j+1); %oblicznie M
                end
            end
        end
        for i=1:N
            for j=1:D-1
                if i+j<D
                Mp(i,j)=s(i+j)-s(j); %oblicznie Mp
                else
                Mp(i,j)=s(D)-s(j);
                end
            end
        end
    K=inv(M'*eye*M+eye(Nu)*lambda)*M'*eye; % macierz K
    K=K(1,:); % do regulacji potrzebujemy tylko 1 kolumny
end % koniec DMC
if strcmp(regulator,'PID')
    r0= K * (1+ 1/(2*Ti) + Td);
    r1= K * (1/(2*Ti)-2*Td-1);
    r2= K * Td;
end %koniec PID
%% petla regulacji 
Jy=0;Ju=0; %zmienne do oceny jakosci regulacji
ymod=zeros(dl,1);
e=zeros(dl,1);
uk=zeros(dl,1);
for k=t+2:dl
    
    if strcmp(regulator,'DMC') || strcmp(regulator,'DMCzo')
            % przesuniecie warto?ci poprzednich sterowania
            for i=D-1:-1:2
               dUp(i)=dUp(i-1);
            end
        %symulacja obiektu
        ymod(k)=tau*(b1*uk(k-t)+b2*uk(k-t-1))-a1*ymod(k-1)-a2*ymod(k-2);
        zak(k)=ymod(k);
        %zaklocenia
        ymod(k)=ymod(k)+zaklocenia*(rand-0.5)*0.1;
        zak(k)=ymod(k)-zak(k);
        %warto?? zmiany sterowania
        dU=sum(K)*(yzad(k)-ymod(k))-K*Mp*dUp;

            if strcmp(regulator,'DMCzo') && ogrdu == 1
                if dU(1) >dumax
                    dU(1)=dumax;
                end
                  if dU(1) <dumin
                    dU(1)=dumin;
                  end
            end
                
        %uaktualnienie poprzednich warto?ci sterowania
        dUp(1)=dU(1);
        
        %warto?? sterowania
        uk(k)=uk(k-1)+dU(1);
        
          if strcmp(regulator,'DMCzo') && ogru == 1
                if uk(k) >umax
                    uk(k)=umax;
                end
                  if uk(k) <umin
                    uk(k)=umin;
                  end
         end
        Jy=Jy+(yzad(k)-ymod(k))^2;
        Ju=Ju+(uk(k)-uk(k-1))^2;
    end % koniec DMC
    if strcmp(regulator,'PID')
        
        %symulacja obiektu
        ymod(k)=b1*uk(k-t)+b2*uk(k-t-1)-a1*ymod(k-1)-a2*ymod(k-2);

        % uchyb
        e(k)=ymod(k)-yzad(k);
        
        %warto?? sterowania
        uk(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+uk(k-1);
     Jy=Jy+(yzad(k)-ymod(k))^2;
     Ju=Ju+(uk(k)-uk(k-1))^2;
    end % koniec PID
    
end
%% prezentacja wynikow
subplot(2,1,1);
stairs(ymod);
hold on;
stairs(yzad,'g--');
% stairs(zak,':');
% title('sygnal wyjsciowy');
if strcmp(regulator,'DMC') || strcmp(regulator,'DMCzo')
title(strcat('lambda=',num2str(lambda),' D=',num2str(D),' N=',num2str(N),' Nu=',num2str(Nu),' Ju=',num2str(Ju),' Jy=',num2str(Jy)));
% title(strcat('ogranicznie sygnalu sterujacego=',num2str(umax),' ograniczenie przyrostu=',num2str(dumax)));
else
title(strcat('K= ',num2str(K),' Td=',num2str(Td),' Ti=',num2str(Ti)));
end
hold off;
subplot(2,1,2);
stairs(uk);
title('sygnal sterujacy');
hold off;