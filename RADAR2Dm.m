clear all; close all;clc;
%% Zaman parametreleri ve gürültü büyüklükleri
T = 0.01; %Zaman adımı
N = 1000;  %Örnekleme sayısı
sigma_x = 10;
sigma_y = 10;

%% Durum geçiş matrisi A ve Gamma
A = [1 T 0 0;...
     0 1 0 0;...
     0 0 1 T;...
     0 0 0 1];  

Gamma = [T^2/2 0;...
         T 0;...
         0 T^2/2;...
         0 T]*0;
%% Başlangıç durumu (true states)
x = [-500 0 -500 100]';  % Başlangıç durum vektörü [x, x_dot, y, y_dot]
true_states = x;

%% True state'lerin üretimi (gerçek rotayı oluşturma)
for k = 1:N
    % Yeni durumu hesapla (sistem dinamikleri ve rastgele gürültü eklenir)
    x = A*x + Gamma*[randn*sigma_x; randn*sigma_y];                
    true_states = [true_states x];
end

%% Gerçek durumu çizdir (x vs y yörüngesi)
figure;
plot(true_states(1,:), true_states(3,:), 'b-');
xlabel('x ekseni');
ylabel('y ekseni');
title('Gerçek Durum Trajektorisi');
%% Ölçüm gürültüleri
rx_sigma = 50; % r ölçüm hatası
ry_sigma = 2; % theta ölçüm hatası
measurements = []; % Ölçümler
for i = 1:N
    % Gerçek konumu al
    x_pos = true_states(1,i); 
    y_pos = true_states(3,i); 
    
    % Polar ölçümlere çevir
    r = sqrt(x_pos^2 + y_pos^2); 
    theta_radians = atan2d(y_pos, x_pos); % atan2 fonksiyonu kullanılır
    
    % Ölçümlere gürültü ekle
    r_noisy = r + randn * rx_sigma; 
    theta_noisy = theta_radians + randn*ry_sigma;
    
    % Ölçüm verisini kaydet
    measurements = [measurements, [r_noisy; theta_noisy]];
end
figure();
plot(measurements(:),'kx');
%% Kalman filtresi parametreleri
x_est = [-500 0 -500 100]'; % Başlangıç tahmini
P_est = diag([10^2 5^2 10^2 5^2]); % Başlangıç kovaryans matrisi
Q = [T^4/4*sigma_x^2 T^3/2*sigma_x^2 0 0;...
     T^3/2*sigma_x^2 T^2*sigma_x^2 0 0;...
     0 0 T^4/4*sigma_y^2 T^3/2*sigma_y^2;...
     0 0 T^3/2*sigma_y^2 T^2*sigma_y^2];
R = [rx_sigma^2 0; 0 ry_sigma^2]; % Ölçüm gürültüsü kovaryans matrisi

%% Filtre döngüsü
estimated_states = [];

for i = 1:N
    % Adım 1: Tahmin
    x_pred = A*x_est; 
    P_pred = A*P_est*A' + Q;     

    % Tahmini ölçümleri hesapla
    r_pred = sqrt(x_pred(1)^2 + x_pred(3)^2);
    theta_pred = atan2d(x_pred(3), x_pred(1));
    h_x_pred = [r_pred; theta_pred];

    % Ölçüm inovasyonu
    z_actual = measurements(:,i);
    v = z_actual - h_x_pred;
    % Açı farkını [-180, 180] aralığına getir
    v(2) = mod(v(2) + 180, 360) - 180;

    % Ölçüm matrisi H'yi hesapla
    H = [(x_pred(1)/r_pred), 0, (x_pred(3)/r_pred), 0;
         (-x_pred(3)/(r_pred^2)), 0, (x_pred(1)/(r_pred^2)), 0];

    % Kalman kazancı
    S = H*P_pred*H' + R;
    G = P_pred*H' / S;

    % Durum güncellemesi
    x_est = x_pred + G*v;
    P_est = P_pred - G*S*G';

    % Tahmin edilen ve gerçek pozisyonları çiz
    plot(x_est(1), x_est(3), 'rx', 'MarkerSize', 5,'DisplayName','Estimated Values'); % Tahmin edilen konum
    hold on;
    x_meas = z_actual(1) * cosd(z_actual(2));
    y_meas = z_actual(1) * sind(z_actual(2));
    plot(x_meas, y_meas, 'ko', 'MarkerSize', 5,'DisplayName','Measured Values'); % Ölçüm verisi
    plot(true_states(1, i+1), true_states(3, i+1), 'm+', 'MarkerSize', 5,'DisplayName','True Values' ); % Gerçek durum
    pause(0.1);
    grid on;
   
    % Tahmin edilen durumu kaydet
    estimated_states = [estimated_states x_est];
end

    
    % Tahmin edilen durumu kaydet
    

%% x ve y pozisyonları için karşılaştırma ve hata grafikleri
figure;

% x pozisyonu karşılaştırması
subplot(3,2,1);
plot(1:N+1, true_states(1,:), 'b-', 'DisplayName', 'Gerçek x Pozisyonu');
hold on;
plot(1:N, estimated_states(1,:), 'r--', 'DisplayName', 'Tahmin Edilen x Pozisyonu');
xlabel('Zaman');
ylabel('x Pozisyonu');
title('x Pozisyonu Karşılaştırması');
legend('show');
grid on;

% x pozisyonu hatası (gerçek - tahmin)
subplot(3,2,2);
plot(1:N, true_states(1,1:N) - estimated_states(1,:), 'k-', 'DisplayName', 'x Pozisyonu Hatası');
xlabel('Zaman');
ylabel('Hata');
title('x Pozisyonu Hatası');
legend('show');
grid on;

% y pozisyonu karşılaştırması
subplot(3,2,3);
plot(1:N+1, true_states(3,:), 'b-', 'DisplayName', 'Gerçek y Pozisyonu');
hold on;
plot(1:N, estimated_states(3,:), 'r--', 'DisplayName', 'Tahmin Edilen y Pozisyonu');
xlabel('Zaman');
ylabel('y Pozisyonu');
title('y Pozisyonu Karşılaştırması');
legend('show');
grid on;

% y pozisyonu hatası (gerçek - tahmin)
subplot(3,2,4);
plot(1:N, true_states(3,1:N) - estimated_states(3,:), 'k-', 'DisplayName', 'y Pozisyonu Hatası');
xlabel('Zaman');
ylabel('Hata');
title('y Pozisyonu Hatası');
legend('show');
grid on;

% x hızı karşılaştırması
subplot(3,2,5);
plot(1:N+1, true_states(2,:), 'b-', 'DisplayName', 'Gerçek x Hızı');
hold on;
plot(1:N, estimated_states(2,:), 'r--', 'DisplayName', 'Tahmin Edilen x Hızı');
xlabel('Zaman');
ylabel('x Hızı');
title('x Hızı Karşılaştırması');
legend('show');
grid on;

% x hızı hatası (gerçek - tahmin)
subplot(3,2,6);
plot(1:N, true_states(2,1:N) - estimated_states(2,:), 'k-', 'DisplayName', 'x Hızı Hatası');
xlabel('Zaman');
ylabel('Hata');
title('x Hızı Hatası');
legend('show');
grid on;

%% y hızı karşılaştırması ve hatası
figure;

% y hızı karşılaştırması
subplot(2,1,1);
plot(1:N+1, true_states(4,:), 'b-', 'DisplayName', 'Gerçek y Hızı');
hold on;
plot(1:N, estimated_states(4,:), 'r--', 'DisplayName', 'Tahmin Edilen y Hızı');
xlabel('Zaman');
ylabel('y Hızı');
title('y Hızı Karşılaştırması');
legend('show');
grid on;

% y hızı hatası (gerçek - tahmin)
subplot(2,1,2);
plot(1:N, true_states(4,1:N) - estimated_states(4,:), 'k-', 'DisplayName', 'y Hızı Hatası');
xlabel('Zaman');
ylabel('Hata');
title('y Hızı Hatası');
legend('show');
grid on;