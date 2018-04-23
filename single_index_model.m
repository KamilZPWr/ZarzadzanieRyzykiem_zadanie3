close all

index_price = importdata('oil_index_month_average_from_01-04-08_to_01_02_18.txt');
oil_data = importdata('oreln_oil_average_price_from_01-04-08_to_01-02-18.txt');
brent_price = oil_data(:,1);
ural_price = oil_data(:,2);

R_brent =(brent_price(2:end)-brent_price(1:end-1))./brent_price(1:end-1)*100;
R_ural = (ural_price(2:end)-ural_price(1:end-1))./ural_price(1:end-1)*100;

R_brent = R_brent(9:end);  
R_ural = R_ural(9:end);

var_brent = var(R_brent);
var_ural = var(R_ural);

std_brent = sqrt(var_brent);
std_ural = sqrt(var_ural);

%cov_brent_ural = max(cov(R_brent, R_ural)); 

R_index = (index_price(2:end)-index_price(1:end-1))./index_price(1:end-1)*100;
R_index = R_index(9:end);

beta_brent = sum((R_brent - mean(R_brent)).*(R_index-mean(R_index)))/sum((R_index-mean(R_index)).^2)
beta_ural = sum((R_ural - mean(R_ural)).*(R_index-mean(R_index)))/sum((R_index-mean(R_index)).^2)

alfa_brent = mean(R_brent) - beta_brent.*mean(R_index)
alfa_ural =  mean(R_ural) - beta_brent.*mean(R_index)

var_random_fractor_brent = var(R_brent - (alfa_brent + beta_brent.*R_index))
var_random_fractor_ural = var(R_ural - (alfa_ural + beta_ural.*R_index))

sharps_index_brent = (mean(R_brent) - mean(R_index))/sqrt(var(R_brent))
sharps_index_ural = (mean(R_ural) - mean(R_index))/sqrt(var(R_ural))

coefficient_of_determination_brent = beta_brent*var(R_index)/var(R_brent)
coefficient_of_determination_ural = beta_ural*var(R_index)/var(R_ural)

weight_brent = 0.86
weight_ural = 0.14

% 1a)
beta_p0 = weight_brent*beta_brent + weight_ural*beta_ural;
E_p_real = weight_brent*alfa_brent + weight_ural*alfa_ural + mean(R_index) * beta_p0
Var_p_real = beta_p0^2 * var(R_index) + weight_brent*var_random_fractor_brent + weight_ural*var_random_fractor_ural

% 1b)
i = 0:0.1:1
figure
hold on
title('Mo¿liwe konstrukcje portfela')
xlabel('Odchylenie standardowe')
ylabel('Oczekiwana stopa zwrotu')
for i = i
    beta_p = i*beta_brent + (1-i)*beta_ural;
    E_p = i*alfa_brent + (1-i)*alfa_ural + mean(R_index) * beta_p
    Var_p = beta_p^2 * var(R_index) + i.^2*var_random_fractor_brent + (1-i).^2*var_random_fractor_ural
    plot(sqrt(Var_p), E_p, 'rx')
end


% 1c)
%weight_brent_min_risk = (var_ural - std_brent*std_ural*cov_brent_ural)/(var_brent + var_ural - 2*std_ural*std_brent*cov_brent_ural)
weight_brent_min_risk = (var(R_index)*beta_ural*(beta_ural-beta_brent)+var_random_fractor_ural)/(var(R_index)*(beta_brent-beta_ural).^2+var_random_fractor_brent+var_random_fractor_ural)
weight_ural_min_risk = 1 - weight_brent_min_risk

figure
w = 0:0.001:1;
sigma_p_w = (w*beta_brent+(1-w)*beta_ural).^2*var(R_index)+w.^2*var_random_fractor_brent+(1-w).^2*var_random_fractor_ural;
sigma_p_wbmr = (weight_brent_min_risk*beta_brent+(1-weight_brent_min_risk)*beta_ural).^2*var(R_index)+weight_brent_min_risk.^2*var_random_fractor_brent+(1-weight_brent_min_risk).^2*var_random_fractor_ural
plot(w,sigma_p_w)
hold on
plot(weight_brent_min_risk, sigma_p_wbmr,'.','markersize',20)

% 1d)

fixedrate = 0.8;
weight_brent_fixedrate = (fixedrate-mean(R_ural))/(mean(R_brent)-mean(R_ural))
weight_ural_fixedrate = 1-weight_brent_fixedrate

figure
w = 0:0.001:1;
sigma_p_w = (w*beta_brent+(1-w)*beta_ural).^2*var(R_index)+w.^2*var_random_fractor_brent+(1-w).^2*var_random_fractor_ural;
sigma_p_wbmr = (weight_brent_min_risk*beta_brent+(1-weight_brent_min_risk)*beta_ural).^2*var(R_index)+weight_brent_min_risk.^2*var_random_fractor_brent+(1-weight_brent_min_risk).^2*var_random_fractor_ural
plot(w,sigma_p_w)
hold on
plot(weight_brent_min_risk, sigma_p_wbmr,'.')

% 1e)

riskfree_rate = 0.175;

w = 0:0.001:1;
sharps_index = 0;
market_portfolio_brent_weight = 0;
for i = w
    beta_p = i*beta_brent + (1-i)*beta_ural;
    E_p = i*alfa_brent + (1-i)*alfa_ural + mean(R_index) * beta_p;
    Var_p = beta_p^2 * var(R_index) + i.^2*var_random_fractor_brent + (1-i).^2*var_random_fractor_ural;
    
    if (sharps_index < (E_p - riskfree_rate)/sqrt(Var_p))
        sharps_index = (E_p - riskfree_rate)/sqrt(Var_p)
        market_portfolio_brent_weight = i;
        market_portfolio_ural_weight = 1-i;
    end
end
sharps_index
market_portfolio_brent_weight
market_portfolio_ural_weight

w = 0:0.001:1;
beta_p = w*beta_brent + (1-w)*beta_ural;
E_p = w*alfa_brent + (1-w)*alfa_ural + mean(R_index) * beta_p;
Var_p = beta_p.^2 * var(R_index) + w.^2*var_random_fractor_brent + (1-w).^2*var_random_fractor_ural;
sharpe_index = (E_p - riskfree_rate)./sqrt(Var_p);
figure 
plot(w,sharpe_index)

