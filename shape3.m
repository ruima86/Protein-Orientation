function dy = shape3(u,y,para,fSharp,ag1,ag2,C01,kappa0,kappa2,p,ap,theta)

C02=para(2);


dy = zeros(9,1);  
psi = y(1);   
dpsi = y(2);
r = y(3);
z = y(4);
alpha = y(5);
beta = y(6);
h = y(7);
gamma = y(8);
area = y(9);


dy(1)=dpsi;

dy(2)=[2.*beta.*h.^2.*r.^(-1).*cos(psi).*(2.*kappa0+ap.*kappa2.*sin( ...
  theta).^4.*tanh(((-1).*ag1+area).*fSharp)+(-1).*ap.*kappa2.*sin( ...
  theta).^4.*tanh(((-1).*ag2+area).*fSharp)).^(-1)+(-2).*dpsi.*h.* ...
  kappa0.*r.^(-1).*cos(psi).*(2.*kappa0+ap.*kappa2.*sin(theta).^4.* ...
  tanh(((-1).*ag1+area).*fSharp)+(-1).*ap.*kappa2.*sin(theta).^4.* ...
  tanh(((-1).*ag2+area).*fSharp)).^(-1)+h.^2.*p.*r.*cos(psi).*(2.* ...
  kappa0+ap.*kappa2.*sin(theta).^4.*tanh(((-1).*ag1+area).*fSharp)+( ...
  -1).*ap.*kappa2.*sin(theta).^4.*tanh(((-1).*ag2+area).*fSharp)).^( ...
  -1)+2.*alpha.*h.^2.*r.^(-1).*sin(psi).*(2.*kappa0+ap.*kappa2.*sin( ...
  theta).^4.*tanh(((-1).*ag1+area).*fSharp)+(-1).*ap.*kappa2.*sin( ...
  theta).^4.*tanh(((-1).*ag2+area).*fSharp)).^(-1)+2.*h.^2.*kappa0.* ...
  r.^(-2).*cos(psi).*sin(psi).*(2.*kappa0+ap.*kappa2.*sin(theta) ...
  .^4.*tanh(((-1).*ag1+area).*fSharp)+(-1).*ap.*kappa2.*sin(theta) ...
  .^4.*tanh(((-1).*ag2+area).*fSharp)).^(-1)+ap.*C02.*fSharp.*h.^2.* ...
  kappa2.*r.*sech(((-1).*ag1+area).*fSharp).^2.*sin(theta).^2.*(2.* ...
  kappa0+ap.*kappa2.*sin(theta).^4.*tanh(((-1).*ag1+area).*fSharp)+( ...
  -1).*ap.*kappa2.*sin(theta).^4.*tanh(((-1).*ag2+area).*fSharp)).^( ...
  -1)+(-1).*ap.*C02.*fSharp.*h.^2.*kappa2.*r.*sech(((-1).*ag2+area) ...
  .*fSharp).^2.*sin(theta).^2.*(2.*kappa0+ap.*kappa2.*sin(theta) ...
  .^4.*tanh(((-1).*ag1+area).*fSharp)+(-1).*ap.*kappa2.*sin(theta) ...
  .^4.*tanh(((-1).*ag2+area).*fSharp)).^(-1)+(-1).*ap.*fSharp.* ...
  h.^2.*kappa2.*cos(theta).^2.*sech(((-1).*ag1+area).*fSharp).^2.* ...
  sin(psi).*sin(theta).^2.*(2.*kappa0+ap.*kappa2.*sin(theta).^4.* ...
  tanh(((-1).*ag1+area).*fSharp)+(-1).*ap.*kappa2.*sin(theta).^4.* ...
  tanh(((-1).*ag2+area).*fSharp)).^(-1)+ap.*fSharp.*h.^2.*kappa2.* ...
  cos(theta).^2.*sech(((-1).*ag2+area).*fSharp).^2.*sin(psi).*sin( ...
  theta).^2.*(2.*kappa0+ap.*kappa2.*sin(theta).^4.*tanh(((-1).*ag1+ ...
  area).*fSharp)+(-1).*ap.*kappa2.*sin(theta).^4.*tanh(((-1).*ag2+ ...
  area).*fSharp)).^(-1)+(-1).*ap.*dpsi.*fSharp.*h.*kappa2.*r.*sech(( ...
  (-1).*ag1+area).*fSharp).^2.*sin(theta).^4.*(2.*kappa0+ap.* ...
  kappa2.*sin(theta).^4.*tanh(((-1).*ag1+area).*fSharp)+(-1).*ap.* ...
  kappa2.*sin(theta).^4.*tanh(((-1).*ag2+area).*fSharp)).^(-1)+ap.* ...
  dpsi.*fSharp.*h.*kappa2.*r.*sech(((-1).*ag2+area).*fSharp).^2.* ...
  sin(theta).^4.*(2.*kappa0+ap.*kappa2.*sin(theta).^4.*tanh(((-1).* ...
  ag1+area).*fSharp)+(-1).*ap.*kappa2.*sin(theta).^4.*tanh(((-1).* ...
  ag2+area).*fSharp)).^(-1)+(-1).*ap.*C02.*h.^2.*kappa2.*r.^(-1).* ...
  cos(psi).*cos(theta).^2.*tanh(((-1).*ag1+area).*fSharp).*(2.* ...
  kappa0+ap.*kappa2.*sin(theta).^4.*tanh(((-1).*ag1+area).*fSharp)+( ...
  -1).*ap.*kappa2.*sin(theta).^4.*tanh(((-1).*ag2+area).*fSharp)).^( ...
  -1)+ap.*h.^2.*kappa2.*r.^(-2).*cos(psi).*cos(theta).^4.*sin(psi).* ...
  tanh(((-1).*ag1+area).*fSharp).*(2.*kappa0+ap.*kappa2.*sin(theta) ...
  .^4.*tanh(((-1).*ag1+area).*fSharp)+(-1).*ap.*kappa2.*sin(theta) ...
  .^4.*tanh(((-1).*ag2+area).*fSharp)).^(-1)+ap.*C02.*h.^2.*kappa2.* ...
  r.^(-1).*cos(psi).*sin(theta).^2.*tanh(((-1).*ag1+area).*fSharp).* ...
  (2.*kappa0+ap.*kappa2.*sin(theta).^4.*tanh(((-1).*ag1+area).* ...
  fSharp)+(-1).*ap.*kappa2.*sin(theta).^4.*tanh(((-1).*ag2+area).* ...
  fSharp)).^(-1)+(-1).*ap.*dpsi.*h.*kappa2.*r.^(-1).*cos(psi).*sin( ...
  theta).^4.*tanh(((-1).*ag1+area).*fSharp).*(2.*kappa0+ap.*kappa2.* ...
  sin(theta).^4.*tanh(((-1).*ag1+area).*fSharp)+(-1).*ap.*kappa2.* ...
  sin(theta).^4.*tanh(((-1).*ag2+area).*fSharp)).^(-1)+ap.*C02.* ...
  h.^2.*kappa2.*r.^(-1).*cos(psi).*cos(theta).^2.*tanh(((-1).*ag2+ ...
  area).*fSharp).*(2.*kappa0+ap.*kappa2.*sin(theta).^4.*tanh(((-1).* ...
  ag1+area).*fSharp)+(-1).*ap.*kappa2.*sin(theta).^4.*tanh(((-1).* ...
  ag2+area).*fSharp)).^(-1)+(-1).*ap.*h.^2.*kappa2.*r.^(-2).*cos( ...
  psi).*cos(theta).^4.*sin(psi).*tanh(((-1).*ag2+area).*fSharp).*( ...
  2.*kappa0+ap.*kappa2.*sin(theta).^4.*tanh(((-1).*ag1+area).* ...
  fSharp)+(-1).*ap.*kappa2.*sin(theta).^4.*tanh(((-1).*ag2+area).* ...
  fSharp)).^(-1)+(-1).*ap.*C02.*h.^2.*kappa2.*r.^(-1).*cos(psi).* ...
  sin(theta).^2.*tanh(((-1).*ag2+area).*fSharp).*(2.*kappa0+ap.* ...
  kappa2.*sin(theta).^4.*tanh(((-1).*ag1+area).*fSharp)+(-1).*ap.* ...
  kappa2.*sin(theta).^4.*tanh(((-1).*ag2+area).*fSharp)).^(-1)+ap.* ...
  dpsi.*h.*kappa2.*r.^(-1).*cos(psi).*sin(theta).^4.*tanh(((-1).* ...
  ag2+area).*fSharp).*(2.*kappa0+ap.*kappa2.*sin(theta).^4.*tanh((( ...
  -1).*ag1+area).*fSharp)+(-1).*ap.*kappa2.*sin(theta).^4.*tanh((( ...
  -1).*ag2+area).*fSharp)).^(-1)];

dy(3)=[h.*cos(psi)];

dy(4)=[(-1).*h.*sin(psi)];

dy(5)=[gamma.*h+(-1).*C01.*dpsi.*kappa0+(1/2).*dpsi.^2.*h.^(-1).*kappa0+ ...
  (1/2).*C01.^2.*h.*kappa0+h.*p.*r.*sin(psi)+(-1/2).*h.*kappa0.*r.^( ...
  -2).*sin(psi).^2+(1/4).*ap.*C02.^2.*h.*kappa2.*tanh(((-1).*ag1+ ...
  area).*fSharp)+(-1/4).*ap.*h.*kappa2.*r.^(-2).*cos(theta).^4.*sin( ...
  psi).^2.*tanh(((-1).*ag1+area).*fSharp)+(-1/2).*ap.*C02.*dpsi.* ...
  kappa2.*sin(theta).^2.*tanh(((-1).*ag1+area).*fSharp)+(1/4).*ap.* ...
  dpsi.^2.*h.^(-1).*kappa2.*sin(theta).^4.*tanh(((-1).*ag1+area).* ...
  fSharp)+(-1/4).*ap.*C02.^2.*h.*kappa2.*tanh(((-1).*ag2+area).* ...
  fSharp)+(1/4).*ap.*h.*kappa2.*r.^(-2).*cos(theta).^4.*sin(psi) ...
  .^2.*tanh(((-1).*ag2+area).*fSharp)+(1/2).*ap.*C02.*dpsi.*kappa2.* ...
  sin(theta).^2.*tanh(((-1).*ag2+area).*fSharp)+(-1/4).*ap.* ...
  dpsi.^2.*h.^(-1).*kappa2.*sin(theta).^4.*tanh(((-1).*ag2+area).* ...
  fSharp)];

dy(6)=[0];

dy(7)=0;

dy(8)=[(-1/4).*ap.*fSharp.*h.^(-1).*kappa2.*r.^(-1).*(sech(((-1).*ag1+ ...
  area).*fSharp).^2+(-1).*sech(((-1).*ag2+area).*fSharp).^2).*(h.*( ...
  C02.*r+(-1).*cos(theta).^2.*sin(psi))+(-1).*dpsi.*r.*sin(theta) ...
  .^2).^2];

dy(9)=[h.*r];
