clear;
close all;
r0=0.001;
rb=3;    %radius of the initial sphere

p=0;  %-0.555; pressure (out-in)
%p=-2*sigma_e/rb;  % Young-Laplace equation for a sphere
kappa0=1;      % kappa_m bending rigidity of the base
kappa2 = 2;    % kappa_p bending rigidity of the anisotropic proteins

psi0=pi;        % psi(u=1) fixed hinge
h0 = pi*rb;     % guess of the initial arclength 
fSharp=100;     % slope of the tanhx

ag1start = 0.5;  % position of the proteins A_1
ag2start = 0.68; % position of the proteins A_2
ap=1;
ag1 = ag1start;
ag2 = ag2start;
H=2*rb;  %guess of the initial membrane height

C01 = 0;  % spontaneous curvature of the base
theta=0;  % tilt angle theta (begin)
thetaend = 1.57; % tilt angle theta (end)
C02begin = 0;  % spontaneous curvature of the anisotropic proteins (begin)
C02end = 16;   % spontaneous curvature of the anisotropic proteins (end)

C02=C02begin;

nmax = 20;  
RTol = 1e-5;     


yeq = @(u,y,para) shape1(u,y,para,fSharp,ag1,ag2,C01,C02,kappa0,kappa2,p,ap,theta);
ybc = @(ya,yb,para) twobc1(ya,yb,para,r0,rb,psi0,C01,C02,fSharp,ag1,ag2,kappa0,kappa2,p,ap,theta);
yinit = @(u) guess(u,r0,p,rb);
solinit = bvpinit(linspace(0,1,1000),yinit,[h0]) ;  
opts = bvpset('RelTol',RTol,'AbsTol',1e-10,'NMax',50000);
sol = bvp5c(yeq,ybc,solinit,opts);

if sol.stats.maxerr > RTol
    return; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'increase theta'
for theta = 0:0.001:thetaend
    theta
    yeq = @(u,y,para)shape1(u,y,para,fSharp,ag1,ag2,C01,C02,kappa0,kappa2,p,ap,theta);
    ybc = @(ya,yb,para)twobc1(ya,yb,para,r0,rb,psi0,C01,C02,fSharp,ag1,ag2,kappa0,kappa2,p,ap,theta);
    soltemp = bvp5c(yeq,ybc,sol,opts); 

    if soltemp.stats.maxerr < RTol
        j=j+1;
        sol=soltemp;  

    else        
        'Too Large theta'
        break;
    end

    figure(1);
    r = sol.y(3,:);
    z = sol.y(4,:);
    area11=sol.y(9,:);
    areaC0=ag2-ag1;
    position=(ag2+ag1)/2;
    Hprotrusion0=sol.y(4,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dpsi= sol.y(2,:);
    h=sol.parameters(1);
    psi=sol.y(1,:);
    area=sol.y(9,:);
    u=sol.x;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalenergy=(1/2).*h.*p.*r.^2.*sin(psi)+(1/2).*h.*kappa0.*r.*((-1).*C01+dpsi.* ...
  h.^(-1)+r.^(-1).*sin(psi)).^2+(1/4).*ap.*h.*kappa2.*r.*((-1).*C02+ ...
  r.^(-1).*cos(theta).^2.*sin(psi)+dpsi.*h.^(-1).*sin(theta).^2) ...
  .^2.*(tanh(((-1).*ag1+area).*fSharp)+(-1).*tanh(((-1).*ag2+area).* ...
  fSharp));
bendenergy=(1/2).*h.*kappa0.*r.*((-1).*C01+dpsi.*h.^(-1)+r.^(-1).*sin(psi)) ...
  .^2+(1/4).*ap.*h.*kappa2.*r.*((-1).*C02+r.^(-1).*cos(theta).^2.* ...
  sin(psi)+dpsi.*h.^(-1).*sin(theta).^2).^2.*(tanh(((-1).*ag1+area) ...
  .*fSharp)+(-1).*tanh(((-1).*ag2+area).*fSharp));
pressureenergy=(1/2).*h.*p.*r.^2.*sin(psi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,idx1]=min(abs(area11-ag1));   
    [~,idx2]=min(abs(area11-ag2));  
    plot(r(1:idx1),z(1:idx1),'k',-r(1:idx1),z(1:idx1),'k','LineWidth',3);
    hold on
    plot(r(idx1:idx2),z(idx1:idx2),'g',-r(idx1:idx2),z(idx1:idx2),'g','LineWidth',3);
    hold on
    plot(r(idx2:end),z(idx2:end),'k',-r(idx2:end),z(idx2:end),'k','LineWidth',3);
    axis equal;  
   title(['C0_{2} = ',num2str(C02),',','areaC0 = ',num2str(areaC0),newline,'ag_{1} = ',num2str(ag1),',','ag_{2} = ',num2str(ag2),',','position = ',num2str(0.5*(ag1+ag2)),',','H = ',num2str(sol.y(4,1)),',','theta = ',num2str(theta)]);
    xlabel('r','FontSize',26,'FontName','Times');
    ylabel('z','FontSize',26,'FontName','Times');
    set(gca,'FontName','Times New Roman','FontSize',26,'LineWidth',2.5)
    set(figure(1), 'Color', 'white');
    drawnow;
    hold off; 
    areasphere=sol.y(9,end);

end 
thetamax=theta;
'vary curvature'
varc = C02begin:0.01:C02end;
fH = [];
dynamicCURVATURESOL=[];
j = 0;
for C02 = varc
    usefulC02=C02;
    yeq = @(u,y,para)shape1(u,y,para,fSharp,ag1,ag2,C01,C02,kappa0,kappa2,p,ap,theta);
    ybc = @(ya,yb,para)twobc1(ya,yb,para,r0,rb,psi0,C01,C02,fSharp,ag1,ag2,kappa0,kappa2,p,ap,theta);

    soltemp = bvp5c(yeq,ybc,sol,opts); 

    if soltemp.stats.maxerr < RTol
        j=j+1;
        sol=soltemp;  

    else 
        usefulC02=C02-(varc(2)-varc(1)); 
        break;
    end

    figure(1);
    r = sol.y(3,:);
    z = sol.y(4,:);
    area11=sol.y(9,:);
    areaC0=ag2-ag1;
    position=(ag2+ag1)/2;
    Hprotrusion0=sol.y(4,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dpsi= sol.y(2,:);
    h=sol.parameters(1);
    psi=sol.y(1,:);
    area=sol.y(9,:);
    u=sol.x;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalenergy=(1/2).*h.*p.*r.^2.*sin(psi)+(1/2).*h.*kappa0.*r.*((-1).*C01+dpsi.* ...
  h.^(-1)+r.^(-1).*sin(psi)).^2+(1/4).*ap.*h.*kappa2.*r.*((-1).*C02+ ...
  r.^(-1).*cos(theta).^2.*sin(psi)+dpsi.*h.^(-1).*sin(theta).^2) ...
  .^2.*(tanh(((-1).*ag1+area).*fSharp)+(-1).*tanh(((-1).*ag2+area).* ...
  fSharp));
bendenergy=(1/2).*h.*kappa0.*r.*((-1).*C01+dpsi.*h.^(-1)+r.^(-1).*sin(psi)) ...
  .^2+(1/4).*ap.*h.*kappa2.*r.*((-1).*C02+r.^(-1).*cos(theta).^2.* ...
  sin(psi)+dpsi.*h.^(-1).*sin(theta).^2).^2.*(tanh(((-1).*ag1+area) ...
  .*fSharp)+(-1).*tanh(((-1).*ag2+area).*fSharp));
pressureenergy=(1/2).*h.*p.*r.^2.*sin(psi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fH(j,:) = [usefulC02,Hprotrusion0,ag1,ag2,simps(u,totalenergy),simps(u,bendenergy),simps(u,pressureenergy),theta];
    dynamicCURVATURESOL=[dynamicCURVATURESOL sol];
    [~,idx1]=min(abs(area11-ag1)); 
    [~,idx2]=min(abs(area11-ag2));   
    plot(r(1:idx1),z(1:idx1),'k',-r(1:idx1),z(1:idx1),'k','LineWidth',3);
    hold on
    plot(r(idx1:idx2),z(idx1:idx2),'g',-r(idx1:idx2),z(idx1:idx2),'g','LineWidth',3);
    hold on
    plot(r(idx2:end),z(idx2:end),'k',-r(idx2:end),z(idx2:end),'k','LineWidth',3);
    axis equal;  
    %ylim([0 6]);      
    %xlim([-4 4]);
   title(['C0_{2} = ',num2str(usefulC02),',','areaC0 = ',num2str(areaC0),newline,'ag_{1} = ',num2str(ag1),',','ag_{2} = ',num2str(ag2),',','position = ',num2str(0.5*(ag1+ag2)),',','H = ',num2str(sol.y(4,1)),',','theta = ',num2str(theta)]);
    xlabel('r','FontSize',26,'FontName','Times');
    ylabel('z','FontSize',26,'FontName','Times');
    set(gca,'FontName','Times New Roman','FontSize',26,'LineWidth',2.5)
    set(figure(1), 'Color', 'white');
    drawnow;
    hold off; 
    areasphere=sol.y(9,end);

end 
C02=usefulC02; 


'vary H'
h0=sol.parameters(1);
C02=fH(end,1);   
Hinitial=fH(end,2);  
sol.parameters=[h0,C02]; 
%if  fH(end-1,2)<fH(end,2)
%    varc = Hinitial:0.001:10; 
%else
%    varc = Hinitial:-0.001:0; 
%end
varc = Hinitial:0.01:10;  
for H = varc  
    yeq = @(u,y,para) shape3(u,y,para,fSharp,ag1,ag2,C01,kappa0,kappa2,p,ap,theta);
    ybc = @(ya,yb,para)twobc3(ya,yb,para,r0,rb,psi0,C01,fSharp,ag1,ag2,kappa0,kappa2,H,p,ap,theta);
    opts = bvpset('RelTol',1e-6,'AbsTol',1e-10,'NMax',50000);
    soltemp = bvp5c(yeq,ybc,sol,opts);    
    if soltemp.stats.maxerr < RTol
           j=j+1;
          sol=soltemp;       
    figure(1);
    r = sol.y(3,:);
    z = sol.y(4,:);
    area11=sol.y(9,:);
    areaC0=ag2-ag1;
    position=(ag2+ag1)/2;
    Hprotrusion0=sol.y(4,1);
    usefulC02=sol.parameters(2);
    C02= usefulC02;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    dpsi= sol.y(2,:);
    h=sol.parameters(1);
    psi=sol.y(1,:);
    area=sol.y(9,:);
    u=sol.x;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalenergy=(1/2).*h.*p.*r.^2.*sin(psi)+(1/2).*h.*kappa0.*r.*((-1).*C01+dpsi.* ...
  h.^(-1)+r.^(-1).*sin(psi)).^2+(1/4).*ap.*h.*kappa2.*r.*((-1).*C02+ ...
  r.^(-1).*cos(theta).^2.*sin(psi)+dpsi.*h.^(-1).*sin(theta).^2) ...
  .^2.*(tanh(((-1).*ag1+area).*fSharp)+(-1).*tanh(((-1).*ag2+area).* ...
  fSharp));
bendenergy=(1/2).*h.*kappa0.*r.*((-1).*C01+dpsi.*h.^(-1)+r.^(-1).*sin(psi)) ...
  .^2+(1/4).*ap.*h.*kappa2.*r.*((-1).*C02+r.^(-1).*cos(theta).^2.* ...
  sin(psi)+dpsi.*h.^(-1).*sin(theta).^2).^2.*(tanh(((-1).*ag1+area) ...
  .*fSharp)+(-1).*tanh(((-1).*ag2+area).*fSharp));
pressureenergy=(1/2).*h.*p.*r.^2.*sin(psi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fH(j,:) = [usefulC02,Hprotrusion0,ag1,ag2,simps(u,totalenergy),simps(u,bendenergy),simps(u,pressureenergy),theta];
    dynamicCURVATURESOL=[dynamicCURVATURESOL sol];
    [~,idx1]=min(abs(area11-ag1));  
    [~,idx2]=min(abs(area11-ag2));  
    plot(r(1:idx1),z(1:idx1),'k',-r(1:idx1),z(1:idx1),'k','LineWidth',3);
    hold on
    plot(r(idx1:idx2),z(idx1:idx2),'g',-r(idx1:idx2),z(idx1:idx2),'g','LineWidth',3);
    hold on
    plot(r(idx2:end),z(idx2:end),'k',-r(idx2:end),z(idx2:end),'k','LineWidth',3);
    axis equal;  
%     ylim([0 6]);      
    xlim([-4 4]);

    title(['C0_{2} = ',num2str(usefulC02),',','areaC0 = ',num2str(areaC0),newline,'ag_{1} = ',num2str(ag1),',','ag_{2} = ',num2str(ag2),',','position = ',num2str(0.5*(ag1+ag2)),',','H = ',num2str(sol.y(4,1)),',','theta = ',num2str(theta)]);
    xlabel('r','FontSize',26,'FontName','Times');
    ylabel('z','FontSize',26,'FontName','Times');
    set(gca,'FontName','Times New Roman','FontSize',26,'LineWidth',2.5)
    set(figure(1), 'Color', 'white');
    drawnow;
    hold off; 
    areasphere=sol.y(9,end);

    else
        break;
    end

end


'vary curvature 2'
h0=sol.parameters(1);
C02=fH(end,1);   
sol.parameters=[h0]; 
if  fH(end-1,1)<fH(end,1)
    varc = C02:0.001:C02end; 
else
    varc = C02:-0.001:0;  
end
for C02 = varc
    usefulC02=C02;
    yeq = @(u,y,para)shape1(u,y,para,fSharp,ag1,ag2,C01,C02,kappa0,kappa2,p,ap,theta);
    ybc = @(ya,yb,para)twobc1(ya,yb,para,r0,rb,psi0,C01,C02,fSharp,ag1,ag2,kappa0,kappa2,p,ap,theta);

    soltemp = bvp5c(yeq,ybc,sol,opts); 

    if soltemp.stats.maxerr < RTol
        j=j+1;
        sol=soltemp;  

    else 
        usefulC02=C02-(varc(2)-varc(1)); 
        break;
    end

    figure(1);
    r = sol.y(3,:);
    z = sol.y(4,:);
    area11=sol.y(9,:);
    areaC0=ag2-ag1;
    position=(ag2+ag1)/2;
    Hprotrusion0=sol.y(4,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dpsi= sol.y(2,:);
    h=sol.parameters(1);
    psi=sol.y(1,:);
    area=sol.y(9,:);
    u=sol.x;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalenergy=(1/2).*h.*p.*r.^2.*sin(psi)+(1/2).*h.*kappa0.*r.*((-1).*C01+dpsi.* ...
  h.^(-1)+r.^(-1).*sin(psi)).^2+(1/4).*ap.*h.*kappa2.*r.*((-1).*C02+ ...
  r.^(-1).*cos(theta).^2.*sin(psi)+dpsi.*h.^(-1).*sin(theta).^2) ...
  .^2.*(tanh(((-1).*ag1+area).*fSharp)+(-1).*tanh(((-1).*ag2+area).* ...
  fSharp));
bendenergy=(1/2).*h.*kappa0.*r.*((-1).*C01+dpsi.*h.^(-1)+r.^(-1).*sin(psi)) ...
  .^2+(1/4).*ap.*h.*kappa2.*r.*((-1).*C02+r.^(-1).*cos(theta).^2.* ...
  sin(psi)+dpsi.*h.^(-1).*sin(theta).^2).^2.*(tanh(((-1).*ag1+area) ...
  .*fSharp)+(-1).*tanh(((-1).*ag2+area).*fSharp));
pressureenergy=(1/2).*h.*p.*r.^2.*sin(psi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fH(j,:) = [usefulC02,Hprotrusion0,ag1,ag2,simps(u,totalenergy),simps(u,bendenergy),simps(u,pressureenergy),theta];
    dynamicCURVATURESOL=[dynamicCURVATURESOL sol];
    [~,idx1]=min(abs(area11-ag1));  
    [~,idx2]=min(abs(area11-ag2));  
    plot(r(1:idx1),z(1:idx1),'k',-r(1:idx1),z(1:idx1),'k','LineWidth',3);
    hold on
    plot(r(idx1:idx2),z(idx1:idx2),'g',-r(idx1:idx2),z(idx1:idx2),'g','LineWidth',3);
    hold on
    plot(r(idx2:end),z(idx2:end),'k',-r(idx2:end),z(idx2:end),'k','LineWidth',3);
    axis equal;  
    %ylim([0 6]);      
    %xlim([-4 4]);
   title(['C0_{2} = ',num2str(usefulC02),',','areaC0 = ',num2str(areaC0),newline,'ag_{1} = ',num2str(ag1),',','ag_{2} = ',num2str(ag2),',','position = ',num2str(0.5*(ag1+ag2)),',','H = ',num2str(sol.y(4,1)),',','theta = ',num2str(theta)]);
    xlabel('r','FontSize',26,'FontName','Times');
    ylabel('z','FontSize',26,'FontName','Times');
    set(gca,'FontName','Times New Roman','FontSize',26,'LineWidth',2.5)
    set(figure(1), 'Color', 'white');
    drawnow;
    hold off; 
    areasphere=sol.y(9,end);

    if  (usefulC02-C02end)>0.01
        break;   
    end

end 
C02=usefulC02; 
'vary H'
h0=sol.parameters(1);
C02=fH(end,1);  
Hinitial=fH(end,2);  
sol.parameters=[h0,C02];
if  fH(end-1,2)<fH(end,2)
    varc = Hinitial:0.001:10;  
else
    varc = Hinitial:-0.001:0;  
end
varc = Hinitial:0.01:10;  
for H = varc  
    yeq = @(u,y,para) shape3(u,y,para,fSharp,ag1,ag2,C01,kappa0,kappa2,p,ap,theta);
    ybc = @(ya,yb,para)twobc3(ya,yb,para,r0,rb,psi0,C01,fSharp,ag1,ag2,kappa0,kappa2,H,p,ap,theta);
    opts = bvpset('RelTol',1e-6,'AbsTol',1e-10,'NMax',50000);
    soltemp = bvp5c(yeq,ybc,sol,opts);     
    if soltemp.stats.maxerr < RTol
           j=j+1;
          sol=soltemp;   
    figure(1);
    r = sol.y(3,:);
    z = sol.y(4,:);
    area11=sol.y(9,:);
    areaC0=ag2-ag1;
    position=(ag2+ag1)/2;
    Hprotrusion0=sol.y(4,1);
    usefulC02=sol.parameters(2);
    C02= usefulC02;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    dpsi= sol.y(2,:);
    h=sol.parameters(1);
    psi=sol.y(1,:);
    area=sol.y(9,:);
    u=sol.x;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalenergy=(1/2).*h.*p.*r.^2.*sin(psi)+(1/2).*h.*kappa0.*r.*((-1).*C01+dpsi.* ...
  h.^(-1)+r.^(-1).*sin(psi)).^2+(1/4).*ap.*h.*kappa2.*r.*((-1).*C02+ ...
  r.^(-1).*cos(theta).^2.*sin(psi)+dpsi.*h.^(-1).*sin(theta).^2) ...
  .^2.*(tanh(((-1).*ag1+area).*fSharp)+(-1).*tanh(((-1).*ag2+area).* ...
  fSharp));
bendenergy=(1/2).*h.*kappa0.*r.*((-1).*C01+dpsi.*h.^(-1)+r.^(-1).*sin(psi)) ...
  .^2+(1/4).*ap.*h.*kappa2.*r.*((-1).*C02+r.^(-1).*cos(theta).^2.* ...
  sin(psi)+dpsi.*h.^(-1).*sin(theta).^2).^2.*(tanh(((-1).*ag1+area) ...
  .*fSharp)+(-1).*tanh(((-1).*ag2+area).*fSharp));
pressureenergy=(1/2).*h.*p.*r.^2.*sin(psi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fH(j,:) = [usefulC02,Hprotrusion0,ag1,ag2,simps(u,totalenergy),simps(u,bendenergy),simps(u,pressureenergy),theta];
    dynamicCURVATURESOL=[dynamicCURVATURESOL sol];
    [~,idx1]=min(abs(area11-ag1));  
    [~,idx2]=min(abs(area11-ag2));   
    plot(r(1:idx1),z(1:idx1),'k',-r(1:idx1),z(1:idx1),'k','LineWidth',3);
    hold on
    plot(r(idx1:idx2),z(idx1:idx2),'g',-r(idx1:idx2),z(idx1:idx2),'g','LineWidth',3);
    hold on
    plot(r(idx2:end),z(idx2:end),'k',-r(idx2:end),z(idx2:end),'k','LineWidth',3);
    axis equal;  
    %ylim([0 6]);      
    xlim([-4 4]);

    title(['C0_{2} = ',num2str(usefulC02),',','areaC0 = ',num2str(areaC0),newline,'ag_{1} = ',num2str(ag1),',','ag_{2} = ',num2str(ag2),',','position = ',num2str(0.5*(ag1+ag2)),',','H = ',num2str(sol.y(4,1)),',','theta = ',num2str(theta)]);
    xlabel('r','FontSize',26,'FontName','Times');
    ylabel('z','FontSize',26,'FontName','Times');
    set(gca,'FontName','Times New Roman','FontSize',26,'LineWidth',2.5)
    set(figure(1), 'Color', 'white');
    drawnow;
    hold off; 
    areasphere=sol.y(9,end);

    if  (usefulC02-C02end)>0.01
        break;
    end
    
    else
        break;
    end

end


%save the data
 DIR = pwd; 
filename = [DIR,'/Initial_kappa0_',num2str(kappa0),...
    '_kappa2_',num2str(kappa2),...
    '_ap_',num2str(ap),...
    '_fSharp_',num2str(fSharp),...
    '_thetamax_',num2str(thetamax),...
    '_p_',num2str(p),...
    '_C02begin_',num2str(C02begin),...
    '_C02end_',num2str(C02end),...
    '_C02max_',num2str(C02),...
    '_ag1begin_',num2str(ag1start),...
    '_ag2begin_',num2str(ag2start),...
    '_nmax_',num2str(nmax),...
    '_rb_',num2str(rb),...
    '_RTol_',num2str(RTol),...
    '_date_',date,'.mat'];
sol_all = cell(1,2);  

    sol_all{1,1} = fH; 
    sol_all{1,2} =  dynamicCURVATURESOL;


    save(filename,'sol_all');


  