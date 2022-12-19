clear;
close all;
n=250;
l_FGF=300;
l_cell=5;
l_RA=10;
cell_type=zeros(n,1)+3;
phase=zeros(n,1);
psm=(1:110);
cell_type(psm)=2;
proliferation=2.28;
FGF_0=1;
FGF=(FGF_0*(exp(-(flip(1:psm(end)))*l_cell/l_FGF)));
w(psm)=(1.391-2.619*exp(-14.94*FGF));
t=0;
dt=1e-3;
time_since_division=0;
im=zeros(250,430,3);
cnt=1;
dif_func=zeros(250,1);
somitogenesis=zeros(250,1);
mesp2=zeros(250,1);
while t<200/4
    psm=find(cell_type==2);
    FGF=(FGF_0*(exp(-(flip(1:psm(end)))*l_cell/l_FGF)));
    if t<100/4
        FGF=(FGF_0*(exp(-(flip(1:psm(end)))*l_cell/l_FGF)));
    else
        FGF=(FGF_0*(exp(-(flip(1:psm(end)))*l_cell/l_FGF)));
    end
    RA=(100*(exp(-((psm-psm(1)))*l_cell/l_RA)))';
    dif_func(psm)=FGF(psm)./((RA+1));
    w(psm)=(1.391-2.619*exp(-14.94*FGF(psm)'));
    phase(psm)=phase(psm)+(w(psm)*dt)';
    
    new_somite=find(dif_func<0.08&cell_type==2);
    for j=1:length(new_somite)
        i=new_somite(j);
        if (mod(phase(i),2*pi)>0)&&(mod(phase(i),2*pi)<pi/2)
            somitogenesis(i)=1;
            mesp2(i)=mesp2(i)+dt;
        end
        if (mod(phase(i),2*pi)>=pi/2)&&somitogenesis(i)==1
            somitogenesis(i)=0;
            cell_type(i)=1;
        end
    end
    if time_since_division>1/proliferation
        cell_type(psm(end)+1)=2;
        phase(psm(end)+1)=phase(psm(end));
        time_since_division=0;
    end    
    if mod(cnt,100)==0
    im(:,cnt/100,1)=mesp2;
    im(psm,cnt/100,2)=(sin(phase(psm))+1)/2;        
    end
    cnt=cnt+1;
    t=t+dt;
    time_since_division=time_since_division+dt;
end
figure(1);
plot(linspace(0,50,500),movmean(sum(im(:,:,1))*5,2));hold on;
save('Control.mat','im');
figure(2);
subplot(1,2,1);
imshow(im);
title('Control');
axis square;
figure(3);
plot(linspace(1250,0,250),im(:,500,1))



clear;
n=250;
l_FGF=300;
l_cell=5;
l_RA=10;
cell_type=zeros(n,1)+3;
phase=zeros(n,1);
psm=(1:110);
cell_type(psm)=2;
proliferation=2.28;
FGF_0=1;
FGF=(FGF_0*(exp(-(flip(1:psm(end)))*l_cell/l_FGF)));
w(psm)=(1.391-2.619*exp(-14.94*FGF));
t=0;
dt=1e-3;
time_since_division=0;
im=zeros(250,430,3);
cnt=1;
dif_func=zeros(250,1);
somitogenesis=zeros(250,1);
mesp2=zeros(250,1);
while t<200/4
    psm=find(cell_type==2);
    FGF=(FGF_0*(exp(-(flip(1:psm(end)))*l_cell/l_FGF)));
    if t<100/4
        FGF=(FGF_0*(exp(-(flip(1:psm(end)))*l_cell/l_FGF)));
    else
        FGF=FGF_0*ones(psm(end),1)';
    end
    RA=(100*(exp(-((psm-psm(1)))*l_cell/l_RA)))';
    dif_func(psm)=FGF(psm)./((RA+1));
    w(psm)=(1.391-2.619*exp(-14.94*FGF(psm)'));
    phase(psm)=phase(psm)+(w(psm)*dt)';
    
    new_somite=find(dif_func<0.08&cell_type==2);
    for j=1:length(new_somite)
        i=new_somite(j);
        if (mod(phase(i),2*pi)>0)&&(mod(phase(i),2*pi)<pi/2)
            somitogenesis(i)=1;
            mesp2(i)=mesp2(i)+dt;
        end
        if (mod(phase(i),2*pi)>=pi/2)&&somitogenesis(i)==1
            somitogenesis(i)=0;
            cell_type(i)=1;
        end
    end
    if time_since_division>1/proliferation
        cell_type(psm(end)+1)=2;
        phase(psm(end)+1)=phase(psm(end));
        time_since_division=0;
    end    
    if mod(cnt,100)==0
    im(:,cnt/100,1)=mesp2;
    im(psm,cnt/100,2)=(sin(phase(psm))+1)/2;        
    end
    cnt=cnt+1;
    t=t+dt;
    time_since_division=time_since_division+dt;
end
figure(1);
hold on;
plot(linspace(0,50,500),movmean(sum(im(:,:,1))*5,2));hold on;
legend('Control','FGF');
ylabel('Somitic Mesoderm Length (mm)');
xlabel('Time (hr)')
axis square;
save('FGF.mat','im');
figure(2);
subplot(1,2,2);
imshow(im);
title('FGF');
axis square;
figure(3);
hold on;
plot(linspace(1250,0,250),im(:,500,1));
legend('Control','FGF');
axis([600 1250 0.6 1.6]);
ylabel('MESP2 Signal (a.u.)');
xlabel('P-A Position (\mum)')

axis square;