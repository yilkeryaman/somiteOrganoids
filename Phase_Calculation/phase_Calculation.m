I=imread('.\Kymograph.tif');
sz=size(I);
t=(1:sz(2))*0.25;
cnt=1;
for i=1:sz(1)
    y=I(i,:);
    [p]=sineFit(t,smooth(detrend(double(y)))');
    if p(3)>0.21 && p(3)<0.24
        data(stack).signal(cnt,:)=double(y);
        data(stack).position(cnt)=i;
        cnt=cnt+1;
    end
end
signal=[];
phase=[];
sz=size(data(stack).signal);
for i=1:sz(1)
    temp=data(stack).signal(i,:);
    temp=temp-movmean(temp,24);
    temp=temp./movstd(temp,24);
    data(stack).detrended_signal(i,:)=temp;
end
for i=1:sz(1)
    data(stack).detrended_signal(i,:)=smooth(data(stack).detrended_signal(i,:));
end
for k=1:sz(2)
    signal(:,k)=movmean(data(stack).detrended_signal(:,k),10);
end

for i=1:sz(1)
    y=signal(i,:);
    y=y-movmean(y,24);
    y=y./movstd(y,24);
    phase(i,:)=unwrap(angle(hilbert(y(:,1:end))));
end
data(stack).phase=phase;
close all;

range=152:172;
sz=size(data(stack).phase);
y=unwrap(mean(data(stack).phase(:,range),2));
x=linspace(0,1,sz(1));
xx=linspace(0,1,300);
phase=spline(x,y,xx);
phase=phase-mean(phase);
hold on;
plot(linspace(0,1,300),phase);
xlabel('Relative P-A Position (a.u.)');
ylabel('Phase (rad)');
axis square;