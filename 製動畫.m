%https://octave-online.net/
clear;clc;close all;
color=['r-o','g-o','b-o','k-o','y-o','m-o','c-o','r-s','g-s','b-s','y-s','m-s','c-s'];
N=21; num=14;
data=load('moive_data.txt');
image=zeros(N,N);
bound=1:N;
for t=1:num
	%figure(t);
	for i=1:N
	    image(N+1-i,:)=data((t-1)*N+i,:);
	end
	[c,handle] = contour(image, 20);
	clabel(c, handle);
    M(t)=getframe;
end
movie(M, 1);