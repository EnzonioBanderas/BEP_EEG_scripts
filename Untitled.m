% nIter=1e5;
% time=zeros(nIter,1);
% for ii=1:nIter
%     tic

N=10;
a=100*rand(N);
check=a;

l=eye(N);
u=zeros(N);

for i=1:N-1
    if a(i,i)==0
        break
    else
        for k=(i+1):N
            a(k,i)=a(k,i)/a(i,i);
            for j=(i+1):N
                a(k,j)=a(k,j)-a(k,i)*a(i,j);
            end
        end
    end
end

u(triu(true(N),0))=a(triu(true(N),0));
l(tril(true(N),-1))=a(tril(true(N),-1));

sum(sum(l*u==check))
imagesc(l*u-check)
colorbar
% time(ii)=toc;

% end
% histogram(time)


      % Display the first 10 Zernike functions
      x = -1:0.01:1;
      [X,Y] = meshgrid(x,x);
      [theta,r] = cart2pol(X,Y);
      idx = r<=1;
      z = nan(size(X));
      n = [0  1  1  2  2  2  3  3  3  3];
      m = [0 -1  1 -2  0  2 -3 -1  1  3];
      Nplot = [4 10 12 16 18 20 22 24 26 28];
      y = zernfunIS(n,m,X(idx),Y(idx));
      figure('Units','normalized')
      for k = 1:10
          z(idx) = y(:,k);
          subplot(4,7,Nplot(k))
          pcolor(x,x,z), shading interp
          set(gca,'XTick',[],'YTick',[])
          axis square
          title(['Z_{' num2str(n(k)) '}^{' num2str(m(k)) '}'])
      end