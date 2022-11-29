data = load('E:\BB Project\TetAlyze\generative.mat');  % please change the path here
figure(3)
s = sliceViewer(data);
% clc
% close all
% run('E:\BB data new\3. big1-1 Refeed experiment\N=3\B1868\24h\B1868 t=24 Poc1-mChy Re-feed-sld - 1009\Alignment.m');

% figname = fullfile('E:\BB data new\3. big1-1 Refeed experiment\N=3\B1868\8h\B1868 t=8 Poc1-mChy Re-feed-sld - 1010\Alignment.fig');
% if (exist(figname,'file'))
%     h1=openfig(figname, 'invisible');
% end
% ax = gca; 
% h = findobj(gca,'Type','line'); 
% 
% 
% X = [];
% Y = [];
% Z = [];
% for i=1:length(h)   
%     X = [X h(i).XData];
%     Y = [Y h(i).YData];
%     Z = [Z h(i).ZData];
% end
% 
% 
% data = vertcat(X,Y)';
% [coeff,score,latent] = pca(data);
% 
% 
% w1 = (max(score(1:end, 1))-min(score(1:end, 1)))/2;
% w2 = (max(score(1:end, 2))-min(score(1:end, 2)))/2;
% a = score(1:end, 1)/w1;
% b = score(1:end, 2)/w2;
% r2 = (score(1:end, 1)/w1).^2 + (score(1:end, 2)/w2).^2;
% % figure(2)
% figure(4)
% f=fit(Z'/max(Z), r2,'poly4')
% plot(f, Z'/max(Z), r2)
% % plot(r2, Z/max(Z), 'o')
% % plot3(score(1:end, 1), score(1:end, 2), Z,'o')
% % axis equal