function h = myOrthoviews(im,pos,tt,newfig,line_color,line_style,linewidth,varargin)
%% function Orthoviews(im,pos,tt,newfig)
%
% Display a orthoviews of a 3D image with gray colormap, removing the grid 
% and keeping the proportions of the image when displaying. 
% - im is the image to display
% - pos is a vector containing the position to extract planes in (i,j,k)
%   coordinates (i.e. (row,colunm,frame))
% - tt is for a global tilte
% - newfig is a boolean true to create a new figure (default true)
%
% -- Example
% load mri.mat; D=squeeze(D);
% Orthoviews(D);
% Orthoviews(D,[30,70,10]);
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
im = gather(im);
sz=size(im);%im=double(im)/double(max(im(:)));
minC = double(min(im(:))); 
maxC = double(max(im(:)));
%maxC = 1.46;
margin=10;
if nargin <2 || isempty(pos), pos=floor(sz/2); end
if nargin <3, tt=' '; end
if nargin <4, newfig=1; end
if nargin < 5, line_color = 'b';end
if nargin < 6, line_style = '--';end
if nargin < 7, linewidth = 1;end
if nargin > 9, invColor = varargin{1};else invColor = false; end
if newfig, figure; end
assert(length(sz)==3,'Wrong length of pos vector');
assert(all(pos>0) && all(pos<=sz),'Given position out of image range');
imm=ones(sz(1)+sz(3)+margin,sz(2)+sz(3)+margin)*0.93;
imm(1:sz(1),1:sz(2))=im(:,:,pos(3));
imm(sz(1)+margin+1:end,1:sz(2))=transpose(squeeze(im(pos(1),:,:)));
imm(1:sz(1),sz(2)+margin+1:end)=squeeze(im(:,pos(2),:));
imagesc(imm); axis image; axis off;

caxis([minC,maxC]);colorbar
if invColor
    colormap(flipud(gray));
else
    colormap gray;
end

line([pos(2) pos(2)],[1 sz(1)],'Color',line_color,'LineWidth',linewidth,'LineStyle',line_style);
line([1 sz(2)],[pos(1) pos(1)],'Color',line_color','LineWidth',linewidth,'LineStyle',line_style);
line([sz(2)+margin+pos(3) sz(2)+margin+pos(3)],[1 sz(1)],'Color',line_color,'LineWidth',linewidth,'LineStyle',line_style);
line([sz(2)+margin+1 sz(2)+margin+sz(3)+1],[pos(1) pos(1)],'Color',line_color,'LineWidth',linewidth,'LineStyle',line_style);
line([1 sz(2)],[sz(1)+margin+pos(3) sz(1)+margin+pos(3)],'Color',line_color,'LineWidth',linewidth,'LineStyle',line_style);
line([pos(2) pos(2)],[sz(1)+margin+1 sz(1)+margin+sz(3)],'Color',line_color,'LineWidth',linewidth,'LineStyle',line_style);
text(-5,pos(1),['I=',num2str(pos(1))],'FontSize',14,'HorizontalAlignment','right')
text(pos(2),-7,['J=',num2str(pos(2))],'FontSize',14,'HorizontalAlignment','center')
text(sz(2)+margin+pos(3),-7,['K=',num2str(pos(3))],'FontSize',14,'HorizontalAlignment','center')
text(-5,sz(1)+margin+pos(3),['K=',num2str(pos(3))],'FontSize',14,'HorizontalAlignment','right')

ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on','FontSize',14)
title(tt);
h = gcf;
end