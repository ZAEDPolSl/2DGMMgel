function scale = plot_gmm(x,y,data,gmm,opts)
% Plot 2D gel image versus GMM decomposition. 

if nargin < 5; opts = default_gmm_opts(); end
[m,n] = size(data);
ploty = zeros(m,n);
[x_vec,y_vec,stats] = par_vec(x,y);
len_row = stats.len_row; len_col = stats.len_col;
coord_vec = [x_vec',y_vec'];
alfa = gmm.alpha; center = gmm.center; covar = gmm.covar;
for c=1:gmm.KS
    tmp = alfa(c)*norm_pdf2(coord_vec,center(c,:),covar(:,:,c));
    ploty = ploty + reshape(tmp,len_row,len_col);
end
scale = sum(255-data(:))/sum(ploty(:));

if opts.show
%     datad = double(data);
%     sumD = sum(sum(255-datad));
%     ploty_norm = (255 - ploty/(sum(sum(ploty))) * (sumD));
%     ploty_norm = round(255-(ploty-min(min(ploty)))/(max(max(ploty))-min(min(ploty))) * max(max(datad)));
    figure; imshow(data,[],'Border','tight','InitialMagnification','fit')
    hold on; 
%     [~,lay1] = contour(scale*ploty,'Color','blue');% set(lay1,'AlphaData',0.9);
%     ch = get(lay1,'child'); alpha(ch,0.001)
    for a=1:gmm.KS
%         plot(gmm.center(a,2)-min(y)+1,gmm.center(a,1)-min(x)+1,'r.','MarkerSize',10)
        plot_2DGauss(gmm.center(a,:)-[min(x)-1,min(y)-1],rot90(gmm.covar(:,:,a),2));
    end
    
    
    % lay1 = imshow(ploty_norm); 
    % figure; contourf(flipud(ploty)); 
%     figure; subplot(1,2,1); title('Original image')
%     imshow(data,[],'Border','tight','InitialMagnification','fit')
%     subplot(1,2,2); title('2D Gaussian GMM model')
%     imshow(255-scale*ploty,[],'Border','tight','InitialMagnification','fit')
    
%     data2 = double(255-imadjust(data));
%     ploty_tmp = double(imadjust(uint8(255-ploty)));
%     data_diff = ploty_tmp - data2;
%     data_diff = scale*ploty - (255-double(data));
%     data_diff_3(:,:,1) = -min(data_diff,0);
%     data_diff_3(:,:,2) = 0*ones(size(data_diff));
%     data_diff_3(:,:,3) = max(data_diff,0);
%     data_diff_3 = uint8(data_diff_3);
%     figure; imshow(data_diff_3,'Border','tight','InitialMagnification','fit')
%     title('Red-lower, Blue-higher')
%     Z = imabsdiff(scale*ploty,double(255-data));
%     figure; imshow(Z,[])
%     imshowpair(scale*ploty,double(255-data));
    %     figure; subplot(1,2,1);imshow(255-data_diff_3(:,:,1),[]); subplot(1,2,2);imshow(255-data_diff_3(:,:,3),[]); 
end

% function y = norm_pdf(x,center,covar)
% den = (6.283185307179585 * sqrt(norm(covar)));
% y = exp(-0.5*(((x - center)/covar)*(x - center)'))/den;

function y = norm_pdf2(x,center,covar)
den = (6.283185307179585 * sqrt(det(covar)));
x_centr = bsxfun(@minus,x,center);
% x_centr = (x-repmat(center,length(x),1));
y = exp(-0.5*sum((x_centr/covar).*x_centr,2))/den;
y = y';