function plotGreyScale(x, c4n, n4e, folderName)
    
    reducedImg = figure;
    axis off;
    axis equal;
    
    plotVolumenfrac(x,c4n,n4e);
    
    F = getframe(reducedImg);
    img = F.cdata;
    imwrite(img, [folderName,'\reducedImg.png']);
end

function plotVolumenfrac(x, c4n, n4e)

    X1 = c4n(n4e(:,1),1);
    Y1 = c4n(n4e(:,1),2);
    X2 = c4n(n4e(:,2),1);
    Y2 = c4n(n4e(:,2),2);
    X3 = c4n(n4e(:,3),1);
    Y3 = c4n(n4e(:,3),2);
    X = [X1'; X2'; X3'];
    Y = [Y1'; Y2'; Y3'];
    Z=zeros(3,size(n4e,1));
    
    colormap gray;
    patch(X,Y,Z,x','EdgeColor','none');
    drawnow;
end