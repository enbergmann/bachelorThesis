% function plotGreyScale(x, c4n, n4e, folderName)
function plotGrayscale(c4n, n4e, x)
  % TODO maybe vectorizable, but this is very similar to plotCR now
  % on the other hand
    X1 = c4n(n4e(:, 1), 1);
    X2 = c4n(n4e(:, 2), 1);
    X3 = c4n(n4e(:, 3), 1);
    Y1 = c4n(n4e(:, 1), 2);
    Y2 = c4n(n4e(:, 2), 2);
    Y3 = c4n(n4e(:, 3), 2);
    X = [X1'; X2'; X3'];
    Y = [Y1'; Y2'; Y3'];
    
    axis off;
    axis equal;
    colormap gray;
    patch(X, Y, x', 'EdgeColor', 'none');
    drawnow;
    
    %TODO what happens here, is it relevant? think about it
    %F = getframe(reducedImg);
    %img = F.cdata;
    %imwrite(img, [folderName,'\reducedImg.png']);
end
