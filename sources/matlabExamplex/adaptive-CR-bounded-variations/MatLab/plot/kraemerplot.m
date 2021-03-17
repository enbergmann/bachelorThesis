function kraemerplot(x, n4e, c4n, s4e, mu1, mu2, lambda, folderName)
% Plots the volume fraction of a given solution cf. CC, Bartels
% modified by Tran Ngoc Tien

% Initalisation
nrElems=size(n4e,1);
t1=sqrt(2*lambda*mu1/mu2);
t2=mu2/mu1*t1;

absDx4e=zeros(nrElems,1);

for j = 1 : nrElems
    nodes = n4e(j,:);   % nodes of this element
    coords = c4n(nodes,:); % coordinates for the nodes
    sides=s4e(j,:);
    grads = [coords';1 1 1]\[-2 0; 0 -2; 0 0]; % gradients for CR basis
    grads = grads([3 1 2],:); % reorder to fit DoF numbering
    xloc = x(sides);
    absDx4e(j) = sqrt(sum((xloc(1)*grads(1,:) + xloc(2)*grads(2,:) + xloc(3)*grads(3,:)).^2));
end

volumeFrac = figure;
plotVolumenfrac(c4n,n4e,absDx4e,t1,t2);
axis off;
saveas(volumeFrac, [folderName,'\volumeFrac.png'],'png');
img = imread([folderName,'\volumeFrac.png']);
image(img);

end

function plotVolumenfrac(c4n,n4e,t,t1,t2)

    X1 = c4n(n4e(:,1),1);
    Y1 = c4n(n4e(:,1),2);
    X2 = c4n(n4e(:,2),1);
    Y2 = c4n(n4e(:,2),2);
    X3 = c4n(n4e(:,3),1);
    Y3 = c4n(n4e(:,3),2);
    X = [X1'; X2'; X3'];
    Y = [Y1'; Y2'; Y3'];
    % Set the colour each triangle is filled with.
    Z=zeros(3,size(n4e,1));
 

for i=1:length(t)
    a=t(i);
    if (a<=t1)
    z(i)=0;
    elseif (a<t2 && a>t1)
    z(i)=(t2-a)/(t2-t1);
    elseif (a>=t2)
    z(i)=1;
    end
end

A=z;


patch(X,Y,Z,A,'EdgeColor','black');
drawnow;
end