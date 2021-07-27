function ms = mixed_matrix(c4n,n4e)

nC = size(c4n,1); % Anzahl Knoten
nE = size(n4e,1); % Anzahl Elemente

  % (counter, counter_max) Hilfsvariable zum Erstellen von I, J und X.
ctr = 0; 
ctr_max = 6*nE;

  % Initialisieren der Hilfsmatrizen, die benötigt werden um die gemischte
  % Matrix zu erstellen.
  
  % I und J Matrizen, um die lokalen Knotennummern eines Dreiecks mit
  % den globalen Knotennummern der Triangulierung identifizieren und die 
  % globale gemischte Matrix zu erstellen.
I = zeros(ctr_max,1); 
J = zeros(ctr_max,1); 
  % Gemischte Matrix.
X = zeros(ctr_max,1);

for j = 1:nE % 
    X_T = [ones(1,3);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,2);eye(2)];
    vol_T = det(X_T)/2;
    for k = 1:3
        for ell = 1:2
            ctr = ctr+1;
            I(ctr) = n4e(j,k); 
            J(ctr) = (j-1)*2+ell;
            X(ctr) = vol_T*grads_T(k,ell);
        end
    end
end
ms = sparse(I,J,X,nC,2*nE);
