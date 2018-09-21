function [s,m] = fe_matrices(c4n,n4e)
  % return s: Steifigkeitsmatrix
  % return m: Massematrix

nC = size(c4n,1); % Anzahl Knoten
nE = size(n4e,1); % Anzahl Elemente

  % Massematrix eines Dreiecks mit Fläche 1 um in der folgenden Schleife nur
  % das Volumen des betrachteten Dreicks multiplizieren zu müssen.
m_loc = (ones(3)+eye(3))/12; 

  % (counter, counter_max) Hilfsvariablen zum Erzeugen von I, J, X_s und X_m
ctr = 0;
ctr_max = 9*nE;

  % Initialisierung der Hilfsmatrizen zum Erstellen der Masse- und Steifig-
  % keitsmatrizen.

  % I und J zum identifizieren der lokalen Knotennummern eines Dreicks mit
  % den globalen Knotennummern der Triangulierung.
I = zeros(ctr_max,1); 
J = zeros(ctr_max,1);  
  % X_s zum Speichern der Einträge der lokalen Steifigkeitsmatrix eines
  % Dreicks.
X_s = zeros(ctr_max,1); 
  % X_m zum Speichern der Einträge der lokalen Massematrix eines
  % Dreicks.
X_m = zeros(ctr_max,1);
  % Vektor mit den Flächeninhalten der einzelnen Dreiecke der Triangulierung.
vol_T = zeros(nE,1);

for j = 1:nE % Das j-te Dreieck T=conv{P_1,P_2,P_3} wird betrachtet.
    
      % Hilfsmatrix zum Berechnen der Gradienten in grads_T und des Volumens
      % des Dreiecks.
    X_T = [ones(1,3);c4n(n4e(j,:),:)'];
      % Für k=1,2,3 stehen in der k-ten Zeile von grads_T die beiden 
      % Komponenten des Gradienten der nodalen Basisfunktion zur Ecke P_k.
    grads_T = X_T\[zeros(1,2);eye(2)];

      % Das Volumen von T wird berechnet und gespeichert.
    vol_T(j) = det(X_T)/2;
    
      % Berechnen der Einträge der Hilfsmatrizen.
    for m = 1:3
        for n = 1:3
            ctr = ctr+1; 
            
            I(ctr) = n4e(j,m); 
            J(ctr) = n4e(j,n);
            
            X_s(ctr) = vol_T(j)*grads_T(m,:)*grads_T(n,:)';
            X_m(ctr) = vol_T(j)*m_loc(m,n);
        end
    end
end

s = sparse(I,J,X_s,nC,nC); % Steifigkeitsmatrix
m = sparse(I,J,X_m,nC,nC); % Massematrix