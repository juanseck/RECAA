function [Poblacion] = generar_poblacion(SmartCells_no,lb,ub,dim)
%Genera una poblacion aleatoria de SmartCells_no smart-cells
%Cada smart-cell es un vector donde cada elementio está entre lb y ub de forma aleatoria, con tantos elementos como dim
if length(lb)==1
    Poblacion=(rand(SmartCells_no,dim)*(ub-lb))+lb;
else
    Poblacion=zeros(SmartCells_no,dim);
    for i=1:dim
        Poblacion(:,i)=(rand(SmartCells_no,1)*(ub(i)-lb(i)))+lb(i);
    end
end
end

