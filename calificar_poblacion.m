function [calificacion] = calificar_poblacion(poblacion,fcosto)
num_sol=size(poblacion,1);
calificacion=zeros(1,num_sol);
for i=1:num_sol
    calificacion(i)=fcosto(poblacion(i,:));
end
end

