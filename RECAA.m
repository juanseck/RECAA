%  Reversible Elementary Cellular Automata Algorithm (RECAA)
%
%  Source codes version 1.0
%
%  Developed in MATLAB R2015a(7.08)
%
%  Author and programmer: Juan Carlos Seck Tuoh Mora
%
%       email:   jseck@uaeh.edu.mx
%                juanseck@gmail.com
%
%       Homepage:
%
%  Main paper: A new algorithm inspired on reversible elementary cellular automata for global optimization 
%  Autores: Juan Carlos Seck-Tuoh-Mora, Omar Lopez-Arias, Norberto Hernandez-Romero, 
%  Genaro J. Martinez, Valeria Volpi-Leon
%  IEEE Access, DOI: http://
%_______________________________________________________________________________________________
% You can simply define your cost function in a seperate file and load its handle to fobj
% The initial parameters that you need are:
%__________________________________________
% SmartCells_no = number of smart-cells
% Neighbors_no = number of neighbors
% Max_iteration = maximum number of iterations
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% dim = number of variables to be tunned
% fobj = @YourCostFunction
%
% If all the variables have equal bounds you can just define lb and ub as two single numbers
%
% To run RCAA: [min_value,position_vector,convergence_curve]= RCAA(SmartCells_no,Neighbors_no,Max_iteration,lb,ub,dim,fobj)
%
% Example: 
% [lb,ub,dim,fobj] = benchmark_functions('F10');
% dim = 500;
% [min_value,position_vector,convergence_curve]=RCAA(12,6,500,lb,ub,dim,fobj)
%______________________________________________________________________________________________


function [min_value,position_vector,convergence_curve,smart_cells] = RECAA(SmartCells_no,Neighbors_no,Max_iteration,lb,ub,dim,fcosto,bandera_imp)

%Parámetros del algoritmo     
proporcion=1.5;       %1   1.7
num_dec_inf=1;      %1  2
num_dec_sup=6;      %4  6

%Generate random smart-cells
[smart_cells] = generar_poblacion(SmartCells_no,lb,ub,dim);
new_smart_cells=smart_cells;

%Evalute fitness of smart-cells
[calificacion] = calificar_poblacion(smart_cells,fcosto);
%[calificacion] = calificar_poblacion(smart_cells,fcosto,num_fc);
new_calificacion = calificacion;
convergence_curve=zeros(1,Max_iteration);
convergence_curve(1)=min(calificacion);

%Número de soluciones elitistas
elitist_num=2;


%Optimization cycle
for i=2:Max_iteration
    %Imprime numero de iteracion
    if bandera_imp>0
        if mod(i,bandera_imp)==0
            disp(['RECAA Iteracion: ' num2str(i)]);
        end
    end
    %Elitism, dejar el mejor en la nueva poblacion
    [~,indE]=mink(calificacion,elitist_num);
    new_smart_cells(1:elitist_num,:)=smart_cells(indE,:);
    new_calificacion(1:elitist_num)=calificacion(indE);
    %Mejor calificacion
    mejor_calif=new_calificacion(1);
    %Tomara todas las smart-cells
    for j=elitist_num+1:SmartCells_no
        %Smart-cell actual
        smart=smart_cells(j,:);
        calif=calificacion(j);
        %Vecino 1, smart-cell aleatorio
        indA=randi(SmartCells_no);
        vecino1=smart_cells(indA,:);
        calif1=calificacion(indA);
        %Mejor de cada vecindad
        nueva_calif=inf;        
        %Neighborhoods of non-elitist solutions
        for k=1:Neighbors_no
            %Elegir regla de forma aleatoria
            regla_ale=rand();
            if regla_ale <=1/12
                evolucion = regla_redondeo(smart,calif,mejor_calif,randi([num_dec_inf,num_dec_sup]));
            elseif regla_ale>1/12 && regla_ale<0.2 %0.25
                evolucion = regla_corrimiento(smart,proporcion,lb,ub);
            else
                evolucion = regla_corrimiento_vecino(smart,calif,vecino1,calif1,proporcion,lb,ub);
            end
            %Pasar valores de evolucion a intervalo valido
            ind_Neg=evolucion<lb;
            if sum(ind_Neg)>0
                if length(lb)>1
                    evolucion(ind_Neg)=lb(ind_Neg)+(((ub(ind_Neg)-lb(ind_Neg))/4)*rand);
                else
                    evolucion(ind_Neg)=lb+((ub-lb)/4)*rand;
                end
            end
            ind_Pos=evolucion>ub;
            if sum(ind_Pos)>0
                if length(ub)>1
                    evolucion(ind_Pos)=ub(ind_Pos)-(((ub(ind_Pos)-lb(ind_Pos))/4)*rand);
                else
                    evolucion(ind_Pos)=ub-((ub-lb)/4)*rand;
                end
            end
            %Calificar vecino
            calif_evolucion = fcosto(evolucion);
            %Best of every neighborhood
            if k==1 || calif_evolucion<nueva_calif                
                mejor_vecino=evolucion;
                nueva_calif=calif_evolucion;
            end
        end
        %Ver si se cambia a la smart cell apara conservar diversidad
        if rand < 0.05 || nueva_calif<calif
            smart=mejor_vecino;
            calif=nueva_calif;
        end            
        %Agregar a new_smart_cells
        new_smart_cells(j,:)=smart;
        new_calificacion(j)=calif;
    end
    %Pasar la nueva poblacion y su calificacion a la actual
    smart_cells=new_smart_cells;
    calificacion=new_calificacion;
    convergence_curve(i)=min(calificacion);
end

%Return the best obtained solution
[min_value,indE]=min(calificacion);
position_vector=smart_cells(indE,:);

end

%Reglas

%Reglas de explotación, usan solo al información de la misma smart-cell

%Regla para redondeo en los valores de la smart_cell
function [evolucion] = regla_redondeo(smart,calif,mejor_calif,num_dec)
evolucion=smart;
%suma de calificaciones
suma=calif+mejor_calif;
%ponderacion de calificación del vecino
pond=1-(calif/suma);
%Hacer redondeo
for i=1:length(smart)
    if rand<=pond
        evolucion(i)=round(evolucion(i),num_dec);
    end
end
end

%Regla para simular una regla reversible ECA(2,1) tomando la posicion a la
%izquierda, la misma o a la derecha de cada elemento de la smart-cell como
%factor de cambio, usando frontera periódica
function [evolucion] = regla_corrimiento(smart,proporcion,lb,up)
evolucion=smart;
%Elegir aleatoriamente corrimiento
tipo_corr=randi([-1,1]);
%tipo_corr=randi([-7,7]);
%Elegir si se toma valor original o complemento 
if tipo_corr==0 
    tipo_com=1;
else
    tipo_com=randi([0,1]);
end
%Hacer corrimiento para tomar vecino adecuado
corrimiento=circshift(smart,tipo_corr);
%Tomar complemento si es necesario
if tipo_com==1
    aux=corrimiento-lb;
    corrimiento=up-aux;
end
%nueva distancia
r=(rand*proporcion);
dist=(smart-corrimiento)*r;
evolucion=evolucion-dist;
end

%Reglas para exploración, usa información de la smart-cell y un vecino

%Regla para simular una regla reversible ECA(2,1) tomando la posicion a la
%izquierda, la misma o a la derecha de un vecino de la smart-cell como
%factor de cambio, usando frontera periódica
function [evolucion] = regla_corrimiento_vecino(smart,cal_smart,vecino,cal_vecino,proporcion,lb,up)
evolucion=smart;
%Elegir aleatoriamente corrimiento
tipo_corr=randi([-1,1]);
%Elegir si se toma valor original o complemento 
tipo_com=randi([0,1]);
%Hacer corrimiento para tomar vecino adecuado
corrimiento=circshift(vecino,tipo_corr);
%Tomar complemento si es necesario
if tipo_com==1
    aux=corrimiento-lb;
    corrimiento=up-aux;
end
%suma de calificaciones
suma=cal_smart+cal_vecino;
%ponderacion de calificación del vecino, si es mucho mejor, pond será
%grande
pond=1-(cal_vecino/suma);
%Si se selecciona al vecino, se suma al elemento del vector multiplicado
%por un aletorio r entre -dist_mayor/2 y dist_mayor/2
r=(rand*proporcion)-(proporcion/2);
for i=1:length(smart)
    if rand<=pond
        evolucion(i)=evolucion(i)+(r*corrimiento(i));
    end
end
end





