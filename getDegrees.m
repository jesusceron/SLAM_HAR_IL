function [ Positions_final, error_final ] = getDegrees(StrideLengths, Thetas, initial_coordinates, final_coordinates)

Positions = [];
Positions(1,:) = initial_coordinates;
%Hallo la diferencia entre el angulo real y el hallado en estas
%actividades
    
    for iDegree=1:360
        
        for iStride=1:size(StrideLengths,1)
            
            if iStride==1
                Positions(iStride+1,1)=StrideLengths(iStride)*cos(Thetas(iStride)+iDegree*(pi/180)) + initial_coordinates(1,1); % X
                Positions(iStride+1,2)=StrideLengths(iStride)*sin(Thetas(iStride)+iDegree*(pi/180)) + initial_coordinates(1,2); % Y
                
                if iStride == 1
                    
                    error_aux(iDegree) = sqrt(((final_coordinates(1)-Positions(iStride+1,1))^2)+((final_coordinates(2)-Positions(iStride+1,2))^2));
                    
                    if iDegree>1
                        if error_aux(iDegree)<error_final
                            error_final = error_aux(iDegree);
                            Positions_final = Positions;
                        end
                    else
                        error_final  = error_aux(iDegree);
                        Positions_final = Positions;
                    end
                    
                end
                    
                    
                
            else
                % acumulo con SL y thetas
                Positions(iStride+1,1)=Positions(iStride+1-1,1) + StrideLengths(iStride)*cos(Thetas(iStride)+iDegree*(pi/180)); % X
                Positions(iStride+1,2)=Positions(iStride+1-1,2) + StrideLengths(iStride)*sin(Thetas(iStride)+iDegree*(pi/180)); % Y
                
                %calcular error de la trayectoria de cada actividad (coordenadas finales menos las verdaderas finales)
                if iStride==size(StrideLengths,1)
                    
                    error_aux(iDegree) = sqrt(((final_coordinates(1)-Positions(iStride+1,1))^2)+((final_coordinates(2)-Positions(iStride+1,2))^2));
                    
                    if iDegree>1
                        if error_aux(iDegree)<error_final
                            error_final = error_aux(iDegree);
                            Positions_final = Positions;
                        end
                    else
                        error_final  = error_aux(iDegree);
                        Positions_final = Positions;
                    end
                    
                end
                
            end
        end
    end
end