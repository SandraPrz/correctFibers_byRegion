
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AUTHORSHIP -------------------------------------------------------- %%%
%                  Author:     Sandra Perez Herrero
%           Creation date:     28/11/2024
%       Last modification:     28/11/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Este código tiene como objetivo modificar de forma aislada las fibras en
% una determinada zona de la malla etiquetada. Lo que hacemos es usar un
% vector definido según su punto incial - final (fibre_data_vectors) y la
% dirección en la que queremos poner las fibras (fibre_data_align =
% transversal/longitudinal con respecto al vector introducido y la normal
% de la superficie). El código pone las fibras tanto en el endocardio como
% en el epicardio de la región deseada y una vez etiquetadas estas nuevas
% fibras en estas mallas de superficie, para cada centroide de elemento de
% la malla se le mapea la fibra más cercana asignada en alguna de estas
% superficies.
% PREPARE THE INPUTS FOR FIBERS MODIFICATION
% 0) Definir el vector en paraview: sources > line > introducir point1 y point2
% 1) malla volumetrica total con las etiquetas necesarias para Elvira
% 2) paraview > cargar la malla .vtk en hexaedros ya etiquetada > threshold > RegionID (donde quieres cambiar las fibras, también puedes hacer un clip activando crinkle clip siempre) ---> guardar como .vtk (esto es la malla volumétrica de la region)
% 3) a partir de la anterior > threshold (surfnodes = 1) > filters > extract surface > filters > triangulate > connectivity > threshold > guardar por separado epi y endo como .vtk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% STEP2: CHANGE THE FIBERS AT THE REGION
% Read surface meshes
epiPathVtk='C:\Users\Sandra\Desktop\laura_fibers\epi.vtk'; % malla de superficie correspondiente al epicardio de la zona que queremos modificar
endoPathVtk='C:\Users\Sandra\Desktop\laura_fibers\endo.vtk';% malla de superficie correspondiente al endocardio de la zona que queremos modificar

% Read voluetric meshes
fullPathInriaMesh='D:\RecursosModelado3D\SANDRA\03_ME_LO_HAN_PASADO\auricula\RECONSTRUCCION_AURICULAR_CI2B\RECONSTRUCCION_AURICULAR_CI2B\Mallas\ps6691641_HEXA.vtk'; % malla completa en hexas con todas las etiquetas
regionPathVtk='C:\Users\Sandra\Desktop\laura_fibers\regionID_28.vtk';% malla volumétrica con los hexahedros que contienen las fibras que quremos modificar. Esta malla proviene de la total y tiene todas las etiquetas que luego usamos para Elvira.

fullPathResultado='C:\Users\Sandra\Desktop\laura_fibers\'; %path para guardar los resultados

%dirección del vector de paraview
fibre_data_vectors = [8,	8, 8,	-4.90147693947432,	-11.1722269500000,	85.4954385345663];
%Cargar la alineación de las fibras
fibre_data_align='transversal' %puede tomar los valores transversal o longitudinal (direcciones respecto al vector introducido)


correctFibers_byRegion_f(epiPathVtk,endoPathVtk,regionPathVtk,fullPathInriaMesh,fullPathResultado,fibre_data_vectors,fibre_data_align)


%%

function correctFibers_byRegion_f(epiPathVtk,endoPathVtk,regionPathVtk,fullPathInriaMesh,fullPathResultado,fibre_data_vectors,fibre_data_align)

% almacenar la información de las mallas de superficie en una estructura
% tipo mesh
mesh_epi = vtk2meshReader_f(epiPathVtk,'tri');
mesh_endo = vtk2meshReader_f(endoPathVtk,'tri');



%%  Calcular "Fibre" endocardios y epicardio

%asignacion de fibras
mesh_epi=assign_Fibre_f(mesh_epi,fibre_data_vectors,fibre_data_align);
mesh_endo=assign_Fibre_f(mesh_endo,fibre_data_vectors,fibre_data_align);


%Convertimos a vtk para asegurarnos de que se etiquetan bien
mesh2vtkWriter_f(mesh_epi,epiPathVtk,'epi_fibers.vtk');
mesh2vtkWriter_f(mesh_endo,endoPathVtk,'endo_fibers.vtk');

%% 3.3. Transferir Fibers desde las skin

region_mesh_hex=vtk2meshReader_f(regionPathVtk,'hex');

%create Fiber in cells

region_mesh_hex = transfer_skin_hex_f(region_mesh_hex,mesh_epi,mesh_endo,"Fiber",'cell');
mesh2vtkWriter_f(region_mesh_hex,fullPathResultado,'region_mesh_hex.vtk');



region_mesh_hex=vtkToolbox2Ci2b_vtkReader(strcat(fullPathResultado,'region_mesh_hex.vtk'),1);
region_mesh_hex_original=vtkToolbox2Ci2b_vtkReader(regionPathVtk,1);

valIdx_region_bad   = find(strcmp({region_mesh_hex.CellData.Name}, 'Fiber'));
valIdx_region_good   = find(strcmp({region_mesh_hex_original.CellData.Name}, 'Fiber'));
region_mesh_hex_original.CellData(valIdx_region_good).Data=region_mesh_hex.CellData(valIdx_region_bad).Data;
%vtkToolbox2Ci2b_vtkWriter(region_mesh_hex_original,fullPathResultado,'region_mesh_hex_original.vtk');


%% 3.4 Transferir a la malla completa solo las nuevas fibras
mesh_hex=vtkToolbox2Ci2b_vtkReader(fullPathInriaMesh,1);

%%
% detectar posiciones de las variables por nombre
idIdx    = find(strcmp({region_mesh_hex_original.CellData.Name}, 'IDRef'));
valIdx   = find(strcmp({region_mesh_hex_original.CellData.Name}, 'Fiber'));
targetIdx= find(strcmp({mesh_hex.CellData.Name}, 'IDRef'));
outIdx   = find(strcmp({mesh_hex.CellData.Name}, 'Fiber'));

% arrays fuente
keys  = region_mesh_hex_original.CellData(idIdx).Data(:);   % ids donde modificar las fibras
vals  = region_mesh_hex_original.CellData(valIdx).Data(:,:);  % valores que quieres asignar



% escribir de vuelta
transform(:,1)=keys;
transform(:,2:4)=vals;

for i=1:size(transform,1)
    mesh_hex.CellData(outIdx).Data(transform(i,1),:) = transform(i,2:4);
end

%%
vtkToolbox2Ci2b_vtkWriter(mesh_hex,fullPathResultado,'atria_final.vtk')


end



function mesh = vtk2meshReader_f(fullPath,vtktype)

    %transform VTK 5.1 file to struct
    %struct=vtk2structReader(fullPath);
    struct=vtkToolbox2Ci2b_vtkReader(fullPath,1);
    
    %transform struct to mesh
    mesh=struct2meshReader(struct,vtktype);
    
end

function mesh2vtkWriter_f(mesh,fullPath,fileName)

    %mesh to struct
    struct=mesh2structWriter(mesh);
    
    %write struct to vtk
    vtkToolbox2Ci2b_vtkWriter(struct,fullPath,fileName);
end


function mesh = assign_Fibre_f(mesh,fibres_data,fibres_alignment)
    
    %obtain the values of "OrganID" for cells
    organid_data=get_data_values_f(mesh,"OrganID",'cell');
    
    %delete duplicated values
    organid_values=unique(organid_data,'rows');
    
    %sort values from lower to higher
    organid_values=sortrows(organid_values,1);
    
    %find which is the index of "Normal" for cells
    normal_data=get_data_values_f(mesh,"Normal",'cell');

    for i=1:length(organid_values)
        
        %find the cells whose OrganID value is equal to the currently
        %examined OrganID value
        target_cells=mesh.cells(organid_data==organid_values(i),1);
        
        %obtain the values of the director vector based on "fiber_data"
        director_vector_beg=fibres_data(1:3);
        director_vector_end=fibres_data(4:6);
        director_vector=director_vector_end-director_vector_beg;
        
        %obtain value of "Normal" of "target_cells"
        normal_vectors=normal_data(target_cells,:);
        
        %obtain type of alignment
        alignment=fibres_alignment;
        
        %check if normal is a null vector
        if isequal(director_vector,[0 0 0]) %if it is
            
            %save fiber as a null vector
            fibres=repmat([0 0 0],size(normal_vectors,1),1);
            
        else %if it is not
            
            %create a repetition matrix with the values of 
            %"direction_vector" in each row and with as many rows as 
            %"target_cells" are
            director_vector=repmat(director_vector,size(normal_vectors,1),1);
            
            %depending on the type of alignment
            switch alignment
                
                case "longitudinal"
                    
                    %get the fiber as the cross product of "normal_vector" 
                    %times "director_vector"
                    fibres=cross(normal_vectors,director_vector,2);
                    
                    %get the cross product of "normal_vectors" times
                    %"fibres" so the new "fibres" value represent vectors
                    %that are parallel to the surface of the cell
                    fibres=cross(normal_vectors,fibres,2);
                    
                    %get the norm of each fibres value
                    norms=cellfun(@norm,num2cell(fibres,2));
            
                    %make each fibre unitary
                    fibres=fibres./norms;
    
                case "transversal"
                    
                    %get the fiber as the cross product of 
                    %"director_vector" times "normal_vector"
                    fibres=cross(director_vector,normal_vectors,2);
                    
                    %get the norm of each fibres value
                    norms=cellfun(@norm,num2cell(fibres,2));
        
                    %make each fibre unitary
                    fibres=fibres./norms;
    
                otherwise
                    
                    %raise an error
                    errordlg(strcat("alignment value for OrganID ", string(organid_values(i)), " is not 'longitudinal' nor 'transversal'"),'assign_fibres()');
                    
                    %save fiber as a null vector
                    fibres=repmat([0 0 0],size(normal_vectors,1),1);
                    
            end
            
        end
        
        %write the corresponding OrganID value to the vertices
        mesh=write_data_f(mesh,target_cells,fibres,"Fiber",'cell');
    end
end

function values = get_data_values_f(mesh,data_name,pointcell)

    %depending on the values of "pointcell"
    switch pointcell
        
        case 'point'
            %get the information describing the position of the data in
            %"mesh"
            [data_index,data_width]=get_data_info_f(mesh,data_name,'point');
            
            %search the information and store it in values
            values=mesh.points(:,(data_index-1)*10+1:(data_index-1)*10+data_width);
            
            
        case 'cell'
            %get the information describing the position of the data in
            %"mesh"
            [data_index,data_width]=get_data_info_f(mesh,data_name,'cell');
            
            %search the information and store it in values
            values=mesh.cells(:,(data_index-1)*10+1:(data_index-1)*10+data_width);
            
        otherwise
            
            %raise an error
            errordlg('pointcell values must be cell or point','get_data_values()');
                
            %save "values" as empty
            values=[];
                
            %end the function
            return;
    end
end
    

function [data_index,data_width] = get_data_info_f(mesh,name,pointcell)

    %depending on the value of "pointcell"
    switch pointcell
        
        case 'cell'
            
            %find which is the index of "name" in "cell_data"
            data_index = find(mesh.cell_data_text(:,2) == name);
    
            %check if "name" was found
            if isempty(data_index) %if it was not
                
                %raise an error
                errordlg('The specified data structure does not exist in cells','get_data_info()');
                
                %save "data_with" as empty
                data_width=[];
                
                %end the function
                return;
                
            else %if it was
                
                %obtain "data_width" from "cell_data_width"
                data_width = mesh.cell_data_width(data_index);
            end
        
        
        case 'point'
            %find which is the index of  "name" in "point_data"
            data_index = find(mesh.point_data_text(:,2)==name);
    
            %check if "name" was found
            if isempty(data_index) %if it was not
                
                %raise an error
                errordlg('The specified data structure does not exist in points','get_data_info()');
                
                %save "data_with" as empty
                data_width=[];
                
                %end the function
                return;
                
            else %if it was
                
                %obtain "data_width" from "point_data_width"
                data_width=mesh.point_data_width(data_index);
            end
            
            
        otherwise
            
            %raise an error
            errordlg('pointcell value must be cell or point','get_data_info()');
                
            %save "data_index" as empty
            data_index=[];
            
            %save "data_with" as empty
            data_width=[];
                
            %end the function
            return;
    end
end
    

function mesh=write_data_f(mesh,target_id,data,data_name,pointcell)
    
    %switch according to the value of pointcell
    switch pointcell
        
        case 'cell'
            %find which is the index of "data_name" and its width
            [data_index,data_width]=get_data_info_f(mesh,data_name,'cell');
    
            %if "data_structure" cannot be found, end the function
            if isempty(data_index)
        
                %raise an error
                errordlg('The specified data structure does not exist','write_data()');
        
                %end function
                return
            end
    
            %assign to the target_ids their corresponding value
            mesh.cells(target_id,(data_index-1)*10+1:(data_index-1)*10+data_width)=data;
            
        case 'point'
            %find which is the index of "data_name" and its width
            [data_index,data_width]=get_data_info_f(mesh,data_name,'point');
    
            %if "data_structure" cannot be found, end the function
            if isempty(data_index)
                
                %raise an error
                errordlg('The specified data structure does not exist','write_data()');
        
                %end function
                return
            end
    
            %assign to the target_ids their corresponding value
            mesh.points(target_id,(data_index-1)*10+1:(data_index-1)*10+data_width)=data;
            
        otherwise
            %raise an error
            errordlg('value of pointcell must be cell or point','write_data()');
    end
end

function mesh_hex = transfer_skin_hex_f(mesh_hex,mesh_skin_epi,mesh_skin_right,data_name,pointcell)

    %switch depending on pointcell
    switch pointcell
        case 'point'
            %% Find Closest Point to HexMesh in Skin Meshes%%
            %create 'skin_points' matrix with an ascending vector form 1
            %to the sum of the number of points in all 3 skin meshes in 
            %the first column and with the coordinates of the points of 
            %all three skin meshes in the second, third and fourth columns
            skin_points(:,1) = [1:sum([mesh_skin_epi.num_points,mesh_skin_right.num_points])]';
            skin_points(:,2:4) = [mesh_skin_epi.points(:,2:4);mesh_skin_right.points(:,2:4)];
    
            %get the closest point in the skin meshes to each centroid in 
            %"mesh_hex"
            [closest_point,~,~] = find_closest_point(mesh_hex.points,skin_points);
    
            %create a matrix that contains in the first column an integer 
            %refering to the skin mesh that is the closest to the 
            %"mesh_hex" point (1 for epi, 2 for right) 
            %and in the second column its adjusted closest cell value
            closest_point_mesh(closest_point(:,1)<=mesh_skin_epi.num_points,1) = 1;
            closest_point_mesh(closest_point_mesh(:,1)==1,2) = closest_point(closest_point_mesh(:,1)==1,1);
    
            closest_point_mesh((closest_point(:,1)>mesh_skin_epi.num_points & closest_point(:,1)<=(mesh_skin_epi.num_points+mesh_skin_right.num_points)),1) = 2;
            closest_point_mesh(closest_point_mesh(:,1)==2,2) = closest_point(closest_point_mesh(:,1)==2,1)-double(mesh_skin_epi.num_points);
    
           
            %% Assign values %%
            %get the values of "data_name" in all three skin meshes
            values_epi=get_data_values_f(mesh_skin_epi,data_name,pointcell);
            values_right=get_data_values_f(mesh_skin_right,data_name,pointcell);
            
    
            %get the values to assign to "mesh_hex" based on
            %"closest_point_mesh"
            values_hex(closest_point_mesh(:,1)==1,:) = values_epi(closest_point_mesh(closest_point_mesh(:,1)==1,2),:);
            values_hex(closest_point_mesh(:,1)==2,:) = values_right(closest_point_mesh(closest_point_mesh(:,1)==2,2),:);
           
    
            %write the values in "mesh_hex"
            mesh_hex=write_data_(mesh_hex,mesh_hex.cells(:,1),values_hex,data_name,pointcell);
        
        case 'cell'
            %% Find Closest Cell to HexMesh in Skin Meshes%%
            %create 'skin_centroids' matrix with an ascending vector form 1
            %to the sum of the number of cells in all 3 skin meshes in 
            %the first column and with the coordinates of the centroids of 
            %all three skin meshes in the second, third and fourth columns
            skin_centroids(:,1) = [1:sum([mesh_skin_epi.num_cells,mesh_skin_right.num_cells])]';
            skin_centroids(:,2:4) = [mesh_skin_epi.centroids(:,2:4);mesh_skin_right.centroids(:,2:4)];
    
            %get the closest cell in the skin meshes to each centroid in 
            %"mesh_hex"
            [closest_cell,~,~] = find_closest_point(mesh_hex.centroids,skin_centroids);
    
            %create a matrix that contains in the first column an integer 
            %refering to the skin mesh that is the closest to the 
            %"mesh_hex" point (1 for epi, 2 for right) 
            %and in the second column its adjusted closest cell value
            closest_cell_mesh(closest_cell(:,1)<=mesh_skin_epi.num_cells,1) = 1;
            closest_cell_mesh(closest_cell_mesh(:,1)==1,2) = closest_cell(closest_cell_mesh(:,1)==1,1);
    
            closest_cell_mesh((closest_cell(:,1)>mesh_skin_epi.num_cells & closest_cell(:,1)<=(mesh_skin_epi.num_cells+mesh_skin_right.num_cells)),1) = 2;
            closest_cell_mesh(closest_cell_mesh(:,1)==2,2) = closest_cell(closest_cell_mesh(:,1)==2,1)-double(mesh_skin_epi.num_cells);
    
            closest_cell_mesh(closest_cell(:,1)>(mesh_skin_epi.num_cells+mesh_skin_right.num_cells),1) = 3;
            closest_cell_mesh(closest_cell_mesh(:,1)==3,2) = closest_cell(closest_cell_mesh(:,1)==3,1)-double(mesh_skin_epi.num_cells+mesh_skin_right.num_cells);
    
            %% Assign values %%
            %get the values of "data_name" in all three skin meshes
            values_epi=get_data_values_f(mesh_skin_epi,data_name,pointcell);
            values_right=get_data_values_f(mesh_skin_right,data_name,pointcell);
           
    
            %get the values to assign to "mesh_hex" based on
            %"closest_point_mesh"
            values_hex(closest_cell_mesh(:,1)==1,:) = values_epi(closest_cell_mesh(closest_cell_mesh(:,1)==1,2),:);
            values_hex(closest_cell_mesh(:,1)==2,:) = values_right(closest_cell_mesh(closest_cell_mesh(:,1)==2,2),:);
            
    
            %write the values in "mesh_hex"
            mesh_hex=write_data_f(mesh_hex,mesh_hex.cells(:,1),values_hex,data_name,pointcell);
            
        otherwise
            %raise an error
            errordlg("pointcell value must be 'point' or 'cell'",'transfer_skin_hex()');
            
            %end the function
            return;
    end
end


function mesh = add_data_f(mesh,type,name,format,width,pointcell)

    %depending on the value of "pointcell"
    switch pointcell
        
        case 'cell'
            %To add the new field to the end of the previous ones, "index" is set
            %as an inmediately superior value to the height of "cell_data_text"
            index=size(mesh.cell_data_text, 1)+1; %height(mesh.cell_data_text)+1;
    
            %add the data description values to "cell_data_text"
            mesh.cell_data_text(index,1)=type;
            mesh.cell_data_text(index,2)=name;
            mesh.cell_data_text(index,3)=format;
            mesh.cell_data_width(index,1)=width;
            
            %put the values reserved in "cells" for the data field to nan
            mesh.cells(:,(index-1)*10+1:(index)*10)=nan;
            
            %put the values occupied in "cells" by the data field to zero
            mesh.cells(:,(index-1)*10+1:(index-1)*10+width)=0;
            
            
        case 'point'
            %To add the new field to the end of the previous ones, "index" is set
            %as an inmediately superior value to the height of "point_data_text"
            index=size(mesh.point_data_text,1)+1; % height(mesh.point_data_text)+1;
    
            %add the data description values to "point_data_text"
            mesh.point_data_text(index,1)=type;
            mesh.point_data_text(index,2)=name;
            mesh.point_data_text(index,3)=format;
            mesh.point_data_width(index,1)=width;
            
            %put the values reserved in "points" for the data field to nan
            mesh.points(:,(index-1)*10+1:(index)*10)=nan;
            
            %put the values occupied in "points" by the data field to zero
            mesh.points(:,(index-1)*10+1:(index-1)*10+width)=0;
            
        case 'both'
            %To add the new field to the end of the previous ones, "index" is set
            %as an inmediately superior value to the height of "cell_data_text"
            index=size(mesh.cell_data_text,1)+1; %height(mesh.cell_data_text)+1;
    
            %add the data description values to "cell_data_text"
            mesh.cell_data_text(index,1)=type;
            mesh.cell_data_text(index,2)=name;
            mesh.cell_data_text(index,3)=format;
            mesh.cell_data_width(index,1)=width;
            %put the values reserved in "cells" for the data field to nan
            mesh.cells(:,(index-1)*10+1:(index)*10)=nan;
            %put the values occupied in "cells" by the data field to zero
            mesh.cells(:,(index-1)*10+1:(index-1)*10+width)=0;
            
            %To add the new field to the end of the previous ones, "index" is set
            %as an inmediately superior value to the height of "point_data_text"
            index=size(mesh.point_data_text,1)+1; %height(mesh.point_data_text)+1;
    
            %add the data description values to "point_data_text"
            mesh.point_data_text(index,1)=type;
            mesh.point_data_text(index,2)=name;
            mesh.point_data_text(index,3)=format;
            mesh.point_data_width(index,1)=width;
            %put the values reserved in "points" for the data field to nan
            mesh.points(:,(index-1)*10+1:(index)*10)=nan;
            %put the values occupied in "points" by the data field to zero
            mesh.points(:,(index-1)*10+1:(index-1)*10+width)=0;
            
        otherwise
            
            %raise an error
            errordlg('pointcell values must be cell, point or both','add_field()');
                
            %end the function
            return;
    end
end