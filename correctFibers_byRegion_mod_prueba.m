
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AUTHORSHIP -------------------------------------------------------- %%%
%                  Author:     Sandra Perez Herrero
%           Creation date:     02/10/2025
%       Last modification:     02/10/2025
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
% 3) a partir de la anterior > threshold (si existe surfnodes = 1 // si no existe filters > extractregionsurface > threshold (CellFacesIds) = 5 ... cuidado con las cositas pequeñas flotando) > filters > extract surface > filters > triangulate > connectivity > threshold > guardar por separado epi y endo como .vtk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% STEP2: CHANGE THE FIBERS AT THE REGION
% Read surface meshes
%epiPathVtk='C:\Users\Sandra\Desktop\laura_fibers\14\epi.vtk'; % malla de superficie correspondiente al epicardio de la zona que queremos modificar
epiPathVtk='.\epi.vtk'; % malla de superficie correspondiente al epicardio de la zona que queremos modificar
endoPathVtk='.\endo.vtk';% malla de superficie correspondiente al endocardio de la zona que queremos modificar

% Read voluetric meshes
fullPathInriaMesh='.\malla14_ajustada_HB_SC.vtk'; % malla completa en hexas con todas las etiquetas
regionPathVtk='.\paredla.vtk';% malla volumétrica con los hexahedros que contienen las fibras que quremos modificar. Esta malla proviene de la total y tiene todas las etiquetas que luego usamos para Elvira.

fullPathResultado='.\'; %path para guardar los resultados

%dirección del vector de paraview
fibre_data_vectors = [10.890044203249941,	18.03776294022534, -21.0107529328481,	-6.992776459971346, 10.834065316547747,	-23.140891348396593];
%Cargar la alineación de las fibras
fibre_data_align_epi='transversal' %puede tomar los valores transversal o longitudinal (direcciones respecto al vector introducido)
fibre_data_align_endo='longitudinal' %puede tomar los valores transversal o longitudinal (direcciones respecto al vector introducido)

correctFibers_byRegion_f(epiPathVtk,endoPathVtk,regionPathVtk,fullPathInriaMesh,fullPathResultado,fibre_data_vectors,fibre_data_align_epi,fibre_data_align_endo)


%%

function correctFibers_byRegion_f(epiPathVtk,endoPathVtk,regionPathVtk,fullPathInriaMesh,fullPathResultado,fibre_data_vectors,fibre_data_align_epi,fibre_data_align_endo)

% almacenar la información de las mallas de superficie en una estructura
% tipo mesh
mesh_epi = vtk2meshReader_f(epiPathVtk,'tri');
mesh_endo = vtk2meshReader_f(endoPathVtk,'tri');



%%  Calcular "Fibre" endocardios y epicardio

%asignacion de fibras
mesh_epi=assign_Fibre_f(mesh_epi,fibre_data_vectors,fibre_data_align_epi);
mesh_endo=assign_Fibre_f(mesh_endo,fibre_data_vectors,fibre_data_align_endo);


%Convertimos a vtk para asegurarnos de que se etiquetan bien
mesh2vtkWriter_f(mesh_epi,epiPathVtk,'epi_fibers.vtk');
mesh2vtkWriter_f(mesh_endo,endoPathVtk,'endo_fibers.vtk');

%% 3.3. Transferir Fibers desde las skin

region_mesh_hex=vtk2meshReader_f(regionPathVtk,'hex');

%create Fiber in cells

region_mesh_hex = transfer_skin_hex_f(region_mesh_hex,mesh_epi,mesh_endo,"Fibres",'cell');
mesh2vtkWriter_f(region_mesh_hex,fullPathResultado,'region_mesh_hex.vtk');



region_mesh_hex=vtkToolbox2Ci2b_vtkReader(strcat(fullPathResultado,'region_mesh_hex.vtk'),1);
region_mesh_hex_original=vtkToolbox2Ci2b_vtkReader(regionPathVtk,1);

valIdx_region_bad   = find(strcmp({region_mesh_hex.CellData.Name}, 'Fibres'));
valIdx_region_good   = find(strcmp({region_mesh_hex_original.CellData.Name}, 'Fibres'));
region_mesh_hex_original.CellData(valIdx_region_good).Data=region_mesh_hex.CellData(valIdx_region_bad).Data;
vtkToolbox2Ci2b_vtkWriter(region_mesh_hex_original,fullPathResultado,'region_mesh_hex_original.vtk');


%% 3.4 Transferir a la malla completa solo las nuevas fibras
mesh_hex=vtkToolbox2Ci2b_vtkReader(fullPathInriaMesh,1);

%%
% detectar posiciones de las variables por nombre
idIdx    = find(strcmp({region_mesh_hex_original.CellData.Name}, 'IDRef'));
valIdx   = find(strcmp({region_mesh_hex_original.CellData.Name}, 'Fibres'));
targetIdx= find(strcmp({mesh_hex.CellData.Name}, 'IDRef'));
outIdx   = find(strcmp({mesh_hex.CellData.Name}, 'Fibres'));

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
    mesh=struct2meshReader_f(struct,vtktype);
    
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
    normal_data=get_data_values_f(mesh,"Normals",'cell');

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
        mesh=write_data_f(mesh,target_cells,fibres,"Fibres",'cell');
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
    
    function mesh = struct2meshReader_f(struct,type)
    %% Check if struct is suitable %%    
    % Check if all cells are of the same type
    if (range(cell2mat(struct.Cells(:,1)))~=0) %if it is not
        
        %raise an error
        errordlg('Different types of cells detected in struct','struct2meshReader()')
        
        %end the function
        return
    end
    
    % Check if cells correspond to the "type" specified
    switch type
        case 'tri'
            % Check if cells are triangles
            if cell2mat(struct.Cells(1,1))~=3 %if they are not
                
                %raise an error
                errordlg('Type tri specified, but cells are not triangles','struct2meshReader()');
        
                %end the function
                return    
            end
            
        case 'hex'
            % Check if cells are hexahedra
            if cell2mat(struct.Cells(1,1))~=8 %if they are not
                
                %raise an error
                errordlg('Type hex specified, but cells are not hexahedra','struct2meshReader()');
        
                %end the function
                return    
            end
    end
    %% Dataset %%
    %Check if there is "Dataset" in "struct"
    if ~isempty(struct.Dataset) %if there is
        
        %save "dataset"
        mesh.dataset=struct.Dataset;
    else
        
        %raise an error
        errordlg('No "Dataset" specified in struct','struct2meshReader()');
        
        %end the program
        return
    end
    
    %% Number of Points %%
    %Check if there is "NumPoints" in "struct"
    if ~isempty(struct.NumPoints) %if there is
        
        %save "num_points" as "NumPoints"
        mesh.num_points=struct.NumPoints;
    else
        
        %raise an error
        errordlg('No "NumPoints" specified in struct','struct2meshReader()');
        
        %end the program
        return
    end
    
    %% Points %%
    %Check if there are "Points" in "struct"
    if ~isempty(struct.Points)
        
        %initialize points as a NaN matrix with 10 columns and num_points
        %rows
        mesh.points=nan(mesh.num_points,10);
        
        %the first column of points will be and ascending vector from 1 to
        %num_points and the second,third and fourth columns will contain the
        %coordinates of each point
        mesh.points(:,1:4)=[double([1:mesh.num_points]') struct.Points];
        
    else
        
        %raise an error
        errordlg('No "Points" specified in struct','struct2meshReader()');
        
        %end the program
        return
        
    end
    
    %% Cell Types: Vertices %%
    %Depending on the type specified
    switch type
        case 'tri'
            
            %set "cell_type_vert" to 3
            mesh.cell_type_vert=3;
            
        case 'hex'
            %set "cell_type_vert" to 8
            mesh.cell_type_vert=8;
    end
    
    %% Cell Types: Edges %%
    %Check if there is "CellTypes" in "struct"
    if any(string(fieldnames(struct))=='CellTypes') %if there is
        
        %save "cell_types_edges"
        mesh.cell_type_edges=struct.CellTypes;
    else
        mesh.cell_type_edges=[];
    end
    
    %% Number of Cells %%
    %Check if there is "NumCells" in "struct"
    if ~isempty(struct.NumCells) %if there is
        
        %save "num_cells"
        mesh.num_cells=struct.NumCells;
    else
        
        %raise an error
        errordlg('No "NumCells" specified in struct','struct2meshReader()');
        
        %end the program
        return
    end
    
    %% Cells %%
    %Check if there are "Cells" in "struct"
    if ~isempty(struct.Cells) %if there are
        
        %initialize cells as a NaN matrix with 10 columns and num_cells
        %rows
        mesh.cells=nan(mesh.num_cells,10);
        
        %the first column of cells will be an ascending vector from 1 to
        %num_cells and the rest will contain its vertices adjusted to +1 to
        %match the values of IDRef
        mesh.cells(:,1:mesh.cell_type_vert+1)=[[1:struct.NumCells]' cell2mat(struct.Cells(:,2)')'+1];
            
    else
        
        %raise an error
        errordlg('No "NumCells" specified in struct','struct2meshReader()');
        
        %end the program
        return
    end
    
    %% Edges %%
    %check if type is 'tri'
    if (isequal(type,'tri')) %if it is
        
        %compute the edges of the struct
        mesh.edges=extract_edges_f(mesh);
        
    else %if it is not
        
        %do not compute the edges
        mesh.edges=[];
    end
    
    %% Centroids %%
    %compute centroids of struct
    mesh.centroids=extract_centroids_f(mesh);
    
    %% Point Data %%
    %Check if there is "PointData" in "struct"
    if ~isempty(struct.PointData) %if there is
        
        %create a copy of PointData in point_info
        point_info=struct.PointData;
        
        %create a string array of the "Name" fields in "point_info" for better
        %data management
        point_info_names=string({point_info.Name});
        
        %Check if "IDRef" already exists in point_info_names
        if any(point_info_names=="IDRef")
                        
            %delete its entry from point_info
            point_info(point_info_names=="IDRef")=[];
            
        end
        
        %move all data in point_info to an inmediately superior index
        point_info(2:end+1)=point_info;
        
        %add IDRef data to the first entry of point_info
        point_info(1).Type='SCALARS';
        point_info(1).Name='IDRef';
        point_info(1).Format='double';
        point_info(1).Data=[1:mesh.num_points]';
        
        %initialize "point_data_text" as a nan matrix with the same number
        %of rows as "point_info" datasets are and 3 columns and convert it
        %to a string matrix
        mesh.point_data_text=nan(size(point_info,2),3);
        mesh.point_data_text=string(mesh.point_data_text);
        
        %initialize "point_data_width" as a nan matrix with the same number
        %of rows as datasets are "point_info" and 1 column
        mesh.point_data_width=nan(size(point_info,2),1);
        
        %add to "points" 10 times the number of "point_info"-1 fields as columns
        %of nan values
        mesh.points=[mesh.points, nan(mesh.num_points,10*(length(point_info)-1))];
        
        %for each field in "point_info"
        for i=1:size(point_info,2)
            
            %save all data related in "point_data_text"
            mesh.point_data_text(i,:)=[string(point_info(i).Type) string(point_info(i).Name) string(point_info(i).Format)];
            
            %save the width of the data in "point_data_text"
            mesh.point_data_width(i,:)=size(point_info(i).Data,2);
            
            %save its data in "points" with the first column on the (i-1)*10+1
            %column
            mesh.points(:,(i-1)*10+1:(i-1)*10+mesh.point_data_width(i))=point_info(i).Data;
        end
        
    else %if there is no "PointData"
        
        %initialize "point_data_text" as a nan matrix with 1 row and 3
        %columns and convert it to a string matrix
        mesh.point_data_text=nan(1,3);
        mesh.point_data_text=string(mesh.point_data_text);
        
        %initialize "point_data_width" as a nan matrix 1 row and 1 column
        mesh.point_data_width=nan(1,1);
    
        %save all data related defining "IDRef"
        mesh.point_data_text(1,1)='SCALARS';
        mesh.point_data_text(1,2)='IDRef';
        mesh.point_data_text(1,3)='double';
            
        %save the width of the data in "cell_data_text"
        mesh.point_data_width(1)=1;
        
    end
    
    %% Cell Data %%
    %Check if there is "CellData" in "struct"
    if ~isempty(struct.CellData) %if there is
        
        %create a copy of CellData in cell_info
        cell_info=struct.CellData;
        
        %create a string array of the Name fields in cell_info for better
        %data management
        cell_info_names=string({cell_info.Name});
        
        %Check if "IDRef" already exists in cell_info_names
        if any(cell_info_names=="IDRef")
                        
            %delete its entry from cell_info
            cell_info(cell_info_names=="IDRef")=[];
            
        end
        
        %move all data in cell_info to an inmediately superior index
        cell_info(2:end+1)=cell_info;
        
        %add IDRef data to the first entry of cell_info
        cell_info(1).Type='SCALARS';
        cell_info(1).Name='IDRef';
        cell_info(1).Format='double';
        cell_info(1).Data=[1:mesh.num_cells]';
        
        %initialize "cell_data_text" as a nan matrix with the same number of
        %rows as "cell_info" datasets are and 3 columns and convert it to a
        %string matrix
        mesh.cell_data_text=nan(length(cell_info),3);
        mesh.cell_data_text=string(mesh.cell_data_text);
        
        %initialize "cell_data_width" as a nan matrix with the same number of
        %rows as "CellData" datasets are in the struct and 1 column
        mesh.cell_data_width=nan(size(struct.PointData,2),1);
        
        %add to "cells" 10 times the number of "cell_info"-1 fields as columns
        %of nan values
        mesh.cells=[mesh.cells, nan(mesh.num_cells,10*(size(cell_info,2)-1))];
        
        %for each field in "cell_info"
        for i=1:length(cell_info)
            
            %save all data related in "cell_data_text"
            mesh.cell_data_text(i,:)=[string(cell_info(i).Type) string(cell_info(i).Name) string(cell_info(i).Format)];
            
            %save the width of the data in "cell_data_text"
            mesh.cell_data_width(i,:)=size(cell_info(i).Data,2);
            
            %save its data in "cells" with the first column on the (i-1)*10+1
            %column
            mesh.cells(:,(i-1)*10+1:(i-1)*10+mesh.cell_data_width(i))=cell_info(i).Data;
        end
        
    else %if there is no "CellData"
        
        %initialize "cell_data_text" as a nan matrix with 1 row and 3
        %columns and convert it to a string matrix
        mesh.cell_data_text=nan(1,3);
        mesh.cell_data_text=string(mesh.cell_data_text);
        
        %initialize "cell_data_width" as a nan matrix 1 row and 1 column
        mesh.cell_data_width=nan(1,1);
        
        %save all data related defining "IDRef"
        mesh.cell_data_text(1,1)='SCALARS';
        mesh.cell_data_text(1,2)='IDRef';
        mesh.cell_data_text(1,3)='double';
            
        %save the width of the data in "cell_data_text"
        mesh.cell_data_width(1)=1;
        
    end
    
    
    %% Field Data %%
    %Check if there is "Field" data in "struct"
    if ~isempty(struct.Field)
        
        %for each field in Field data
        for i=1:size(struct.Field,2)
            
            %check if it is corresponding to cell data or to point data
            %depending on the height of the "Data" field
            switch size(struct.Field(i).Data,1) % height(struct.Field(i).Data)
                case mesh.num_cells
                    data_field='cell';
                    target_id=[1:mesh.num_cells];
                    
                case mesh.num_points
                    data_field='point';
                    target_id=[1:mesh.num_points];
                    
                otherwise
                    data_field=[];
            end
            
            %check if it is a NORMAL or a SCALAR depending on the width of
            %the "Data" field
            switch size(struct.Field(i).Data,2)
                case 1
                    data_type='SCALARS';
                    data_width=1;
                    
                case 2
                    data_type='NORMALS';
                    data_width=3;
                    
                otherwise
                    data_type=[];
            end
            
            %if data_type and data_field can be detected
            if (~isempty(data_type) & ~isempty(data_field))
                
                %add data description
                mesh=add_data(mesh,data_type,struct.Field(i).Name,struct.Field(i).Format,data_width,data_field);
                
                %write data values
                mesh=write_data(mesh,target_id,struct.Field(i).Data,struct.Field(i).Name,data_field);
                
            end
            
        end
    end
    
end
 
function edges_unique = extract_edges_f(mesh)
    
    %check if "mesh" is a triangular mesh
    if mesh.cell_type_vert~=3 %if it is not
        
        %raise an error
        errordlg('mesh is not a triangular mesh','extract_edges()')
        
        %end the function without computing edges
        edges_unique=[];
        return;
    end
    
    %initialize "edges" as a nan matrix with 3 times as much rows as cells 
    %are in the mesh and 4 columns
    edges=nan(3*size(mesh.cells,1),4);

    %for each edge in a cell (in triangular cells)
    for i=1:size(mesh.cells,1)
        
        %save in "edges" the index and the current combination of edges,
        %leaving the fourth column empty
        edges((i-1)*3+1:(i-1)*3+3,1:3)=[repmat(i,3,1) nchoosek(mesh.cells(i,2:4),2)];
    end

    %save the indices of the rows whose end-point value in the second
    %columns is greater than its end-point value in the third column
    edges_swap=edges(:,2)>edges(:,3);
    
    %change the order of the vertices to get the smaller value on
    %the first column and the greater in the second one
    aux=edges(edges_swap,2);
    edges(edges_swap,2)=edges(edges_swap,3);
    edges(edges_swap,3)=aux;
    
    %delete duplicated edges and save its repetition indices
    [edges_unique,~,edges_unique_index]=unique(edges(:,2:3),'row');
    
    %save unique edges on columns 3 and 4
    edges_unique(:,3:4)=edges_unique(:,1:2);
    
    %add the column vector indicating edge repetition to the fourth column
    %of "edges"
    edges(:,4)=edges_unique_index;

    %divide "edges" based on the value of "edges(:,4)" and apply the
    %function "extract_1" to each group
    neighbours=splitapply(@extract_1,edges,edges(:,4));
    
    %transform "neighbours" from cell array to numerical matrix
    neighbours=cell2mat(neighbours);
    
    %add the identifiers of neighbours to the first and second columns of
    %edges_unique
    edges_unique(:,1:2)=neighbours;
    
    
end

%this is an auxiliary function that returns a cell containing the fisrt two
%values of the first column of a matrix "v" when its height is greater than
%1 and a cell containing its first column value and nan if its height is one
function extracted=extract_1(v)
    if size(v,1)>1
        extracted={v(1:2,1)'};
    else
        extracted={[v(1) nan]};
    end
end

function centroids = extract_centroids_f(mesh)

    %initialize "centroids" as a nan matrix with as much rows as cells are
    %in the mesh and 4 columns
    centroids=nan(size(mesh.cells,1),4);
    
    %for each cell
    for i=1:size(mesh.cells,1)
        
        %extract the vertices of the cell, which are allocated between the
        %second and the tenth column of "mesh.cells"
        vertices=mesh.cells(i,2:10);
        
        %delete all NaN elements
        vertices(isnan(vertices))=[];
        
        %extract from "mesh.points" the coordinates of each vertex
        vertices_coord=mesh.points(vertices',2:4);
        
        %"centroid" is a vector containing the average of each coordinate
        centroid=sum(vertices_coord,1)./length(vertices);
        
        %add the IDRef value of the cell and its "centroid" to "centroids"
        centroids(i,:)=[mesh.cells(i,1),centroid];
    end
end


