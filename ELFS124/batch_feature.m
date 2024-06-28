function [fail, succ]=batch_feature(input_folder,output_filename,dt,Tmax)
[tmp, lghPath] = size(input_folder);
if input_folder(lghPath) ~= '\'
    input_folder = strcat(input_folder, '\');
end


drt = dir(input_folder);
[fnum, tmp] = size(drt);

FEA=zeros(fnum-2,48);
Fname=cell(fnum-2,1);
succ = 0;

for i = 1:fnum
    if drt(i).isdir == 1
        continue;
    end
    [tmp,lghName] = size(drt(i).name);
     if (isempty(findstr(lower(drt(i).name), '.off'))) || ((findstr(lower(drt(i).name), '.off') + 3) ~= lghName)
            continue;
     end
    inName = strcat(input_folder, drt(i).name);
    ori_file=inName;
     
    [FEA(i-2,:)] = EFF_fea( ori_file,dt,Tmax); 
    Fname{i-2} = drt(i).name;
    save(output_filename,'FEA', 'Fname','-v7.3');
    succ=succ+1
end
end
