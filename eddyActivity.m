function [eddyPath,eddyNodes,G] = eddyActivity(trakFilePath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    fid = fopen(trakFilePath);
    tline = fgetl(fid);
    Data=[];

    while ischar(tline)
        out = split(tline)';
        if(out{1}~="Frame")
            cellData=[cellfun(@str2num,out,'un',0)];
            out=cell2mat(cellData);
            Data=[Data;[out,zeros(1,15-length(out))]];
        else
            Data=[Data;zeros(1,15)];
        end
        tline = fgetl(fid);
    end
    fclose(fid);



    [y,w] = size(Data);
    G = digraph();
    time = 1;

    i = 1;
    while i <= y
        if (Data(i,1)) == 0
            i = i + 1;
            time = time + 1;
        end
        if Data(i,1) == -1

            eddyname = Data(i,2);
            newnode = append(string(time),'.',string(eddyname));
            G = addnode(G, num2str(newnode));

        elseif Data(i,2) == -1

            if (Data(i,3)) == 0
                if time == 2
                    neweddyname = Data(i,1);
                    newnode = append(string(time-1),'.',string(neweddyname));
                    G = addnode(G, num2str(newnode));
                end
            else
                oldeddyname = Data(i,1);
                oldnode = append(string(time-1),'.',string(oldeddyname));

                for j = 3:w
                    if (Data(i,j)) ~= 0
                        neweddyname = Data(i,j);
                        newnode = append(string(time),'.',string(neweddyname));
                        G = addnode(G, num2str(newnode));
                        G = addedge(G, num2str(oldnode),num2str(newnode));
                    end
                end
            end


        else
            j=find(Data(i,:) == -1);

            neweddyname = Data(i,j+1);
            newnode = append(string(time),'.',string(neweddyname));
            G = addnode(G,num2str(newnode));
            old = j-1;
            while old >= 1
                oldeddyname = Data(i,old);
                oldnode = append(string(time-1),'.',string(oldeddyname));
                G = addedge(G, num2str(oldnode), num2str(newnode));
                old = old - 1;
            end
        end



        i = i + 1;
    end
    figure,plot(G,"Layout","layered");

    eddyPath=mygetpaths(G);
    eddyNodes=table2array(G.Nodes);
%     for activityIndex=1:1:length(eddyPathTemp)
%         eddyPathLine=eddyPathTemp(activityIndex);
%     test2=arrayfun(@(y) eddyNodes(y), test);
    
end

