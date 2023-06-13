function [varargout] = getInputs(varargin)
    for i=1:nargin
        varargout{i} = getInput(varargin{i});
    end
end

function [y] = getInput(x)
    
    chars = ["[","]"];
    if (x=="GI") % ----------- General Inputs -----------------------------
        TGen = readtable("..\1. Inputs\1. General Inputs.csv","Delimiter",",","ReadVariableNames",false);
        GI = string(table2array(TGen));
        GI = GI(:,2);
        newStr = erase(GI,chars);
        y = newStr;
    
    elseif (x=="BI") % -------- Read Building Inputs ----------------------
        T1 = readtable("..\1. Inputs\2. Building Inputs.csv");
        range = T1{:,1};
        rangeStr = string(range);
        T2 = removevars(T1,'Range');
        T2 = erase(string(table2array(T2)),chars);

        all_inputs = zeros(length(rangeStr),5);

        for i=(1:length(rangeStr))
            str = rangeStr(i);
            newStr = erase(str,chars);
            x = split(newStr,'-')';
            y = double(x);
            for j = (y(1):y(length(y)))
                all_inputs(j,:)= str2double(T2(i,:));
            end
            BI = all_inputs;
        end
        y = BI;
 
    elseif (x=="AI") % ----------- Analysis Inputs ------------------------
        TAna = readtable("..\1. Inputs\3. Analysis Inputs.csv","Delimiter",",","ReadVariableNames",false);
        AI = TAna{:,2};
        AI = string(erase(AI,chars));
        y = AI;
        
    elseif (x=="EI") % ----------- Event Inputs ---------------------------
        TEve = readtable("..\1. Inputs\4. Event Inputs.csv","Delimiter",",","ReadVariableNames",false);
        EI = TEve{:,2};
        EI = string(erase(EI,chars));
        y = EI;
        
    elseif (x=="OC") % ----------- Output Control Inputs ------------------    
        OPC = readtable("..\1. Inputs\5. OC Inputs.csv","Delimiter",",","ReadVariableNames",false);
        OPC = string(erase(OPC{:,2},chars));
        y = str2double(OPC);
        
    elseif (x=="VI")
        VI = readtable("..\1. Inputs\6. Vector Inputs.csv","Delimiter",",","ReadVariableNames",false);
        if(isempty(VI) || size(VI,2)<2)
            y = [];   
        else
            for t=1:size(VI,1)
                id = strsplit(string(VI{t,1}),'_');
                VI_2 = string(erase(VI{t,2},chars));
                fieldName = "var_"+string(t);
                y.(fieldName)  = [id(1),id(2),strsplit(VI_2,";")];
            end
        end
    
%     elseif (x=="TI")
%         TTMD = readcell('..\1. Inputs\7. TMD Inputs.xlsx');
%         TI.nTMD = TTMD{1,3};
%         range=3:(2+TI.nTMD);
%         TI.idTMD = (TTMD(range,1));
%         TI.locTMD = cell2mat(TTMD(range,2));
%         TI.dirTMD = (TTMD(range,3));
%         TI.targFreqTMD = (TTMD(range,4));
%         TI.massRatioTMD = cell2mat(TTMD(range,5));
%         y = TI;

    elseif (x=="TI")
        TTMD = readcell('..\1. Inputs\7. TMD Inputs.csv');
        chars = ["[","]"];
        TTMD = string(erase(TTMD,chars));

        TI.nTMD = size(TTMD,1)-1;
        TI.idTMD = num2cell(double(TTMD(2:end,1)));
        TI.locTMD = double(TTMD(2:end,2));
        TI.dirTMD = {TTMD{2:end,3}}';
        TI.targFreqTMD = {TTMD{2:end,4}}';
        TI.massRatioTMD = double(TTMD(2:end,5));
        y = TI;
    end 
end
    
