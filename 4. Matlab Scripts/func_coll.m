function func_coll(varargin) % collection of functions
  if nargin==1
     str = varargin{1};
     fun=str2func(str); fun();
  elseif nargin==2
     str = varargin{1};
     fun=str2func(str); fun(varargin{2}); 
  elseif nargin==3
     str = varargin{1};
     sampleNo = varargin{2};
     rv_vector = varargin{3};
     fun=str2func(str); fun(sampleNo,rv_vector);   
  end
end


function generate_Excitation()
   G = getInputs("GI");
   if(G(3)=="W")          % check if the analysis is for wind or seismic 
      run("..\4. Matlab Scripts\generate_wind_load"); 
      disp("Wind load - Generated!");
   else 
      disp("Seismic Load")
   end
end

function analyze(dir)
    A = getInputs("AI");
      if (A(5)=="0")                  % without TMD
        MDOF_Shear_IMK_Wind(dir);
      else                            % with TMD
        MDOF_Shear_IMK_Wind_TMD(dir);
      end
    
end

function run_()
    
    clearOutputFile();     % clear the output folder   
    [G,V] = getInputs("GI","VI");
    
    % If there are multiple valued parameters
    if(~isempty(V))
       
       for k = 1:numel(fieldnames(V))
          fieldName = "var_"+string(k);
          Vi = V.(fieldName);      
          n = length(Vi);
          inp_name = Vi(1);
          inp_row  = double(Vi(2));
          
          if(inp_name=="BI")
            dir = "..\1. Inputs\2. Building Inputs.csv";
          elseif(inp_name=="AI")
            dir = "..\1. Inputs\3. Analysis Inputs.csv";
          elseif(inp_name=="EI")
            dir = "..\1. Inputs\4. Event Inputs.csv";
          end
          
          for i = 3:n 
               id = strcat(Vi(1),Vi(2),"_",string(Vi(i)));
               id_dir = "..\2. Outputs\"+id;
               mkdir(id_dir);
               
           for ii = 1:double(G(4))     

               if(double(G(4))>1)    
                   id2 = "sim_"+string(ii);
                   id_dir2 = strcat(id_dir,"\",id2);
                   mkdir(id_dir2);
                    
                   p = strcat(id,"\",id2);
                   
                   id_path = strcat("\",p);
                   updateVal(dir,inp_row,Vi(i)); 
                   generate_Excitation();
                   analyze(id_path);                   
               else
                   id_path = strcat("\",id);
                   updateVal(dir,inp_row,Vi(i)); 
                   
                   generate_Excitation();
                   analyze(id_path);
               end
               
              
           end
          end
       end
       
    else
       if(double(G(4))>1) 
        w = waitbar(0,"starting simulation...");
         for i = 1:double(G(4))     
           id = "sim_"+string(i);
           id_dir = strcat("..\2. Outputs\",id);
           mkdir(id_dir);

           id_path = strcat("\",id);
           generate_Excitation();
           analyze(id_path);    
           
           w = waitbar(i/double(G(4)),w,"loading...");
         end
       else
           generate_Excitation();
           analyze("");
       end  
    end
end

function clearOutputFile()
    try 
      rmdir ("..\2. Outputs\*",'s');
    catch 
    end
    try
      delete("..\2. Outputs\*.csv"); 
    catch
    end
end
function updateVal(filePath,row,val)
     filePath = "..\1. Inputs\7. TMD Inputs.xlsx";
     table = readtable(filePath,"Delimiter",",","ReadVariableNames",false);
     str = strcat("[",string(val),"]");
     table{row,2} = {char(str)};
     writetable(table,filePath,'WriteVariableNames',0);
end
function collect_res(id_path)
    %% Collect Results --------------------------------------------------------
    data = readmatrix("..\6. OpenSees\bin\tempfiles\info.csv");
    nS = data(1);
    N = data(2)+1;
    dt = data(3);

    accel_TH = zeros(nS,N);
    disp_TH = zeros(nS,N);
    force_TH = zeros(nS,N);
    vel_TH = zeros(nS,N);
    
    max_absDisp = zeros(nS,1);
    rms_Disp = zeros(nS,1);
    
    max_absAccel = zeros(nS,1);
    rms_Accel = zeros(nS,1);
    
    max_absVel = zeros(nS,1);
    rms_Vel = zeros(nS,1);

%     max_absSF = zeros(nS,1);
%     rms_SF = zeros(nS,1);
    
    per=zeros(nS,1);
    timeVec=0:dt:(N-1)*dt;
    x = getInputs("OC");clc;
    
    writematrix(timeVec,"..\2. Outputs"+id_path+"\O1 TimeVector_TH.csv") % write the time vector as a csv file
    
    for i=1:nS
 
       if x(1)==1
            disp_TH(i,:)=readmatrix("..\6. OpenSees\bin\tempfiles\nodeDisp_"+i+".csv");
            if(i==nS)
              writematrix(disp_TH,"..\2. Outputs"+id_path+"\O2 RelativeDisplacement_TH.csv");
            end
       end

       SF=readmatrix("..\6. OpenSees\bin\tempfiles\shearForce_"+i+".csv");SF=SF(:,8);
       force_TH(i,:) = SF;
       
       if x(3)==1
           vel_TH(i,:)=readmatrix("..\6. OpenSees\bin\tempfiles\nodeVel_"+i+".csv");
            if(i==nS)
              writematrix(vel_TH,"..\2. Outputs"+id_path+"\O4 Velocity_TH.csv");       
            end
       end
       if x(4)==1
           accel_TH(i,:)=readmatrix("..\6. OpenSees\bin\tempfiles\nodeAccel_"+i+".csv");
           if(i==nS)
               writematrix(accel_TH,"..\2. Outputs"+id_path+"\O5 Acceleration_TH.csv");
           end
       end
       
       if x(10)==1 % max abs disp
          MD = readmatrix("..\6. OpenSees\bin\tempfiles\max_abs_disp_"+i+".csv"); MD = MD(3);
          max_absDisp(i,1) = MD;
          if(i==nS)
             writematrix(max_absDisp,"..\2. Outputs"+id_path+"\O11 MaxAbsDisplacement.csv");
          end
       end
       
       if x(11)==1 % rms disp
          RD = readmatrix("..\6. OpenSees\bin\tempfiles\rms_disp_"+i+".csv");
          rms_Disp(i,1)=RD;
          if(i==nS)    
            writematrix(rms_Disp,"..\2. Outputs"+id_path+"\O12 rmsDisplacement.csv");
          end
       end
       
       if x(12)==1 % max abs accel
          MA = readmatrix("..\6. OpenSees\bin\tempfiles\max_abs_acc_"+i+".csv");  MA = MA(3);
          max_absAccel(i,1)= MA;
          if(i==nS)
            writematrix(max_absAccel,"..\2. Outputs"+id_path+"\O13 MaxAbsAcceleration.csv");
          end
       end
       
       if x(13)==1 % rms accel
          RA = readmatrix("..\6. OpenSees\bin\tempfiles\rms_acc_"+i+".csv");    
          rms_Accel(i,1)= RA;
          if(i==nS)    
            writematrix(rms_Accel,"..\2. Outputs"+id_path+"\O14 rmsAcceleration.csv");
          end
       end
       
       if x(14)==1 % max abs vel
          MV = readmatrix("..\6. OpenSees\bin\tempfiles\max_abs_vel_"+i+".csv"); MV = MV(3);
          max_absVel(i,1)=MV;
          if(i==nS)    
            writematrix(max_absVel,"..\2. Outputs"+id_path+"\O15 MaxAbsVelocity.csv");
          end
       end
       
       if x(15)==1 % rms vel
          RV = readmatrix("..\6. OpenSees\bin\tempfiles\rms_vel_"+i+".csv");
          rms_Vel(i,1)=RV;
          if(i==nS)    
            writematrix(rms_Vel,"..\2. Outputs"+id_path+"\O16 rmsVelocity.csv");
          end
       end

    end

    if x(9)==1
      for i=[4:2:3+2*nS;1:nS]
       per(i(2))=data(i(1));
      end
      writematrix(per,"..\2. Outputs"+id_path+"\O10 FundamentalPeriod.csv");
    end
    
    if x(2)==1
       writematrix(force_TH,"..\2. Outputs"+id_path+"\O3 ShearForce_TH.csv");
    end
    if x(16)==1
       writematrix(max(abs(force_TH'))',"..\2. Outputs"+id_path+"\O17 MaxAbsShearForce.csv");
    end
    if x(17)==1
        writematrix(rms(force_TH')',"..\2. Outputs"+id_path+"\O18 rmsAbsShearForce.csv");
    end
        
end