classdef Imp
    %%
    properties
        time
        fs                %Freq of sampling
        NSample=5         %Number of sample for mean
        row               %Number of rows in Ze and Zem
        NFFT              %signal length for caculate FFT
        ImpL              %Impedance signal length
        ind               %Indices according to specific Frequency range
        ind2              %Indices according to specific test time range
        Fnew              %Freq according to specific Frequency range
        start_ind1
        start_ind2
        end_ind1
        end_ind2
        start_ind3
        end_ind3
        start_freq=20000  %Freq of start result impedance signal
        end_freq=100000   %Freq of end result impedance signal
        start_test=0.003  %Freq of start result impedance signal
        end_test=0.017    %Freq of end result impedance signal
        KaiserAlpha=5
        TukeyAlpha=0.5
        NFolder=1
        ResultFolder
        RMSD_FirstArg         %firt sample in each data folder
        SampleFile            %sample path for each interation
        Ze
        RMSD_StartFreq=40000  %start of freq RMSD range for anomaly detection
        RMSD_EndFreq=50000    %end of freq RMSD range for anomaly detection
        RMSD_inx              %index that RMSD caculated
        ZeRMSD                %impedance for RMSD
        FnewRMSD
        metric=0              %RMSD value
        RMSDthresh=0.075       %RMSD threshold that determine signal is normal or not.
        AbnormalSignal=false  %Actives when detects an abnormal signal.
        cropping=false        %For using in window method
    end
    %%
    methods
        %%Changing file name
        function changename(~,datafile)
            global path
            
            % manage paths
            % for id=[0 1 17]
            %   a='C:\Users\mehrab\Desktop\';
            %     b='%d';
            %     b=sprintf(b,id);
            %   datafile='C:\Users\mehrab\Desktop\aaa';
            
            path.results =datafile ;
            currentPath = fullfile(path.results);
            currentruns = dir(currentPath);
            
            % ******************** count data files; renaming*****************
            oldnames = {currentruns.name};
            validnames.idx = find(contains(oldnames, '.osc.txt'));
            oldnamesValid = oldnames(validnames.idx);
            noRuns = numel(validnames.idx);
            for k = 1:numel(oldnamesValid)
                if ~contains(oldnamesValid{k},  'run')
                    movefile(fullfile(currentPath, oldnamesValid{k}),...
                        fullfile(currentPath, ['run', num2str(k) '.txt']));
                end
            end
        end
        %% Loading data without initialization(more computation cost than with initialization
        function [self,A,B]=LoadData(self,datafile)
            % Loading the data
            aa=load(datafile);
            A=aa(2:end,1);
            B=aa(2:end,2);
            self.fs = aa(1,1); %sample rate
            Tinterval=1/self.fs;
            self.time = Tinterval: Tinterval:length(B)*Tinterval;
        end
        %% Loading data with initialization(for decresing computation cost)
        function [A,B]=LoadDataIni(self)
            % Loading the data
            aa=load(self.SampleFile);
            A=aa(2:end,1);
            B=aa(2:end,2);
        end
        %%
        function [self,tnew,Anew,Bnew]=FindSignal(self,input,output)
            % Spectrogram using STFT
            fstart = 0;    % intended start freq.
            fend = 200000;    % intended final freq.
            res = 64;      % resulution
            Tinterval=1/self.fs;
            window = 1024*4 ; % windows length in STFt
            overlap = 512*4; % overlap between consecutive windows
            FTres = fstart:res:fend;
            [~, ~, T, ps] = spectrogram(input, window, overlap,FTres,self.fs, 'yaxis');
            
            % figure(2)
            % spectrogram(A, window, overlap,FFT,fs, 'yaxis');
            
            % Extracting the start and end indicess of one sweep
            tresh1 = 1e-3;
            ps(ps<tresh1) = 0;
            % % ps=1000*ps;
            [row,col,~] = find(ps);
            
            tresh2 = 150;
            abruptChangeCol = [];
            for i =2: length(row)
                if abs(row(i) - row(i-1)) > tresh2
                    abruptChangeCol = [abruptChangeCol, col(i)];
                end
            end
            self.AbnormalSignal=true; %check if signal is incomplete
            C=0.6; %desiered signal length to all data length ratio
            abruptChangeTime = T(abruptChangeCol);
            for i=2:size(abruptChangeTime,2)
                if abruptChangeTime(i)-abruptChangeTime(i-1)>C*size(input,1)*Tinterval
                    timeRange(1)=abruptChangeTime(i-1);
                    timeRange(2)=T(abruptChangeCol(i)-1);
                    self.AbnormalSignal=false;
                    break
                end
            end
            switch true
                case ~self.AbnormalSignal
                    startIndex = find (abs(self.time-timeRange(1))<Tinterval);
                    endIndex = find (abs(self.time-timeRange(2))<Tinterval);
                    tnew=self.time(startIndex:endIndex)';
                    Anew=input(startIndex:endIndex);
                    Bnew=output(startIndex:endIndex);
                    %Find number of FFT for Impedance FFT
                    L1=length(Anew);
                    L2=length(Bnew);
                    self.NFFT=max([2^21 2^nextpow2(L1) 2^nextpow2(L2)]);
                otherwise  %If signal is abnormal
                    fileID = fopen('ErrorReport.txt','a'); % Open the file for appending
                    fprintf(fileID,'%s\n',['Uncomplete test duration  in ',self.SampleFile,' file']);
                    fclose(fileID); % Close the file
                    tnew=[];
                    Anew=[];
                    Bnew=[];
                    
            end
        end
        %%
        function [adm,V1,V2]=Admittance(self,Anew,Bnew)
            %---calculate impedance--------------------
            V1=(fft(Anew,self.NFFT));
            V2=(fft(Bnew,self.NFFT));
            %Ze=((V1-V2)./V2);
            adm=(V2./(V1-V2));
            L=size(adm,1);
            adm = adm(1:floor(L/2));
            % V1=V1(1:floor(L/2));
            % V2=V2(1:floor(L/2));
            %  FZe=smooth(FZe,10);
        end
        %%
        function [Zem,V1m,V2m]=MeanImp(self,V1,V2)
            switch nargin
                case 4
                    Zem=mean(self.Ze,2);
                    V1m=mean(V1,2);
                    V2m=mean(V2,2);
                case 3
                    Zem=mean(self.Ze,2);
                    V1m=mean(V1,2);
                    
                case 2
                    Zem=mean(self.Ze,2);
                    %                     [ImpNew,FNew]=CutImp(self,Zem);
                otherwise
                    error('input arguments are not correct')
            end
        end
        %% Finding variables that are constant in every Sample.
        function [self,In,Out]=initial(self,FilePath)
            %             datafile=[FilePath,int2str(self.start_ind1),'\run',int2str(self.start_ind2),'.txt'];
            datafile=FilePath;
            [self,input,output]=LoadData(self,datafile);     %Find fs,time
            [self,tnew,~,~]=FindSignal(self,input,output);   %find NFFT and test duration
            Tinterval=1/self.fs;
            %% Finding specific range of test that doesnt have noise.
            testL=length(tnew);
            TestTime=0:Tinterval:(testL-1)*Tinterval;
            self.ind2=find(TestTime>self.start_test&TestTime<self.end_test);
            %             Anew=zeros(length(self.ind2),1);
            %             Bnew=Anew;
            %% finding indices of specific frequency range of Impedance
            self.ImpL=self.NFFT/2;                           %Number of impedance length after fft
            f = self.fs*(0:floor(self.ImpL)-1)/(2*self.ImpL);
            self.ind=find(f>self.start_freq&f<self.end_freq);
            self.RMSD_inx=find(f>self.RMSD_StartFreq&f<self.RMSD_EndFreq);
            self.Fnew=f(self.ind);
            self.FnewRMSD=f(self.RMSD_inx);
            %%
            self.row=length(self.ind);                            %Num of Re and Rem rows
            %             col2=floor((self.start_ind1-self.end_ind1)*(self.start_ind2-self.end_ind2)/5);
            self.Ze=zeros(self.row,self.NSample);
            In=self.Ze;
            Out=self.Ze;
            %             Zem=zeros(self.row,col2);
        end
        %%
        function plotZe(self)
            L=size(self.Ze,1);
            %             f = self.fs*(0:floor(L)-1)/(2*L);
            f=self.Fnew;
            %figure
            plot(f,self.Ze)
            % xlabel("Frequency(Hz)")
            % ylabel("Conductance")
            xlim([20000 100000])
        end
        %% Cut signal for specific frequency range
        function [ImpNew,FNew]=CutImp(self,Zem)
            L=size(Zem,1);
            f = self.fs*(0:floor(L)-1)/(2*L);
            self.ind=find(f>self.start_freq&f<self.end_freq);
            ImpNew=Zem(self.ind);
            FNew=f(self.ind);
            %save("pristine1.mat","admitance3","FNew","-append")
        end
        %%
        function PlotData(self,A,B,Anew,Bnew,tnew)
            %Ploting input and output signal in two seprate figure
            %Each figure plot signal and modified signal together
            %% Plot new signal
            figure
            plot (self.time , A)
            xticks=0:0.1:1;
            xticks = 1000*xticks;
            xticklabels(xticks);
            hold on
            plot (tnew, Anew)
            figure
            plot (self.time , B)
            %clear xticks
            xticks = 1000*xticks;
            xticklabels(xticks);
            hold on
            plot (tnew, Bnew)
            %% specrogram of New signal
            fstart = 0;           % intended start freq.
            fend = 200000;        % intended final freq.
            res = 64;             % resulution
            window = 1024*8 ; % windows length in STFt
            overlap = 1024*6; % overlap between consecutive windows
            FTres = fstart:res:fend;
            figure
            spectrogram(Anew, window, overlap,FTres,self.fs, 'yaxis')
            figure
            spectrogram(Bnew, window, overlap,FTres,self.fs, 'yaxis')
            
        end
        %%
        function [WinIn,WinOut]=windowingData(self,input,output,win)
            %Used for appling window before calculating the FFT of signals
            %For plottig data add time as third function arg
            switch win
                case 'kaiser'
                    win=kaiser(length(input),self.KaiserAlpha);
                case 'tukey'
                    win=tukeywin(length(input),(self.TukeyAlpha-0.2*self.cropping));
            end
            WinIn=win.*input;
            WinOut=win.*output;
        end
        %%
        function [self,FileNum,inx1,inx2]=analisys(self,input,output,inx1,inx2,FileNum,opt)
            
            arguments
                %%
                self
                input
                output
                inx1
                inx2
                FileNum
                opt.window='rectangular'
                opt.cropping=false  %To crop signal and start and end of the signal
                opt.visualization=false %For plotting In/Out signal in time domain and spectrogram
                opt.A  %Insert it if Visualization is true...
                %                 Input signal without finding start and end of the signal
                opt.B  %Insert it if Visualization is true...
                %                 Output signal without finding start and end of the signal
                opt.tnew
                opt.InOut=false %to saving Input output signal besides admittance
                opt.In  %For ussing in code
                opt.Out %For ussing in code
                opt.RMSD=false %to save first sample of each damage case as reference signal and ...
                %compare it with other signals in the same case
                opt.Sfile
                
            end
            
            
            %%      cropping signal
            %         For cropping signal in time domain use following lins...
            %            (ind2 is indecies of cropped test time):
            if(opt.cropping)
                input=input(self.ind2);
                output=output(self.ind2);
                opt.tnew=opt.tnew(self.ind2);
                self.cropping=true;
            end
            %                                  tnew=tnew(self.ind2); %for using plot data function active this line
            %%
            if(isequal(opt.window,'kaiser') || isequal(opt.window,'tukey'))
                %     Appling window
                %For appling window before calculating the FFT of them
                [input,output]=self.windowingData(input,output,opt.window);
                %         [input,output,tnew]=windowingData(self,input,output,'kaiser',tnew); %for using PlotdData function active this
                %                                                                    line instead of pervious line
            end
            %%     Plotting test data
            if(opt.visualization)
                self.PlotData(opt.A,opt.B,input,output,opt.tnew)
                
            end
            %% Compute admittance and assign freq range to signals
            [ZeTot,In,Out]=self.Admittance(input,output);
            [self.Ze(:,inx1)]=ZeTot(self.ind);
            [self.ZeRMSD]=ZeTot(self.RMSD_inx);
            if(opt.InOut)
                [opt.In(:,inx1)]=In(self.ind);
                [opt.Out(:,inx1)]=Out(self.ind);
            end
            if(opt.RMSD)
                if(inx1==1)
                    self.RMSD_FirstArg=self.ZeRMSD;
                end
                self.metric=self.RMSD();
            end
            %% check for unusual signal
            %To compare reference signal  with other signals in the same case
            switch true
                case self.metric<self.RMSDthresh %If signal is normal the rest of the code runs.
                    %%     Compute impdance average every 'NSalmple' sample
                    if rem(inx1,self.NSample)==0
                        [Zem]=self.MeanImp(self.Ze);
                        %                 writematrix(Zem,[self.ResultFolder{FileNum},'\admittance',int2str(inx2),'.csv'])
                        %                         Zem=round(Zem, 6);
                        Zem=single(Zem);
                        save([self.ResultFolder{FileNum},'\Admittance',int2str(inx2),'.mat'],'Zem');
                        if(opt.InOut)
                            [Inm,Outm]=self.MeanImp(opt.In,opt.Out);
                            Inm=single(Inm);
                            Outm=single(Outm);
                            save([self.ResultFolder{FileNum},'\Input',int2str(inx2),'.mat'],'Inm');
                            save([self.ResultFolder{FileNum},'\Output',int2str(inx2),'.mat'],'Outm');
                            
                        end
                        
                        FileNum=FileNum+1;
                        %             ZeTotal(:,:,s)=Ze;  % To get all Ze signal in one variables
                        inx2=inx2+1;
                        %                 self.Ze=zeros(self.row,self.NSample);
                        inx1=0;
                    end
                    inx1=inx1+1;
                    
                otherwise  %If signal is abnormal
                    fileID = fopen('ErrorReport.txt','a'); % Open the file for appending
                    fprintf(fileID,'%s\n',['Noisy signal in ',opt.Sfile,' file']);
                    fclose(fileID); % Close the file
                    self.AbnormalSignal=true;
            end
            
            
        end
        %%
        %         function datafile=createfolderType2(~,parentfolder,varargin)
        %             %Create folder for saving data in analisys function. Number of
        %             %folder is depend on number of input arg.
        %
        %             for i=1:(nargin-2)
        %                 datafile=[parentfolder varargin{i}];
        %                 mkdir(datafile)
        %             end
        %         end
        %%
        %         function datafile=createfolder(self,parentfolder,FileVector)
        %             %Creating folder with names contain in FileVector entries. folders are created in all of
        %             %parentfolder entries.
        %             %result parentfolder is a cell array with length(FileVector)*length(parentfolder)
        %             %
        %
        %             inx=1;
        %             datafile=cell(length(FileVector)*length(parentfolder),1);
        %             for j=1:length(parentfolder) %loop over pervious folders to creating same folders in...
        %                 %                 each folders
        %
        %                 for i=1:length(FileVector)
        %                     datafile{inx}=[parentfolder{j} '\' FileVector{i}];
        %                     mkdir(datafile{inx})
        %                     inx=inx+1;
        %                     %             end
        %                 end
        %             end
        %         end
        %%
        function self=CreateFolder(self,AnalisysMode,parentfolder,datafile1,FilePath)
            self.ResultFolder=cell(length(AnalisysMode),1);
            for mode=1:length(AnalisysMode)
                self.ResultFolder{mode}=[parentfolder{1},'\',AnalisysMode{mode},datafile1((length(FilePath)+1):end)];
                mkdir(self.ResultFolder{mode});
            end
        end
        %%
        function aa=RMSD(self)
            aa=(sum((real(self.RMSD_FirstArg)-real(self.ZeRMSD)).^2)./sum(real(self.ZeRMSD).^2))^0.5;
        end
        %%
        function save(self)
            fileID = fopen('ErrorReport.txt','a'); % Open the file for appending
            fprintf(fileID,'%s\n',['Uncomplete test duration  in ',self.SampleFile,' file']);
            fclose(fileID); % Close the file
        end
        
        %% Finding RMSD of 2 signal and plot bar graph
        function Y=plotRMSD(self,input,range,step,opt)
            %Plot being sepreated by range and step args
            %Choosing title and label for cases is optional
            %             the first row should be the baseline of the signals
            %             the output Y is the matrix of different RMSD in different
            %             samples and ranges
            %           this function also plot the RMSD results
            arguments
                self
                input
                range
                step
                opt.Title=9999
                opt.Label=999
                opt.plot=true
            end
            %input: matrix of signal that each column is a indiviual signal
            %range: range that RMSD will calculate with respect to step of range
            %step:Each signal devide bt this step in range
            %fs:sampling range
            input=input';
            L1=size(input,2);
            L=size(input,1);
            input=real(input);
            %             f = self.fs*(0:floor(L)-1)/(2*L);
            s=1;
            Y=zeros((L1-1),length(range(1):step:(range(2)-step)));
            Legend=cell((L1-1),1);
            for j=range(1):step:(range(2)-step)
                
                index=find(self.Fnew>j&self.Fnew<(j+step));
                self.RMSD_FirstArg=input(index,1);
                for i=2:L1
                    self.ZeRMSD=input(index,i);
                    Y(i-1,s)=self.RMSD();
                    Legend(i-1)={['case',int2str(i)]};
                end
                s=s+1;
            end
            if(opt.plot)
                for i=1:size(Y,2)
                    cat=(range(1)+i*step)/1000;
                    label(i)={[int2str(cat-step/1000) '-' int2str(cat) 'KHz']};
                end
                X = categorical(label);
                X = reordercats(X,label);
                bar(X,100*Y')
                if(~(isnumeric(opt.Label)))
                    Legend=opt.Label;
                end
                legend(Legend)
                ylabel('RMSD(%)')
                
                if(~(isnumeric(opt.Title)))
                    title(opt.Title)
                end
            end
        end
        %% Finding sample in data(cocatated version)
        function num=RSampleVal(~,sen,loc,sample,Location,Sensor)
            %Location and Sensor are cocatated matrix of of them
            %sen,loc and sample are sensor location and sample witch we
            %want to find the data number in cocatated version of data
            isen=(Sensor==sen);
            iLoc=(Location==loc);
            index=find(isen.*iLoc);
            num=index(sample);
        end
        
        %%
        function dict=LoadMatData(~,pathfile,sev,win)
            %loading data from matrix of data and saving them in
            %dictionary dict respect to their severity.
            
            dict=containers.Map();
            for i=1:length(sev)
                FFile=[pathfile,'/',sev{i},'/',win];
                load(FFile)
                dict(['Adm',sev{i}])=real(Admittance);
                dict(['Sev',sev{i}])=real(Severity);
                dict(['Loc',sev{i}])=real(Location);
                dict(['Sen',sev{i}])=real(Sensor);
                if(i==1)
                    dict('Freq')=Freq;
                end
            end
        end
        %%
        function [num, IndTot]=RSampleVal2(~,sen,loc,sev,Location,Sensor,Severity,sample)
            % Finding the number of samples in the data matrix with the given detail.
            % The 'loc' and 'sev' arguments should be vectors representing the desired ...
            % location and severity.
            % If the number of input arguments 7, the function returns the...
            % range of samples that match the given detail.
            % Otherwise, it only returns the 'sample'th element of the given detail.
            IndTot=[];
            num=[];
            for i=1:length(sev)
                %                 for j=1:length(loc)
                isen=(Sensor==sen);
                iLoc=(Location==loc(i));
                iSev=(Severity==sev(i));
                %     iAdmittance=(Ad
                ind4=find(isen.*iLoc.*iSev);
                ind3=find(isen);
                if(~isempty(ind4))
                    if(nargin==7)
                        num=[ind4(1),ind4(end);num];
                        IndTot=[IndTot ind4'];
                    else
                        num=ind4(sample);
                        isen2=(Sensor==sen);
                        iLoc2=(Location==0);
                        iSev2=(Severity==0);
                        ind3=find(isen2.*iLoc2.*iSev2);
                        IndTot=num-ind3(1)+1;
                    end
                end
                %                 end
            end
        end
        %%
        function mat=RMSDMat(self,sev,loc,sen,dict,sevcell,sample)
            % The RMSDMat function computes the matrix required for plotRMSD.
            %The function takes the following input arguments:
            % 1. self: the object being analyzed.
            % 2. sev: the severity of the tests.
            % 3. loc: the location of the tests.
            % 4. sen: the sensor number of the tests.
            % 5. dict: a dictionary holding the data values.
            % 6. sevcell: the names of the keys that the function should read.
            % 7. sample: the number of test samples taken for each
            % specified tests with same modes.
            %   if 'mean' is inserted, the function will average from all similar test samples
            %	if '0.5-mean' is inserted, the function will average from half of all
            %   similar test samples
            % Note: There are two choices for input arguments. ...
            %The inputs can have a fixed number of severities, locations,...
            %and sensors, and a 2D vector for samples that represents a range of samples.
            %Alternatively, the inputs can be vectors where each element denotes a specific sample.
            %The function returns a matrix, which is required for plotRMSD.
            MInd=0;
            for i=1:length(sevcell)
                Admittance=dict(['Adm',sevcell{i}]);
                Severity=dict(['Sev',sevcell{i}]);
                Location=dict(['Loc',sevcell{i}]);
                Sensor=dict(['Sen',sevcell{i}]);
                
                if(ischar(sample))
                    switch(sample)
                        case('mean')
                            num=self.RSampleVal2(sen,loc(i),sev(i),Location,Sensor,Severity);
                            Inx=num(1):num(2);
                            mat(i,:)=mean(Admittance(Inx,:));
                            
                        case('0.5-mean')
                            num=self.RSampleVal2(sen,loc(i),sev(i),Location,Sensor,Severity);
                            Inx=num(1):num(2);
                            Inx1=Inx(1:floor(end/2));
                            mat((2*i-1),:)=mean(Admittance(Inx1,:));
                            Inx2=Inx((floor(end/2)+1):end);
                            mat(2*i,:)=mean(Admittance(Inx2,:));
                    end
                    
                else
                    if(and(length(sev)~=1,length(loc)~=1))
                        num=self.RSampleVal2(sen,loc(i),sev(i),Location,Sensor,Severity,...
                            sample(i));
                        mat(i,:)=Admittance(num,:);
                    elseif(length(sample)~=1)
                        num(1)=self.RSampleVal2(sen,loc,sev,Location,Sensor,Severity,sample(1));
                        num(2)=self.RSampleVal2(sen,loc,sev,Location,Sensor,Severity,sample(2));
                        MIndNew=((num(2)-num(1))+1)+MInd;
                        mat((MInd+1):MIndNew,:)=Admittance(num(1):num(2),:);
                        MInd=MIndNew;
                    else
                        num=self.RSampleVal2(sen,loc,sev,Location,Sensor,Severity,sample);
                        mat=Admittance(num,:);
                    end
                end
            end
        end
        %%
        function [out1,out2,out3,out4]=SortData(~,in1,in2,in3,in4,dim,MaxNum,condition)
            [p,sortInd]=sort(in1,dim,condition);
            out1=p(1:MaxNum);
            out2=in2(sortInd(1:MaxNum));
            out3=in3(sortInd(1:MaxNum));
            out4=in4(sortInd(1:MaxNum));
        end
        %%
        function [values,rows,columns]=FindMax(~,mat,num)
            abs_A =mat(:); % vectorize values
            [sorted_A, idx] = sort(abs_A,'descend'); % sort in descending order
            values = sorted_A(1:num); % select top num absolute values
            top_10_idx = idx(1:num); % select indices corresponding to top num absolute values
            [rows, columns] = ind2sub(size(mat), top_10_idx); % convert linear indices to subscripts
            
        end
        %%
        function [rows,columns,values]=remove(~,i,rows,columns,values)
            rows(i:(end-1))=rows((i+1):end);
            columns(i:(end-1))=columns((i+1):end);
            values(i:(end-1))=values((i+1):end);
            rows(end)=0;
            columns(end)=0;
            values(end)=0;
        end
        %%
        function result=divideData(self,dict,loc,sen,sev,optRMSD)
            % Finding data for training testing and validation set
            % result: a dictionary of vars in desired range
            % dict:input variables in dict format
            % loc:location:MUST include contnious range
            % sen:choose sensor
            % sev:choose whether function consider multiple severity
            % if the range and step args have value the function gets RMSD
            % if the mean be true the function gets the mean of each case
            %%
            arguments
                %%
                self
                dict
                loc
                sen
                sev
                optRMSD.range='a';
                optRMSD.step='a';
                optRMSD.mean=false
                optRMSD.concat=false
            end
            
            
            jj=1;
            for i=sev
                strsev=int2str(i);
                sevcell=[strsev,'mag'];
                Admittance=dict(['Adm',sevcell]);
                Severity=dict(['Sev',sevcell]);
                Location=dict(['Loc',sevcell]);
                Sensor=dict(['Sen',sevcell]);
                Label=dict(['Label',sevcell]);
                SevVec=(loc~=0)*i;
                %%
                for s=sen
                    [~,Index]=RSampleVal2(self,s,loc,...
                        SevVec,Location,Sensor,Severity);
                    %%
                    admit=Admittance(Index,:);
                    if(optRMSD.mean)
                        admit=mean(Admittance(Index,:));
                    end
                    if(isnumeric(optRMSD.range))
                        Mat=RMSDMat(self,0,0,s,dict,{sevcell},'mean');
                        MatF=[Mat;Admittance(Index,:)];
                        Y=plotRMSD(self,MatF,optRMSD.range,optRMSD.step,plot=false);
                        admit=Y;
                    end
                    %%
                    if(optRMSD.mean)
                        if(jj==1)
                            TLoc=[Location(Index(1))];
                            TSen=[Sensor(Index(1))];
                            
                            TAdm=admit;
                            TLabel=[Label(Index(1),:)];
                            TSev=[Severity(Index(1),:)];
                            jj=2;
                        else
                            TLoc=[TLoc;Location(Index(1))];
                            TSen=[TSen;Sensor(Index(1))];
                            TAdm=[TAdm;admit];
                            TLabel=[TLabel;Label(Index(1),:)];
                            TSev=[TSev;Severity(Index(1),:)];
                        end
                    else
                        %%
                        
                        if(jj==1)
                            TLoc=[Location(Index)];
                            TSen=[Sensor(Index)];
                            
                            TAdm=admit;
                            TLabel=[Label(Index,:)];
                            TSev=[Severity(Index,:)];
                            jj=2;
                        else
                            TLoc=[TLoc;Location(Index)];
                            TSen=[TSen;Sensor(Index)];
                            TAdm=[TAdm;admit];
                            TLabel=[TLabel;Label(Index,:)];
                            TSev=[TSev;Severity(Index,:)];
                        end
                    end
                end
                
            end
            
            %%
            result=containers.Map();
            result('Adm')=TAdm;
            result('Sev')=TSev;
            result('Loc')=TLoc;
            result('Sen')=TSen;
            result('Label')= TLabel;
            
            
        end
        
        %%
        function res=ImpRange(~,Impedance, Freq, rangefreq)
            inx=((rangefreq(1)<Freq).*(Freq<rangefreq(2)));
            res=Impedance(inx);
            
        end
        
        %%
        function result=FeatureConcat(self,dict,loc,sen,sev,opt)
            arguments
                %%
                self
                dict
                loc
                sen
                sev
                opt.tensor='rectangular'
            end
            %%
            %Finding data for training testing and validation set
            %             result: a dictionary of vars in desired range
            % dict:input variables in dict format
            % loc:location:MUST include contnious range
            % sen:choose sensor
            % sev:choose whether function consider multiple severity
            %%
%             ll=true;
            
            Admittance=dict('Adm');
            Severity=dict('Sev');
            Location=dict('Loc');
            Sensor=dict('Sen');
            Label=dict('Label');
            %                 SevVec=(loc~=i)*i;
            %%
%             Finding indices for zeros preallocations
            index=cell(length(sev),length(loc),length(sen));
            Lsamples=cell(length(sev),length(loc));
            Lsamples(:,:)={2e6};
            RowSum=0;
            for ii=1:length(sev)
                SevVec=ones(size(loc))*sev(ii);
                for Lo=1:length(loc)
                    for s=sen
                        [~,ind89]=RSampleVal2(self,s,loc(Lo),...
                            SevVec(Lo),Location,Sensor,Severity);
                        index(ii,Lo,s)={ind89};
                        if(length(ind89)<Lsamples{ii,Lo})
                            Lsamples(ii,Lo)={length(ind89)};
                        end
                        
                    end
                    RowSum=RowSum+Lsamples{ii,Lo};
                end
            end
            ColNum=size(Admittance,2);
            colsev=size(Severity,2);
            colloc=size(Location,2);
%             colsen=size(Sensor,2);
            collabel=size(Label,2);
            FAdm=zeros(1,4*ColNum);
            TAdm=zeros(RowSum,4*ColNum);TSev=zeros(RowSum,colsev);TLoc=zeros(RowSum,colloc);
            TLabel=zeros(RowSum,collabel);
            %%
            GolRow=1;
            for i=1:length(sev)
                %                 strsev=int2str(i);
                %                 sevcell=[strsev,'mag'];
                %%
                for Lo=1:length(loc)
%                     index=containers.Map();
                    %%
                    
                    
                    %%
                    for sam=1:Lsamples{i,Lo} %loop over samples that have same properties except...
%                         sensor number. max of iterator is the min of
%                         all sensor samples
%                         jj=true;
                        for sennum=sen  %loop over all sensors for each sample
                            sampll=index{i,Lo,sennum};
                            admit=Admittance(sampll(sam),:);
                            ColRange=(((sennum-1)*ColNum)+1):sennum*ColNum;
                            FAdm(1,ColRange)=admit;
%                             %                             jj=false;
%                             %                             if(jj)
%                             %                                 FAdm=admit;
%                             %                                 jj=false;
%                             %                             else
%                             %                                 FAdm=[FAdm,admit];
%                             %                             end
                        end
                        llabel=Label(sampll(sam),:);
                        
%                         if(ll)
%                             TAdm= FAdm;
%                             TLabel=llabel;
%                             TSev=Severity(sampll(sam),:);
%                             TLoc=Location(sampll(sam));
%                             ll=false;
%                         else
%                             Lsam=Lsamples{i,Lo};
                            rows=GolRow;
                            TAdm(rows,:)=FAdm;
                            TLabel(rows,:)=llabel;
                            TLoc(rows,:)=Location(sampll(sam));
                            TSev(rows,:)=Severity(sampll(sam));
                            GolRow=GolRow+1;
%                         end
                    end
                end
                %%
            end
            result=containers.Map();
            switch(opt.tensor)  
                case('cubic') %If the tensor arg be cubin, the function return 3D admittance tensor
                    samples=size(TAdm,1);
                    features=size(TAdm,2);
                    sen_num=length(sen);
                    width=features/sen_num;
                    tensor=zeros(samples,sen_num,width);
                    start_inx = (1:width:features);
                    for s=1:sen_num
                        
                        col_range=(start_inx(s):(start_inx(s)+width-1));
                        tensor(:,s,:)=TAdm(:,col_range);
                    end
                    result('Adm')=tensor;
                otherwise
                    result('Adm')=TAdm;
            end
            result('Sev')=TSev;
            result('Loc')=TLoc;
            result('Sen')=ones(size(TLoc));
            result('Label')= TLabel;
        end
        %%
        function res=sortingdict(~,peaks,FeatureSize)
            %the aim is for sorting the data with consideration of points
            %of peaks
%             peaks=dict1('Adm');
            len=size(peaks,2);
            res=zeros(size(peaks));
            s=1;
            for j=0:(FeatureSize-1)
                for i=1:FeatureSize:len
                    res(:,s)=peaks(:,(i+j));
                    s=s+1;
                end
            end
        end
        %%
        function [NormTrain,MeanTrain,VarTrain,NormVal]=normalizing(~,inputTrain,inputVal)
            MeanTrain=mean(inputTrain);
            VarTrain=std(inputTrain);
            NormTrain=(inputTrain-MeanTrain)./(VarTrain);
            if(nargin==3)
            NormVal=(inputVal-MeanTrain)./(VarTrain);
            end
            
        end
        %%
        function [NormTrain,MeanTrain,NormVal]=meancentric(~,inputTrain,inputVal)
            MeanTrain=mean(inputTrain);
            %             VarTrain=std(inputTrain);
            NormTrain=(inputTrain-MeanTrain);
            if(nargin==3)
                NormVal=(inputVal-MeanTrain);
            end
            
        end
        %            dict2=dict1;
        % Function for downsampling with maximum
        %%
        function [downsampled_data_max,fnew] = downsample_max(~,data, downsampling_factor,freq)
            num_samples = size(data, 1);
            num_sensors = size(data, 2);
            num_fetures = size(data, 3);
            width=floor(num_fetures/downsampling_factor);
            
            downsampled_data_max = zeros(num_samples,num_sensors,width);
            indFreq=zeros(width,1);
            itr=1;
            
            for i = 1:downsampling_factor:num_fetures
                start_idx = i;
                end_idx = min(i + downsampling_factor - 1, num_fetures);
                [downsampled_data_max(:,:,itr),Ind]= max(data(:,:,(start_idx:end_idx)),[],3);
                indFreq(itr)=start_idx+Ind-1;
                itr=itr+1;
            end
            fnew=freq(indFreq);
        end
        %%
        function [downsampled_data_max,fnew] = downsample_mean(~,data, downsampling_factor,freq)
            num_samples = size(data, 1);
            num_sensors = size(data, 2);
            num_fetures = size(data, 3);
            width=floor(num_fetures/downsampling_factor);
            
            downsampled_data_max = zeros(num_samples,num_sensors,width);
%             indFreq=zeros(width,1);
            itr=1;
            
            for i = 1:downsampling_factor:num_fetures
                start_idx = i;
                end_idx = min(i + downsampling_factor - 1, num_fetures);
                [downsampled_data_max(:,:,itr)]= mean(data(:,:,(start_idx:end_idx)),3);
                itr=itr+1;
            end
            fnew=freq(1:downsampling_factor:num_fetures);
        end
        
    end
end
%%