clc; clear; close all
%%  Program explenation:
%    This program compute admittance of files with .TXT format and path witch saved in "FilePath".
%    Program saving admittance signals in Ze variable and averaging  them every 'NSample' then
%    saving in Zem. all functions defined as methods in Class Imp. This calss also contain
%    visulization methods that some of them are commented throgh main program.
obj=Imp();

%% Insert data file specification
[obj.start_ind1,obj.end_ind1]=deal(1,500);        %loop over samples indexx
[obj.start_ind2,obj.end_ind2]=deal(0,23);        %loop over file of damage cases inde
% [obj.start_ind3,obj.end_ind3]=deal(1,3);        %loop over file of damage cases inde
FilePath='C:\data\Newfolder\6magnet_test'; %Last folder that all data exists on that
FilePath4Ini='C:\data\Newfolder\6magnet_test\sensor2\0\run1.txt';   %file path for initilization
parentfolder={'C:\data\Newfolder\6magRes'};     %parent file for saving results
%% Initialization(this step is for decreasing computation cost)
[obj,In,Out]=obj.initial(FilePath4Ini);
AnalisysMode={'kaiser','rec','recCut','tukey','tukeyCut'};                              %Intering different analisys modes.
% parentfolder=obj.createfolder(parentfolder,AnalisysMode); %creating folder for different modes.
% SensorName={'sensor1'};                                   %intering different sensor name
% parentfolder=obj.createfolder(parentfolder,SensorName);   %creating folder for different sensor name
% ZeTotal=zeros([size(Ze) 2]);
%%
% for s=obj.start_ind3:obj.end_ind3
for s=2:4 
    %%
    FilePath1=[FilePath,'\sensor',int2str(s),'\'];
    %     for sever=obj.start_ind2:obj.end_ind2
    %         %%
    %         FilePath1=[FilePath2,'\',int2str(sever),'\'];
    %         %%
    for j=obj.start_ind2:obj.end_ind2
        
        %     mkdir ['D:\Newfolder(5)\New folder' '] 'int2str(j)'
        datafile1=[FilePath1,int2str(j)];
        obj.changename(datafile1)
        obj=obj.CreateFolder(AnalisysMode,parentfolder,datafile1,FilePath);
        
        inx1=1;
        inx2=1; %Zem(mean of impedance) counter
        %%
        for i=obj.start_ind1:obj.end_ind1
            obj.SampleFile=[datafile1,'\run',int2str(i),'.txt'];
            [A,B]=obj.LoadDataIni();
            %         [obj,~,Anew,Bnew]=obj.FindSignal(A,B);
            %       Note:for plotting signals use following FindSignal Function instead of
            %            pervious function
            [obj,tnew,Anew,Bnew]=obj.FindSignal(A,B);
            if(obj.AbnormalSignal)
                obj.AbnormalSignal=false;
                continue;
            end
            FileNum=1; %index for self.ResultFolder folder
            if(inx1==1)
                ArrObj(1:length(AnalisysMode))=obj;
            end
            %% First mode
            [ArrObj(1),FileNum]=ArrObj(1).analisys(Anew,Bnew,inx1,inx2,FileNum,...
                A=A,B=B,tnew=tnew,window='kaiser',cropping=true,RMSD=true,Sfile=obj.SampleFile);
            if(ArrObj(1).AbnormalSignal)
                ArrObj(1).AbnormalSignal=false;
                continue;
            end
            %% Second mode
            [ArrObj(2),FileNum]=ArrObj(2).analisys(Anew,Bnew,inx1,inx2,FileNum,...
                A=A,B=B,tnew=tnew,visualization=false...
                ,window='rec',cropping=false,In=In,Out=Out,InOut=true);
            %% Third mode
            [ArrObj(3),FileNum]=ArrObj(3).analisys(Anew,Bnew,inx1,inx2,FileNum,...
                A=A,B=B,tnew=tnew,window='rec',cropping=true);
            %% Forth mode
            [ArrObj(4),FileNum]=ArrObj(4).analisys(Anew,Bnew,inx1,inx2,FileNum,...
                A=A,B=B,tnew=tnew,window='tukey',cropping=false);
            %% fifth mode     
            %in the last mode inx must be updated
            [ArrObj(5),FileNum,inx1,inx2]=ArrObj(5).analisys(Anew,Bnew,inx1,inx2,FileNum,...
                A=A,B=B,tnew=tnew,window='tukey',cropping=true);
            
            
        end
        running=['sensor ',int2str(s),'File ',int2str(j)]
    end
    
end
% end
%% Results visualization
% plot(obj.Fnew,Zem)
%save("impedannce.mat",'obj','Zem','-append')
% obj.plotZe(real(Zem))
% hold on
% plotZe(imag(Ze),fs)
% plotZe(abs(Ze),fs)
% legend("Real","Imag","Abs")
%%
