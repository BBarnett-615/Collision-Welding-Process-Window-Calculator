%Collision-Welding-Calculator
%(C) Blake Barnett 2023
%Collision Welding Process Parameter Calculator for similar and dissimilar
%welding.
%
%This collision welding parameter calculation tool is licensed under a 
%Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License 
%CC BY-NC-SA 4.0 [https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode]
%
%For questions and licensing permissions, contact: 
%barnett.615@osu.edu

classdef Barnett_Weld_Window_v1_2_7_23_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        A1CheckBox                      matlab.ui.control.CheckBox
        A2CheckBox                      matlab.ui.control.CheckBox
        S1CheckBox                      matlab.ui.control.CheckBox
        S2CheckBox                      matlab.ui.control.CheckBox
        E1CheckBox                      matlab.ui.control.CheckBox
        E2CheckBox                      matlab.ui.control.CheckBox
        FlyerMaterialDropDownLabel      matlab.ui.control.Label
        FlyerMaterialDropDown           matlab.ui.control.DropDown
        TargetMaterialDropDownLabel     matlab.ui.control.Label
        TargetMaterialDropDown          matlab.ui.control.DropDown
        FlyerThicknessmmEditFieldLabel  matlab.ui.control.Label
        FlyerThicknessmmEditField       matlab.ui.control.NumericEditField
        TargetThicknessmmEditFieldLabel  matlab.ui.control.Label
        TargetThicknessmmEditField      matlab.ui.control.NumericEditField
        WeldPairSwitch                  matlab.ui.control.RockerSwitch
        VimsEditFieldLabel              matlab.ui.control.Label
        VimsEditField                   matlab.ui.control.NumericEditField
        ImpactAngleDegreesEditFieldLabel  matlab.ui.control.Label
        ImpactAngleDegreesEditField     matlab.ui.control.NumericEditField
        HotspotIntervalmicronsSliderLabel  matlab.ui.control.Label
        HotspotIntervalmicronsSlider    matlab.ui.control.Slider
        WeldLengthmmEditFieldLabel      matlab.ui.control.Label
        WeldLengthmmEditField           matlab.ui.control.NumericEditField
        UITable                         matlab.ui.control.Table
        CalculateButton                 matlab.ui.control.Button
        dxEdit                          matlab.ui.control.NumericEditField
        NFitLabel                       matlab.ui.control.Label
        A1NEdit                         matlab.ui.control.NumericEditField
        A2NEdit                         matlab.ui.control.NumericEditField
        S1NEdit                         matlab.ui.control.NumericEditField
        S2NEdit                         matlab.ui.control.NumericEditField
        E1NEdit                         matlab.ui.control.NumericEditField
        E2NEdit                         matlab.ui.control.NumericEditField
        UpperLimitsLabel                matlab.ui.control.Label
        WeldingParametersLabel          matlab.ui.control.Label
        UIAxes                          matlab.ui.control.UIAxes
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxes3                         matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        T % Description
        Mat1 % Flyer Material String
        Mat2 % Target Material String
        beta %range of impact angles
        Vi %range of impact velocities
        h1 % Flyer thickness (mm)
        h2 % Target thickenss (mm)
        L % Weld length (mm)
        Vexp = 25; %Impact Velocity in the specified experiment
        Betaexp %Impact Angle in the specified experiment
        Up1exp %Experimental flyer particle velocity
        Up2exp %Experimental target particle velocity
        Us1exp %Experimental flyer shock velocity
        Us2exp %Experimental target shock velocity
        Pexp %Experimental Pressure
        ShockTable %Contains shock calculations for weld pair
        dx=50 %Hot Spot Spacing (microns)
        Vtable %Table with columns of Vi,beta,Vc sorted in order from least to greatest Vc
        ULimits %Initialize Upper Limit Selections
        Flyvals %Reporting Values for Flyer
        Targvals %Reporting Values for Target
        I %Index in ShockTable of Experimental Impact Value
        Calc=0; %Tracks whether the Calculate Button has been pushed
        %Upper Limit plot handles and fitting parameters for the Welding Window graph
        hA1
        A1N=.11;
        hA2
        A2N=.11;
        hS1
        S1N=1;
        hS2
        S2N=1;
        hE1
        E1N=1;
        hE2
        E2N=1;
        hZ1
        hZ2
        FS
        TS
    end
    
    methods (Access = private)
        
        function ShockTable = ShockCalcs(app,T,Mat1,Mat2,Vtable)
            C1 = T{Mat1,'C0'};
            S1 = T{Mat1,'S'};
            rho1 = T{Mat1,'Rho'};
            C2 = T{Mat2,'C0'};
            S2 = T{Mat2,'S'};
            rho2 = T{Mat2,'Rho'};
            ImpactVs=Vtable(:,1).*cosd(Vtable(:,2));
            
            Up1=zeros(numel(ImpactVs),1);
            Up2=zeros(numel(ImpactVs),1);
            Us1 = zeros(numel(ImpactVs),1);
            Us2 = zeros(numel(ImpactVs),1);
            P = zeros(numel(ImpactVs),1);
            
            for i = 1:numel(ImpactVs)
                a = (rho2*S2)-(rho1*S1);
                b = (rho2*C2)+(rho1*C1)+(2*rho1*S1*ImpactVs(i));
                c = -rho1*((C1*ImpactVs(i))+(S1*(ImpactVs(i))^2));
                d = b^2-(4*a*c);
        
                if d > 0
                    root1 = ((-1*b)+sqrt(d))/(2*a);
                    root2 = ((-1*b)-sqrt(d))/(2*a);
                elseif d==0
                    root1 = -b/(2*a);
                    root2 = NaN;
                else
                    root1 = NaN;
                    root2 = NaN;
                end
        
                if a == 0
                    Up2(i) = ImpactVs(i)/2;
                elseif 0<root1 && root1<ImpactVs(i)
                    Up2(i) = root1;
                else
                    Up2(i) = root2;
                end
                Up1(i)=ImpactVs(i)-Up2(i);
                P(i)=(rho2*(C2+S2*Up2(i))*Up2(i));
                Us1(i)=C1+S1*Up1(i);
                Us2(i)=C2+S2*Up2(i);
            end 
            ShockTable = [Vtable Up1 Up2 Us1 Us2 P];
        end
        
        function [VminB,JFilter] = MinJet(app,Mat1,Mat2)
             Vc = app.Vtable(:,3);
             V1 = sqrt(app.T{Mat1,'UTS'}/app.T{Mat1,'Rho'});
             V1B = 2*real(asind(V1./(2.*Vc)));
             V2 = sqrt(app.T{Mat2,'UTS'}/app.T{Mat2,'Rho'});
             V2B = 2*real(asind(V2./(2.*Vc)));
             VminB = zeros(numel(V1B),1);
             for i=1:numel(VminB)
                 VminB(i) = max(V1B(i),V2B(i));
             end
             
              Pa1 = zeros(numel(Vc),numel(app.beta));
              Pa2 = zeros(numel(Vc),numel(app.beta));
              Jlim = zeros(numel(Vc),1);
              
              for i =1:numel(Vc)
                for j = 1:numel(app.beta)
                    Up1 = app.ShockTable(i,4);
                    Up2 = app.ShockTable(i,5);
                    Us1 = app.ShockTable(i,6);
                    Us2 = app.ShockTable(i,7);
           
                    Pa1(i,j)= atand(Up1*sqrt(Vc(i)^2-Us1^2)/(Vc(i)^2-(Up1*Us1)));
                    Pa2(i,j)= atand(Up2*sqrt(Vc(i)^2-Us2^2)/(Vc(i)^2-(Up2*Us2)));
                end
              end
              ac1 = real(max(Pa1,[],2));
              ac2 = real(max(Pa2,[],2));
              start1=1;
              start2=1;
              
              for i = numel(Vc):-1:1
                if ac1(i)<=0
                    start1=i;
                    break
                end
              end
              
              for i = numel(Vc):-1:1
                if ac2(i)<=0
                    start2=i;
                    break
                end
              end

              start = min(start1,start2);
              for i = start:numel(Jlim)
                  Jlim(i) = max(ac1(i),ac2(i));
              end
              
              [C,~] = envelope(Jlim,5,'peak');
              for i = 1:numel(Jlim)
                  if Jlim(i)>0
                      EnvelopeStart = i;
                      break
                  end
              end
              JFilter = [Vc(EnvelopeStart:numel(Vc)) C(EnvelopeStart:numel(C))];
        end
        
        function Vtable = Velocities(app,Vi,beta)
            Vc = zeros(numel(Vi),numel(beta));
            Velocities=[];
            for i = 1:numel(Vi)
                for j = 1:numel(beta)
                    Vc(i,j)=Vi(i)./(2.*sind(beta(j)));
                    Velocities = [Velocities;Vi(i) beta(j) Vc(i,j)];   
                end
            end
            Vtable = sortrows(Velocities,3,"ascend");
        end
        
        function A1 = A1Calc(app,T,Mat1,Mat2,h1,h2,N,Vtable)
             Wh1 = .1*h1;
             Vc=Vtable(:,3);
             WC1=.01*T{Mat1,'C0'};
             WTm1=T{Mat1,'Tm'};
             Wk1=.01*T{Mat1,'k'};
             Wcp1=1e9*T{Mat1,'cp'};
             Wrho1=.001*T{Mat1,'Rho'};
             N=double(N);
             D1=((1/N)*(WTm1*WC1)^.5) * ((Wk1*Wcp1*WC1)/(Wrho1*Wh1))^.25;
             delta1=zeros(numel(Vc),1);
             Wh2 = .1*h2;
             WC2=.01*T{Mat2,'C0'};
             WTm2=T{Mat2,'Tm'};
             Wk2=.01*T{Mat2,'k'};
             Wcp2=1e9*T{Mat2,'cp'};
             Wrho2=.001*T{Mat2,'Rho'};
             D2=(1/N)*(WTm2*WC2)^.5 * ((Wk2*Wcp2*WC2)/(Wrho2*Wh2))^.25;
             delta2=zeros(numel(Vc),1);
             delta=zeros(numel(Vc),1);
             for i=1:numel(Vc)
                 delta1(i)=real(asind(D1/(Vc(i)^2)));
                 delta2(i)=real(asind(D2/(Vc(i)^2)));
                 delta(i)=min(delta1(i),delta2(i));
             end
             A1=[Vc delta];
        end
        
        function A2 = A2Calc(app,T,Mat1,Mat2,h1,h2,N,Vtable)
             Wh1 = .1*h1;
             N=double(N);
             Vc=Vtable(:,3);
             WC1=.01*T{Mat1,'C0'};
             WTm=min(T{Mat1,'Tm'},T{Mat2,'Tm'});
             Wk =(2*(.01*T{Mat1,'k'})*.01*T{Mat2,'k'})/(.01*T{Mat1,'k'}+.01*T{Mat2,'k'});
             Wcp= (2*(1e9*T{Mat1,'cp'})*1e9*T{Mat2,'cp'})/(1e9*T{Mat1,'cp'}+1e9*T{Mat2,'cp'});
             Wrho1=.001*T{Mat1,'Rho'};
             Wrho2=.001*T{Mat2,'Rho'};
             D1=(1/N)*(WTm*WC1)^.5 * ((Wk*Wcp*WC1)/(Wrho1*Wh1))^.25;
             delta1=zeros(numel(Vc),1);
             Wh2 = .1*h2;
             WC2=.01*T{Mat2,'C0'};
             D2=(1/N)*(WTm*WC2)^.5 * ((Wk*Wcp*WC2)/(Wrho2*Wh2))^.25;
             delta2=zeros(numel(Vc),1);
             delta=zeros(numel(Vc),1);
             for i =1:numel(Vc)
                 delta1(i)=real(asind(D1/(Vc(i)^2)));
                 delta2(i)=real(asind(D2/(Vc(i)^2)));
                 delta(i)=min(delta1(i),delta2(i));
             end
             A2=[Vc delta];
        end
        
        function S1 = S1Calc(app,T,Mat1,Mat2,h1,N,ShockTable)
            h1=.001*h1;
            f = 0.9;
            Vc=ShockTable(:,3);
            rho1 = T{Mat1,'Rho'};
            Tm1 = T{Mat1,'Tm'};
            k1 = T{Mat1,'k'};
            alpha1=T{Mat1,'alpha'};
            Tm2= T{Mat2,'Tm'};
            k2= T{Mat2,'k'};
            alpha2=T{Mat2,'alpha'};
            k =(2*k1*k2)/(k1+k2);
            alpha = (2*alpha1*alpha2)/(alpha1+alpha2);
            Tm = min([Tm1 Tm2]);
            
            STable1=zeros(height(ShockTable),2);
            for i=1:height(STable1)
                delta1=(Vc(i)/ShockTable(i,6))*cosd(ShockTable(i,2))-tand(ShockTable(i,2));
                STable1(i,1)=Vc(i);
                STable1(i,2)=real(asind((((Tm*k)/(f*rho1*h1^2*Vc(i)*delta1))^.5*(((pi*Vc(i)*h1)/alpha)*delta1)^.25)/(N*Vc(i))));
            end
            
            SRaw=[STable1(:,1) STable1(:,2)];
            
            SFilter=zeros(size(SRaw));
            SFilter(:,1)=SRaw(:,1);
            
            for i = 1:height(SRaw)
                if SRaw(i,2)<=0
                    SFilter(i,2) = 0;
                else
                    SFilter(i,2) = SRaw(i,2);
                end
            end
                        
            SFilter2=zeros(size(SFilter));
            SFilter2(:,1)=SFilter(:,1);
            bench=180;
            for i = 1:height(SRaw)
                if SFilter(i,2)<=bench&&SFilter(i,2)~=0
                    SFilter2(i,2)=SFilter(i,2);
                    bench=SFilter2(i,2);
                else
                    SFilter2(i,2)=bench;
                end
            end            
            S1=[SFilter2(:,1) SFilter2(:,2)];
        end
        
        function S2 = S2Calc(app,T,Mat1,Mat2,h1,h2,ShockTable,N,L,dx)
            h1=.001*h1;
            h2=.001*h2;
            dx=dx*1e-6;
            f = 0.9;
            L = .001*L;
            Vc=ShockTable(:,3);
            rho1 = T{Mat1,'Rho'};
            rho2 = T{Mat2,'Rho'};
            Tm1 = T{Mat1,'Tm'};
            k1 = T{Mat1,'k'};
            cp1=T{Mat1,'cp'};
            Tm2= T{Mat2,'Tm'};
            k2= T{Mat2,'k'};
            cp2=T{Mat2,'cp'};
            Tm = min([Tm1 Tm2]);
            t1=(2*h1)/(app.Us1exp);
            t2=(2*h2)/(app.Us2exp);
            tf=min(t1,t2);
            A = ((rho1*cp1*k1)^.5+(rho2*cp2*k2)^.5);
            B = N*((4*tf*Tm^2*A^2*L^2*pi)/(f^2*rho1^2*h1^2*dx^2))^.25;
            S2=zeros(numel(Vc),2);
            for i = 1:numel(Vc)
                S2(i,1)=Vc(i);
                S2(i,2)=2*asind(B/(2*Vc(i)));
            end
        end
        
        function E1 = E1Calc(app,T,Mat1,~,h1,h2,Vtable,N)
            h1=.001*h1;
            h2=.001*h2;
            Vc=Vtable(:,3);
            E1=[Vc zeros(numel(Vc),1)];
            G1 = T{Mat1,'G'};
            Tm1=T{Mat1,'Tm'};
            k1 = T{Mat1,'k'};
            alpha1=T{Mat1,'alpha'};
            rho1=T{Mat1,'Rho'};
            a = (Tm1*k1)/alpha1;
            b = rho1*(h2/(h1+h2))*sqrt(h1);
            for i=1:height(E1)
                E1(i,2)=N*real(asind((11.8*(Vc(i)^(-5/4))*sqrt((a)/(b)))*(0.5+0.66*(rho1*Vc(i)^2/G1))^.25));
            end
        end
        
        function E2 = E2Calc(app,T,Mat1,Mat2,h1,h2,Vtable,N)
            h1=.001*h1;
            h2=.001*h2;
            Vc=Vtable(:,3);
            E2=[Vc zeros(numel(Vc),1)];
            G1 = T{Mat1,'G'};
            Tm1=T{Mat1,'Tm'};
            Tm2=T{Mat2,'Tm'};
            Tm=min(Tm1,Tm2);
            k1 = T{Mat1,'k'};
            k2 = T{Mat2,'k'};
            k = (2*k1*k2)/(k1+k2);
            alpha1=T{Mat1,'alpha'};
            alpha2=T{Mat2,'alpha'};
            alpha = (2*alpha1*alpha2)/(alpha1+alpha2);
            rho1=T{Mat1,'Rho'};
            a = (Tm*k)/alpha;
            b = rho1*(h2/(h1+h2))*sqrt(h1);
            for i=1:height(E2)
                E2(i,2)=N*real(asind((11.8*(Vc(i)^(-5/4))*sqrt((a)/(b)))*(0.5+0.66*(rho1*Vc(i)^2/G1))^.25));
            end
        end
        
        function [z1,z2,t,solidtimes,maxz,FShock,TShock]=ztCalc(app,T,Mat1,Mat2,h1,h2,Vexp,L,dx)
            h1=.001*h1;
            h2=.001*h2;
            dx=dx*1e-6;
            f = 0.9;
            L = .001*L;
            n=L/dx;   
            rho1 = T{Mat1,'Rho'};
            rho2 = T{Mat2,'Rho'};
            k1 = T{Mat1,'k'};
            cp1=T{Mat1,'cp'};
            k2= T{Mat2,'k'};
            cp2=T{Mat2,'cp'};
            alpha1=T{Mat1,'alpha'};
            alpha2=T{Mat2,'alpha'};
            Tm1=T{Mat1,'Tm'};
            Tm2=T{Mat2,'Tm'};
            A = ((rho1*cp1*k1)^.5+(rho2*cp2*k2)^.5)*sqrt(pi);
            Q = .5*f*rho1*h1*Vexp^2;
            Qn=Q/n;
            dtFS = .01*h1/app.Us1exp;
            out1 = (0:app.Us1exp*dtFS:h1).*1e6;
            t1=(0:dtFS:h1/app.Us1exp)*1e6;
            FShock = [t1' -out1';t1(end)+t1' -h1*1e6+out1'];
            dtTS = .01*h2/app.Us2exp;
            out2 = (0:app.Us2exp*dtTS:h2).*1e6;
            t2=(0:dtTS:h2/app.Us2exp)*1e6;
            TShock = [t2' out2';t2(end)+t2' h2*1e6-out2'];
            dt=.1e-6;
            ttot=max(FShock(end,1),TShock(end,1));
            t=dt:dt:ttot;
            z1=zeros(numel(t),1);
            z2=zeros(numel(t),1);
            for i=1:numel(t)
                z1(i)=real(sqrt(-4*alpha1*t(i)*log((A*sqrt(t(i))*(Tm1-25))/Qn)));
                z2(i)=real(sqrt(-4*alpha2*t(i)*log((A*sqrt(t(i))*(Tm2-25))/Qn)));
            end
            solidtimes=[t(find(z1==0,1,"first")) t(find(z2==0,1,"first"))];
            maxz=zeros(1,2);
            maxz(1,1)=max(abs(z1));
            maxz(1,2)=max(abs(z2));
            
        end
        
        function Fnums = Fcalc(app,T,Mat1,Mat2,dx,Vexp,Betaexp)
            dx = dx*1e-6;
            alpha1 = T{Mat1,'alpha'};
            alpha2 = T{Mat2,'alpha'};
            Vc = Vexp/(2*sind(Betaexp/2));
            Fnums = [alpha1/(dx*Vc) alpha2/(dx*Vc)];
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            %Initialize Variables
            app.Vi = 25:25:max(1025,2*app.Vexp);
            app.beta = 1:.25:30;
            app.Vtable=Velocities(app,app.Vi,app.beta);
            app.T = readtable('AppData.xlsx');
            app.FlyerMaterialDropDown.Items = app.T.Symbol;
            app.Mat1 = app.FlyerMaterialDropDown.Value;
            app.TargetMaterialDropDown.Items = app.T.Symbol;
            app.Mat2 = app.TargetMaterialDropDown.Value;
            app.VimsEditField.Value = 0;
            app.ImpactAngleDegreesEditField.Value = 0;
            app.T.Properties.RowNames = app.T.Symbol;
            app.FlyerThicknessmmEditField.Value = 0;
            app.TargetThicknessmmEditField.Value = 0;
            app.h1 = 0;
            app.h2 = 0;
            app.L = 0;
            app.dx = 0;
            app.ULimits=zeros(6,1);
            Props = string(["Particle Velocity (m/s)";"Shock Velocity (m/s)";"Impact Pressure (GPa)";"Shock Transit (microseconds)";"Longitudinal Transit (microseconds)";"Solidification Time (microseconds)";"Max Melt Extent (microns)";"Fourier Number"]);
            app.Flyvals=zeros(8,1);
            app.Targvals=zeros(8,1);
            app.UITable.Data=table(Props,app.Flyvals,app.Targvals);
            
            %Initialize Plot Limits
             app.UIAxes.XLim = [0 20000];
             app.UIAxes.YLim = [0 45];
             app.UIAxes2.XLim = [0 1000];
             app.UIAxes2.YLim = [0 30];
             app.UIAxes3.XLim = [-10 10];
             app.UIAxes3.YLim = [0 5];
            
        end

        % Value changed function: FlyerMaterialDropDown
        function FlyerMaterialDropDownValueChanged(app, event)
            app.Mat1 = app.FlyerMaterialDropDown.Value;
        end

        % Value changed function: WeldPairSwitch
        function WeldPairSwitchValueChanged(app, event)
            value = app.WeldPairSwitch.Value;
            
            if strcmp(value,'Dissimilar')
                app.TargetMaterialDropDown.Enable = 'on';
            else
                app.Mat2 = app.FlyerMaterialDropDown.Value;
                app.TargetMaterialDropDown.Value = app.FlyerMaterialDropDown.Value;
                app.TargetMaterialDropDown.Enable = 'off';
            end
        end

        % Value changed function: TargetMaterialDropDown
        function TargetMaterialDropDownValueChanged(app, event)
            app.Mat2 = app.TargetMaterialDropDown.Value;
        end

        % Value changed function: FlyerThicknessmmEditField
        function FlyerThicknessmmEditFieldValueChanged(app, event)
            app.h1 = app.FlyerThicknessmmEditField.Value;
        end

        % Value changed function: TargetThicknessmmEditField
        function TargetThicknessmmEditFieldValueChanged(app, event)
            app.h2 = app.TargetThicknessmmEditField.Value;
        end

        % Value changed function: VimsEditField
        function VimsEditFieldValueChanged(app, event)
            app.Vexp= app.VimsEditField.Value;
        end

        % Value changed function: ImpactAngleDegreesEditField
        function ImpactAngleDegreesEditFieldValueChanged(app, event)
            app.Betaexp = app.ImpactAngleDegreesEditField.Value;
%             if app.Betaexp>30
%                 app.UIAxes.YLim = [0 round(1.2*app.Betaexp)];
%                 app.beta=1:.25:round(1.2*app.Betaexp);
%                 app.Vtable=Velocities(app,app.Vi,app.beta);
%             else
%                 app.UIAxes.YLim = [0 30];
%                 app.beta=1:.25:30;
%                 app.Vtable=Velocities(app,app.Vi,app.beta);
%             end
        end

        % Value changed function: WeldLengthmmEditField
        function WeldLengthmmEditFieldValueChanged(app, event)
            app.L = app.WeldLengthmmEditField.Value;
            app.HotspotIntervalmicronsSlider.Limits=[1e-3,max(app.L/10,.5)*1e3];
        end

        % Value changed function: HotspotIntervalmicronsSlider
        function HotspotIntervalmicronsSliderValueChanged(app, event)
            app.dx = app.HotspotIntervalmicronsSlider.Value;
            app.dxEdit.Value=app.dx;
            %delete existing S2 and replace with recalculated curve using
            %new dx value
            if app.ULimits(4,1)==1 && app.Calc>0
                delete(app.hS2)
                S2=S2Calc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.ShockTable,app.S2N,app.L,app.dx);
                app.hS2=plot(app.UIAxes,S2(:,1),S2(:,2),"LineWidth",3,"LineStyle","--","DisplayName","Shock Upper Limit with Periodic Hot Spots");
                delete(app.hZ1)
                delete(app.hZ2)
                [z1,z2,t,solidtimes,maxz,FShock,TShock]=ztCalc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.Vexp,app.L,app.dx);
                app.hZ1=plot(app.UIAxes3,-z1*1E6,t*1E6,"LineWidth",3,"DisplayName",app.Mat1+" Melt Contour","Color",'r');
                app.hZ2=plot(app.UIAxes3,z2*1E6,t*1E6,"LineWidth",3,"DisplayName",app.Mat2+" Melt Contour","Color",'b');
                delete(app.FS)
                delete(app.TS)
                app.FS=plot(app.UIAxes3,FShock(:,2),FShock(:,1),"DisplayName","Flyer Shock","LineWidth",1);
                app.TS=plot(app.UIAxes3,TShock(:,2),TShock(:,1),"DisplayName","Target Shock","LineWidth",1);
                if isempty(maxz)||max(abs(maxz))==0
                    xlim = 10;
                else
                    xlim = max(abs(maxz))*1.2e6;
                end
                app.UIAxes3.XLim=[-xlim xlim];
                ylim=max(FShock(end,1),TShock(end,1));
                app.UIAxes3.YLim=[0 ylim];
                if isempty(solidtimes)
                    stimes = [NaN NaN];
                else
                    stimes = solidtimes.*1e6;
                end
                app.UITable.Data(6,2:3)=num2cell(stimes);
                app.UITable.Data(7,2:3)=num2cell(maxz.*1e6);
                Fnums = Fcalc(app,app.T,app.Mat1,app.Mat2,app.dx,app.Vexp,app.Betaexp);
                app.UITable.Data(8,2:3)=num2cell(Fnums);
            end
        end

        % Value changed function: A1CheckBox
        function A1CheckBoxValueChanged(app, event)
            app.ULimits(1,1) = app.A1CheckBox.Value;
            if app.ULimits(1,1)==1
                app.A1NEdit.Enable=1;
                app.A1NEdit.Visible=1;
                app.A1NEdit.Editable=1;
                app.A1NEdit.Value=app.A1N;
                app.NFitLabel.Visible=1;
            else
                app.A1NEdit.Enable=0;
                app.A1NEdit.Visible=0;
                app.A1NEdit.Editable=0;
            end       
        end

        % Value changed function: A2CheckBox
        function A2CheckBoxValueChanged(app, event)
            app.ULimits(2,1)= app.A2CheckBox.Value;
            if app.ULimits(2,1)==1
                app.A2NEdit.Enable=1;
                app.A2NEdit.Visible=1;
                app.A2NEdit.Editable=1;
                app.A2NEdit.Value=app.A2N;
                app.NFitLabel.Visible=1;
            else
                app.A2NEdit.Enable=0;
                app.A2NEdit.Visible=0;
                app.A2NEdit.Editable=0;
            end       
        end

        % Value changed function: S1CheckBox
        function S1CheckBoxValueChanged(app, event)
            app.ULimits(3,1)= app.S1CheckBox.Value;
            if app.ULimits(3,1)==1
                app.S1NEdit.Enable=1;
                app.S1NEdit.Visible=1;
                app.S1NEdit.Editable=1;
                app.S1NEdit.Value=app.S1N;
                app.NFitLabel.Visible=1;
            else
                app.S1NEdit.Enable=0;
                app.S1NEdit.Visible=0;
                app.S1NEdit.Editable=0;
            end       
        end

        % Value changed function: S2CheckBox
        function S2CheckBoxValueChanged(app, event)
            app.ULimits(4,1) = app.S2CheckBox.Value;
            app.HotspotIntervalmicronsSlider.Value=app.L/10;
            if app.ULimits(4,1)==1
                app.HotspotIntervalmicronsSlider.Enable = 1;
                app.HotspotIntervalmicronsSliderLabel.Enable=1;
                app.dxEdit.Visible=1;
                app.dxEdit.Enable=1;
                app.dxEdit.Editable=1;
                app.S2NEdit.Enable=1;
                app.S2NEdit.Visible=1;
                app.S2NEdit.Editable=1;
                app.S2NEdit.Value=app.S2N;
                app.NFitLabel.Visible=1;
            else
                app.S2NEdit.Enable=0;
                app.S2NEdit.Visible=0;
                app.S2NEdit.Editable=0;
                app.HotspotIntervalmicronsSlider.Enable = 0;
                app.HotspotIntervalmicronsSliderLabel.Enable=0;
                app.dxEdit.Editable=0;
            end
        end

        % Value changed function: E1CheckBox
        function E1CheckBoxValueChanged(app, event)
            app.ULimits(5,1)= app.E1CheckBox.Value;
            if app.ULimits(5,1)==1
                app.E1NEdit.Enable=1;
                app.E1NEdit.Visible=1;
                app.E1NEdit.Editable=1;
                app.E1NEdit.Value=app.E1N;
                app.NFitLabel.Visible=1;
            else
                app.E1NEdit.Enable=0;
                app.E1NEdit.Visible=0;
                app.E1NEdit.Editable=0;
            end       
        end

        % Value changed function: E2CheckBox
        function E2CheckBoxValueChanged(app, event)
            app.ULimits(6,1)= app.E2CheckBox.Value;
            if app.ULimits(3,1)==1
                app.E2NEdit.Enable=1;
                app.E2NEdit.Visible=1;
                app.E2NEdit.Editable=1;
                app.E2NEdit.Value=app.E2N;
                app.NFitLabel.Visible=1;
            else
                app.E2NEdit.Enable=0;
                app.E2NEdit.Visible=0;
                app.E2NEdit.Editable=0;
            end       
        end

        % Value changed function: dxEdit
        function dxEditValueChanged(app, event)
            app.dx = cast(app.dxEdit.Value,"double");
            app.HotspotIntervalmicronsSlider.Value=app.dx;
            %delete existing S2 and replace with recalculated curve 
            if app.ULimits(4,1)==1 && app.Calc>0
                delete(app.hS2)
                S2=S2Calc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.ShockTable,app.S2N,app.L,app.dx);
                app.hS2=plot(app.UIAxes,S2(:,1),S2(:,2),"LineWidth",3,"LineStyle","--","DisplayName","Shock Upper Limit with Periodic Hot Spots");
                delete(app.hZ1)
                delete(app.hZ2)
                delete(app.FS)
                delete(app.TS)
                [z1,z2,t,solidtimes,maxz,FShock,TShock]=ztCalc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.Vexp,app.L,app.dx);
                app.hZ1=plot(app.UIAxes3,-z1*1E6,t*1E6,"LineWidth",3,"DisplayName",app.Mat1+" Melt Contour","Color",'r');
                app.hZ2=plot(app.UIAxes3,z2*1E6,t*1E6,"LineWidth",3,"DisplayName",app.Mat2+" Melt Contour","Color",'b');
                app.FS=plot(app.UIAxes3,FShock(:,2),FShock(:,1),"DisplayName","Flyer Shock","LineWidth",1);
                app.TS=plot(app.UIAxes3,TShock(:,2),TShock(:,1),"DisplayName","Target Shock","LineWidth",1);
                if isempty(maxz)||max(abs(maxz))==0
                    xlim=10;
                else
                    xlim = max(abs(maxz))*1.2e6;
                end
                app.UIAxes3.XLim=[-xlim xlim];
                ylim=max(FShock(end,1),TShock(end,1));
                app.UIAxes3.YLim=[0 ylim];
                app.UITable.Data(6,2:3)=num2cell(solidtimes*1e6);
                app.UITable.Data(7,2:3)=num2cell(maxz*1e6);
                Fnums = Fcalc(app,app.T,app.Mat1,app.Mat2,app.dx,app.Vexp,app.Betaexp);
                app.UITable.Data(8,2:3)=num2cell(Fnums);
            end
            
        end

        % Button pushed function: CalculateButton
        function CalculateButtonPushed(app, event)
            %Initialize values
            app.Calc=app.Calc+1;
            cla(app.UIAxes)
            cla(app.UIAxes2)
            cla(app.UIAxes3)
            app.Mat1=app.FlyerMaterialDropDown.Value;
            app.Mat2=app.TargetMaterialDropDown.Value;
            app.h1=app.FlyerThicknessmmEditField.Value;
            app.h2=app.TargetThicknessmmEditField.Value;
            app.dx=50;
            app.Vexp=app.VimsEditField.Value;
            app.Betaexp=app.ImpactAngleDegreesEditField.Value;
            
            %Calculate Shock Parameters
            app.ShockTable=ShockCalcs(app,app.T,app.Mat1,app.Mat2,app.Vtable);
            %plot Vmin and Jetting Limit Curves
            [VminB,JFilter]=MinJet(app,app.Mat1,app.Mat2);
            plot(app.UIAxes,app.Vtable(:,3),VminB(:,1),"LineWidth",3,"LineStyle",'-',"DisplayName","Minimum Velocity Limit")
            plot(app.UIAxes,JFilter(:,1),JFilter(:,2),"LineWidth",3,"LineStyle",':',"DisplayName","Jetting Limit")
            %plot experimental data point
            Vcexp = app.Vexp/(2*sind(.5*app.Betaexp));
            plot(app.UIAxes,Vcexp,app.Betaexp,"Marker",".","MarkerSize",15,"DisplayName","Experimental Conditions")
            Vtabexp = [app.Vexp app.Betaexp Vcexp];
            ShockExp = ShockCalcs(app,app.T,app.Mat1,app.Mat2,Vtabexp);
            app.Up1exp=ShockExp(1,4);
            app.Up2exp=ShockExp(1,5);
            app.Us1exp=ShockExp(1,6);
            app.Us2exp=ShockExp(1,7);
            app.Pexp=ShockExp(1,8);
            
            %check which upper limit(s) to plot, calculate, and display
            if app.ULimits(1,1)==1
                A1=A1Calc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.A1N,app.Vtable);
                app.hA1=plot(app.UIAxes,A1(:,1),A1(:,2),"LineWidth",3,"LineStyle","--","DisplayName","Original Acoustic Limit");
            end
            
            if app.ULimits(2,1)==1
                A2=A2Calc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.A2N,app.Vtable);
                app.hA2=plot(app.UIAxes,A2(:,1),A2(:,2),"LineWidth",3,"LineStyle","-.","DisplayName","Acoustic Limit with Harmonic Mean Thermal Properties");
            end
            
            if app.ULimits(3,1)==1
                S1=S1Calc(app,app.T,app.Mat1,app.Mat2,app.h1,app.S1N,app.ShockTable);
                app.hS1=plot(app.UIAxes,S1(:,1),S1(:,2),"LineWidth",3,"LineStyle","-","DisplayName","Shock Upper Limit with Moving Hot Zone");
            end
            
            if app.ULimits(4,1)==1
                S2=S2Calc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.ShockTable,app.S2N,app.L,app.dx);
                app.hS2=plot(app.UIAxes,S2(:,1),S2(:,2),"LineWidth",3,"LineStyle","--","DisplayName","Shock Upper Limit with Periodic Hot Spots");
            end
            
            if app.ULimits(5,1)==1
                E1=E1Calc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.Vtable,app.E1N);
                app.hE1=plot(app.UIAxes,E1(:,1),E1(:,2),"LineWidth",3,"LineStyle","-","Marker",'o',"DisplayName","Original Elastic Upper Limit");
            end
            
            if app.ULimits(6,1)==1
                E2=E2Calc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.Vtable,app.E2N);
                app.hE2=plot(app.UIAxes,E2(:,1),E2(:,2),"LineWidth",3,"LineStyle","-","Marker",'.',"DisplayName","Elastic Upper Limit with Harmonic Mean Thermal Properties");
            end
            
            legend(app.UIAxes,"Location","southoutside")
            
            %Plot P-Up, Display relevant data in Table
            plot(app.UIAxes2,app.ShockTable(:,5),app.ShockTable(:,8).*1e-9,"LineWidth",.5,"LineStyle",'-',"Color",'b',"DisplayName","Target")
            plot(app.UIAxes2,app.Vexp-app.ShockTable(:,4),app.ShockTable(:,8).*1e-9,"LineWidth",.5,"LineStyle",'-',"Color",'r',"DisplayName","Flyer")
            legend(app.UIAxes2,"Location","southeastoutside")
            
            app.UITable.Data(1,2:3)={app.Up1exp, app.Up2exp};
            app.UITable.Data(2,2:3)={app.Us1exp, app.Us2exp};
            app.UITable.Data(3,2:3)={app.Pexp*1e-9, app.Pexp*1e-9};
            app.UITable.Data(4,2:3)={(.002*app.h1/app.Us1exp)*1e6, (.002*app.h2/app.Us2exp)*1e6};
            app.UITable.Data(5,2:3)={(.002*app.h1/app.T{app.Mat1,'C1'})*1e6, (.002*app.h2/app.T{app.Mat2,'C1'})*1e6};
            
            %Calculate and plot melt extent vs time., add data to table
            [z1,z2,t,solidtimes,maxz,FShock,TShock]=ztCalc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.Vexp,app.L,app.dx);
            app.hZ1=plot(app.UIAxes3,-z1*1E6,t*1E6,"LineWidth",3,"DisplayName",app.Mat1+" Melt Contour","Color",'r');
            app.hZ2=plot(app.UIAxes3,z2*1E6,t*1E6,"LineWidth",3,"DisplayName",app.Mat2+" Melt Contour","Color",'b');
            app.FS=plot(app.UIAxes3,FShock(:,2),FShock(:,1),"DisplayName","Flyer Shock","LineWidth",1);
            app.TS=plot(app.UIAxes3,TShock(:,2),TShock(:,1),"DisplayName","Target Shock","LineWidth",1);
            if isempty(maxz)||max(abs(maxz))==0
                    xlim = 10;
                else
                    xlim = max(abs(maxz))*1.2e6;
            end
            app.UIAxes3.XLim=[-xlim xlim];
            ylim=max(FShock(end,1),TShock(end,1));
            app.UIAxes3.YLim=[0 ylim];
            app.UITable.Data(6,2:3)=num2cell(solidtimes*1e6);
            app.UITable.Data(7,2:3)=num2cell(maxz*1e6);
            legend(app.UIAxes3,"Location","southwestoutside")  
            Fnums = Fcalc(app,app.T,app.Mat1,app.Mat2,app.dx,app.Vexp,app.Betaexp);
            app.UITable.Data(8,2:3)=num2cell(Fnums);
        end

        % Value changed function: A1NEdit
        function A1NEditValueChanged(app, event)
            app.A1N = app.A1NEdit.Value;
            delete(app.hA1)
            A1=A1Calc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.A1N,app.Vtable);
            app.hA1=plot(app.UIAxes,A1(:,1),A1(:,2),"LineWidth",3,"LineStyle","--","DisplayName","Original Acoustic Limit");
        end

        % Value changed function: A2NEdit
        function A2NEditValueChanged(app, event)
            app.A2N = app.A2NEdit.Value;
            delete(app.hA2)
            A2=A2Calc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.A2N,app.Vtable);
            app.hA2=plot(app.UIAxes,A2(:,1),A2(:,2),"LineWidth",3,"LineStyle","--","DisplayName","Acoustic Limit with Harmonic Mean Thermal Properties");
        end

        % Value changed function: S1NEdit
        function S1NEditValueChanged(app, event)
            app.S1N = app.S1NEdit.Value;
            delete(app.hS1)
            S1a=S1calc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.S1N,app.ShockTable);
            app.hS1=plot(app.UIAxes,S1a(:,1),S1a(:,2),"LineWidth",3,"LineStyle","--","DisplayName","Shock Limit with Moving Hot Zone");
        end

        % Value changed function: S2NEdit
        function S2NEditValueChanged(app, event)
            app.S2N = app.S2NEdit.Value;
            delete(app.hS2)
            S2=S2Calc(app,app.T,app.Mat1,app.Mat2,app.h1,app.h2,app.ShockTable,app.S2N,app.L,app.dx);
            app.hS2=plot(app.UIAxes,S2(:,1),S2(:,2),"LineWidth",3,"LineStyle","--","DisplayName","Shock Limit with Periodic Hot Spots");
        end

        % Value changed function: E2NEdit
        function E2NEditValueChanged(app, event)
            app.E2N = app.E2NEdit.Value;
        end

        % Value changed function: E1NEdit
        function E1NEditValueChanged(app, event)
            app.E1N = app.E1NEdit.Value;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 900 567];
            app.UIFigure.Name = 'MATLAB App';

            % Create A1CheckBox
            app.A1CheckBox = uicheckbox(app.UIFigure);
            app.A1CheckBox.ValueChangedFcn = createCallbackFcn(app, @A1CheckBoxValueChanged, true);
            app.A1CheckBox.Text = 'A1';
            app.A1CheckBox.Position = [302 515 37 22];

            % Create A2CheckBox
            app.A2CheckBox = uicheckbox(app.UIFigure);
            app.A2CheckBox.ValueChangedFcn = createCallbackFcn(app, @A2CheckBoxValueChanged, true);
            app.A2CheckBox.Text = 'A2';
            app.A2CheckBox.Position = [302 486 37 22];

            % Create S1CheckBox
            app.S1CheckBox = uicheckbox(app.UIFigure);
            app.S1CheckBox.ValueChangedFcn = createCallbackFcn(app, @S1CheckBoxValueChanged, true);
            app.S1CheckBox.Text = 'S1';
            app.S1CheckBox.Position = [302 457 37 22];

            % Create S2CheckBox
            app.S2CheckBox = uicheckbox(app.UIFigure);
            app.S2CheckBox.ValueChangedFcn = createCallbackFcn(app, @S2CheckBoxValueChanged, true);
            app.S2CheckBox.Text = 'S2';
            app.S2CheckBox.Position = [302 429 37 22];

            % Create E1CheckBox
            app.E1CheckBox = uicheckbox(app.UIFigure);
            app.E1CheckBox.ValueChangedFcn = createCallbackFcn(app, @E1CheckBoxValueChanged, true);
            app.E1CheckBox.Text = 'E1';
            app.E1CheckBox.Position = [302 401 37 22];

            % Create E2CheckBox
            app.E2CheckBox = uicheckbox(app.UIFigure);
            app.E2CheckBox.ValueChangedFcn = createCallbackFcn(app, @E2CheckBoxValueChanged, true);
            app.E2CheckBox.Text = 'E2';
            app.E2CheckBox.Position = [302 373 37 22];

            % Create FlyerMaterialDropDownLabel
            app.FlyerMaterialDropDownLabel = uilabel(app.UIFigure);
            app.FlyerMaterialDropDownLabel.HorizontalAlignment = 'right';
            app.FlyerMaterialDropDownLabel.Position = [11 509 78 22];
            app.FlyerMaterialDropDownLabel.Text = 'Flyer Material';

            % Create FlyerMaterialDropDown
            app.FlyerMaterialDropDown = uidropdown(app.UIFigure);
            app.FlyerMaterialDropDown.ValueChangedFcn = createCallbackFcn(app, @FlyerMaterialDropDownValueChanged, true);
            app.FlyerMaterialDropDown.Position = [104 509 100 22];

            % Create TargetMaterialDropDownLabel
            app.TargetMaterialDropDownLabel = uilabel(app.UIFigure);
            app.TargetMaterialDropDownLabel.HorizontalAlignment = 'right';
            app.TargetMaterialDropDownLabel.Position = [11 470 85 22];
            app.TargetMaterialDropDownLabel.Text = 'Target Material';

            % Create TargetMaterialDropDown
            app.TargetMaterialDropDown = uidropdown(app.UIFigure);
            app.TargetMaterialDropDown.ValueChangedFcn = createCallbackFcn(app, @TargetMaterialDropDownValueChanged, true);
            app.TargetMaterialDropDown.Position = [104 470 100 22];

            % Create FlyerThicknessmmEditFieldLabel
            app.FlyerThicknessmmEditFieldLabel = uilabel(app.UIFigure);
            app.FlyerThicknessmmEditFieldLabel.HorizontalAlignment = 'right';
            app.FlyerThicknessmmEditFieldLabel.Position = [11 431 121 22];
            app.FlyerThicknessmmEditFieldLabel.Text = 'Flyer Thickness (mm)';

            % Create FlyerThicknessmmEditField
            app.FlyerThicknessmmEditField = uieditfield(app.UIFigure, 'numeric');
            app.FlyerThicknessmmEditField.ValueChangedFcn = createCallbackFcn(app, @FlyerThicknessmmEditFieldValueChanged, true);
            app.FlyerThicknessmmEditField.HorizontalAlignment = 'left';
            app.FlyerThicknessmmEditField.Position = [155 431 35 22];

            % Create TargetThicknessmmEditFieldLabel
            app.TargetThicknessmmEditFieldLabel = uilabel(app.UIFigure);
            app.TargetThicknessmmEditFieldLabel.HorizontalAlignment = 'right';
            app.TargetThicknessmmEditFieldLabel.Position = [11 393 127 22];
            app.TargetThicknessmmEditFieldLabel.Text = 'Target Thickness (mm)';

            % Create TargetThicknessmmEditField
            app.TargetThicknessmmEditField = uieditfield(app.UIFigure, 'numeric');
            app.TargetThicknessmmEditField.ValueChangedFcn = createCallbackFcn(app, @TargetThicknessmmEditFieldValueChanged, true);
            app.TargetThicknessmmEditField.HorizontalAlignment = 'left';
            app.TargetThicknessmmEditField.Position = [155 393 35 22];

            % Create WeldPairSwitch
            app.WeldPairSwitch = uiswitch(app.UIFigure, 'rocker');
            app.WeldPairSwitch.Items = {'Similar', 'Dissimilar'};
            app.WeldPairSwitch.ValueChangedFcn = createCallbackFcn(app, @WeldPairSwitchValueChanged, true);
            app.WeldPairSwitch.Position = [219 396 20 45];
            app.WeldPairSwitch.Value = 'Dissimilar';

            % Create VimsEditFieldLabel
            app.VimsEditFieldLabel = uilabel(app.UIFigure);
            app.VimsEditFieldLabel.HorizontalAlignment = 'right';
            app.VimsEditFieldLabel.Position = [11 355 46 22];
            app.VimsEditFieldLabel.Text = 'Vi (m/s)';

            % Create VimsEditField
            app.VimsEditField = uieditfield(app.UIFigure, 'numeric');
            app.VimsEditField.ValueChangedFcn = createCallbackFcn(app, @VimsEditFieldValueChanged, true);
            app.VimsEditField.HorizontalAlignment = 'left';
            app.VimsEditField.Position = [131 355 45 22];

            % Create ImpactAngleDegreesEditFieldLabel
            app.ImpactAngleDegreesEditFieldLabel = uilabel(app.UIFigure);
            app.ImpactAngleDegreesEditFieldLabel.HorizontalAlignment = 'right';
            app.ImpactAngleDegreesEditFieldLabel.Position = [11 317 132 22];
            app.ImpactAngleDegreesEditFieldLabel.Text = 'Impact Angle (Degrees)';

            % Create ImpactAngleDegreesEditField
            app.ImpactAngleDegreesEditField = uieditfield(app.UIFigure, 'numeric');
            app.ImpactAngleDegreesEditField.ValueChangedFcn = createCallbackFcn(app, @ImpactAngleDegreesEditFieldValueChanged, true);
            app.ImpactAngleDegreesEditField.HorizontalAlignment = 'left';
            app.ImpactAngleDegreesEditField.Position = [155 317 35 22];

            % Create HotspotIntervalmicronsSliderLabel
            app.HotspotIntervalmicronsSliderLabel = uilabel(app.UIFigure);
            app.HotspotIntervalmicronsSliderLabel.HorizontalAlignment = 'center';
            app.HotspotIntervalmicronsSliderLabel.VerticalAlignment = 'top';
            app.HotspotIntervalmicronsSliderLabel.WordWrap = 'on';
            app.HotspotIntervalmicronsSliderLabel.Enable = 'off';
            app.HotspotIntervalmicronsSliderLabel.Position = [353 280 79 43];
            app.HotspotIntervalmicronsSliderLabel.Text = 'Hotspot Interval (microns)';

            % Create HotspotIntervalmicronsSlider
            app.HotspotIntervalmicronsSlider = uislider(app.UIFigure);
            app.HotspotIntervalmicronsSlider.ValueChangedFcn = createCallbackFcn(app, @HotspotIntervalmicronsSliderValueChanged, true);
            app.HotspotIntervalmicronsSlider.Enable = 'off';
            app.HotspotIntervalmicronsSlider.Position = [353 365 132 3];

            % Create WeldLengthmmEditFieldLabel
            app.WeldLengthmmEditFieldLabel = uilabel(app.UIFigure);
            app.WeldLengthmmEditFieldLabel.HorizontalAlignment = 'right';
            app.WeldLengthmmEditFieldLabel.Position = [11 279 104 22];
            app.WeldLengthmmEditFieldLabel.Text = 'Weld Length (mm)';

            % Create WeldLengthmmEditField
            app.WeldLengthmmEditField = uieditfield(app.UIFigure, 'numeric');
            app.WeldLengthmmEditField.ValueChangedFcn = createCallbackFcn(app, @WeldLengthmmEditFieldValueChanged, true);
            app.WeldLengthmmEditField.HorizontalAlignment = 'left';
            app.WeldLengthmmEditField.Position = [123 279 35 22];

            % Create UITable
            app.UITable = uitable(app.UIFigure);
            app.UITable.ColumnName = {'Property'; 'Flyer'; 'Target'};
            app.UITable.RowName = {};
            app.UITable.ColumnSortable = true;
            app.UITable.Position = [299 88 256 166];

            % Create CalculateButton
            app.CalculateButton = uibutton(app.UIFigure, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.CalculateButton.BackgroundColor = [0.6353 0.0784 0.1843];
            app.CalculateButton.FontSize = 16;
            app.CalculateButton.FontWeight = 'bold';
            app.CalculateButton.FontColor = [1 1 1];
            app.CalculateButton.Position = [791 466 100 92];
            app.CalculateButton.Text = 'Calculate!';

            % Create dxEdit
            app.dxEdit = uieditfield(app.UIFigure, 'numeric');
            app.dxEdit.RoundFractionalValues = 'on';
            app.dxEdit.ValueChangedFcn = createCallbackFcn(app, @dxEditValueChanged, true);
            app.dxEdit.Editable = 'off';
            app.dxEdit.HorizontalAlignment = 'left';
            app.dxEdit.Enable = 'off';
            app.dxEdit.Position = [450 293 35 22];
            app.dxEdit.Value = 50;

            % Create NFitLabel
            app.NFitLabel = uilabel(app.UIFigure);
            app.NFitLabel.FontWeight = 'bold';
            app.NFitLabel.Visible = 'off';
            app.NFitLabel.Position = [384 532 35 22];
            app.NFitLabel.Text = 'N, Fit';

            % Create A1NEdit
            app.A1NEdit = uieditfield(app.UIFigure, 'numeric');
            app.A1NEdit.ValueChangedFcn = createCallbackFcn(app, @A1NEditValueChanged, true);
            app.A1NEdit.Editable = 'off';
            app.A1NEdit.HorizontalAlignment = 'left';
            app.A1NEdit.Enable = 'off';
            app.A1NEdit.Visible = 'off';
            app.A1NEdit.Position = [375 515 32 22];

            % Create A2NEdit
            app.A2NEdit = uieditfield(app.UIFigure, 'numeric');
            app.A2NEdit.ValueChangedFcn = createCallbackFcn(app, @A2NEditValueChanged, true);
            app.A2NEdit.Editable = 'off';
            app.A2NEdit.HorizontalAlignment = 'left';
            app.A2NEdit.Enable = 'off';
            app.A2NEdit.Visible = 'off';
            app.A2NEdit.Position = [375 486 32 22];

            % Create S1NEdit
            app.S1NEdit = uieditfield(app.UIFigure, 'numeric');
            app.S1NEdit.ValueChangedFcn = createCallbackFcn(app, @S1NEditValueChanged, true);
            app.S1NEdit.Editable = 'off';
            app.S1NEdit.HorizontalAlignment = 'left';
            app.S1NEdit.Enable = 'off';
            app.S1NEdit.Visible = 'off';
            app.S1NEdit.Position = [375 457 32 22];

            % Create S2NEdit
            app.S2NEdit = uieditfield(app.UIFigure, 'numeric');
            app.S2NEdit.ValueChangedFcn = createCallbackFcn(app, @S2NEditValueChanged, true);
            app.S2NEdit.Editable = 'off';
            app.S2NEdit.HorizontalAlignment = 'left';
            app.S2NEdit.Enable = 'off';
            app.S2NEdit.Visible = 'off';
            app.S2NEdit.Position = [375 429 32 22];

            % Create E1NEdit
            app.E1NEdit = uieditfield(app.UIFigure, 'numeric');
            app.E1NEdit.ValueChangedFcn = createCallbackFcn(app, @E1NEditValueChanged, true);
            app.E1NEdit.Editable = 'off';
            app.E1NEdit.HorizontalAlignment = 'left';
            app.E1NEdit.Enable = 'off';
            app.E1NEdit.Visible = 'off';
            app.E1NEdit.Position = [375 401 32 22];

            % Create E2NEdit
            app.E2NEdit = uieditfield(app.UIFigure, 'numeric');
            app.E2NEdit.ValueChangedFcn = createCallbackFcn(app, @E2NEditValueChanged, true);
            app.E2NEdit.Editable = 'off';
            app.E2NEdit.HorizontalAlignment = 'left';
            app.E2NEdit.Enable = 'off';
            app.E2NEdit.Visible = 'off';
            app.E2NEdit.Position = [375 373 32 22];

            % Create UpperLimitsLabel
            app.UpperLimitsLabel = uilabel(app.UIFigure);
            app.UpperLimitsLabel.FontWeight = 'bold';
            app.UpperLimitsLabel.Position = [281 536 79 22];
            app.UpperLimitsLabel.Text = 'Upper Limits';

            % Create WeldingParametersLabel
            app.WeldingParametersLabel = uilabel(app.UIFigure);
            app.WeldingParametersLabel.FontWeight = 'bold';
            app.WeldingParametersLabel.Position = [18 536 120 22];
            app.WeldingParametersLabel.Text = 'Welding Parameters';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Welding Window')
            xlabel(app.UIAxes, 'Vc (m/s)')
            ylabel(app.UIAxes, 'Impact Angle \beta (Degrees)')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.PlotBoxAspectRatio = [1.79918032786885 1 1];
            app.UIAxes.NextPlot = 'add';
            app.UIAxes.Position = [511 308 277 250];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.UIFigure);
            title(app.UIAxes2, '1-D Approximated P-Up')
            xlabel(app.UIAxes2, 'Up (m/s)')
            ylabel(app.UIAxes2, 'P (GPa)')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.PlotBoxAspectRatio = [1.12981744421907 1 1];
            app.UIAxes2.NextPlot = 'add';
            app.UIAxes2.Position = [1 28 271 241];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.UIFigure);
            title(app.UIAxes3, '      Solidification and Stress Arrival')
            xlabel(app.UIAxes3, 'Distance from weld line z (\mum)')
            ylabel(app.UIAxes3, 'time t (\mus)')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.PlotBoxAspectRatio = [1.18126272912424 1 1];
            app.UIAxes3.NextPlot = 'add';
            app.UIAxes3.Position = [581 27 287 241];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Barnett_Weld_Window_v1_2_7_23_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end