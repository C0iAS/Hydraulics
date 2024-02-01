function hidraulica
clc;clear all;close all;
%%Datos entrada
Tipo=input('digite tipo de tuberia 1=Tuberia Simple , 2=Tuberis en serie :  ');
opcion=input('Digite tipo de problema 1=Comprobacion(Q), 2=Dise�o(D), 3=Potencia(H): ');

%% Tuberia Simple
if Tipo==1
    ks=input('digite el coeficiente de rugosidad: ');
    km=input('digite el coeficiente de perdidas menores: ');
    v=input('digite el coeficiente de viscosidad: ');
    l=input('digite la longitud de la tuberia: ');
    g=9.81;
    if opcion==1
        H=input('digite la cabeza de altura: ');
        d=input('digite  el diametro: ');
        semilla=0.1;
        caudal=@(Q) ecuacion(Q,d,l,H,g,v,ks,km);
        x=fzero(caudal,semilla);
        disp('el caudal es'); disp(x);
    end
    
    if opcion==2
        
        Q=input('digite el caudal: ');
        H=input('digite la altura de cabeza: ');
        hf1=H;
        ds=input('digite un diametro comercial peque�o en pulgadas: ');
        tol=0.001;
        centinela=true;
        while centinela
            cumple=0;
            Cr1=hf1;
            d=ds*0.0254;
            V=-2*(sqrt(2*g*d*hf1)/sqrt(l))*log10((ks/(3.7*d))+(2.51*sqrt(l)*v/(d*sqrt(2*g*d*hf1))));
            Qr=V*pi*d*d/4;
            if  Qr>=Q
                cumple=1;
            else
                ds=ds+2;
            end
            hm=km*V^2/(2*g);
            if cumple==1
                hf2=H-hm;
                Cr2=hf2;
                hf1=hf2;
                errorf=abs(Cr1-Cr2);
                  
                    if errorf<=tol
                        centinela=false;
                        break
                        
                    end
             end
        end
        disp('El diametro es: ');disp(ds);
    end
    
    if opcion==3
        d=input('digite el diametro: ');
        Q=input('digite  el caudal: ');
        efici=input('digite la eficiencia en tanto por uno: ');
        z=input('digite la altura:  ');
        densidad=999.1;
        V=(4*Q/(pi*d^2));
        Re=((V*d)/v);
        ff=((0.25)/log10((ks/3.7*d)+(5.74/Re^0.9))^2);
        hf=(ff*l*V*V)/(2*d*g);
        hm=(km*V^2)/(2*g);
        H=hf+hm+z;
        pot=(1/efici)*densidad*Q*g*H;
        disp(pot);
        
    end
end
%% Tuberia en Serie
if Tipo==2
    %Datos entrada
    ntubos=input('digite el numero de tuber�as: ');
    v=input('digite la viscosidad cinematica: ');
    g=9.81;
    p=999.1;
    
    %leer datos de cada tuberia
    for i=1:ntubos
        disp('tuber�a');disp(i);
        l=input('digite la longitud: ');
        ks=input('digite la rugosidad: ');
        km=input('digite coeficiente de perdidas menores: ');
        q=input('digite el caudal lateral: ');
        Vl(i)=l;
        Vks(i)=ks;
        Vkm(i)=km;
        Vq(i)=q;
        disp(i);
    end
    %comprobacion
    if opcion==1
        %pedir los diametros de las n tuber�as y q lateral
        for i=1:ntubos
            disp('tuber�a');disp(i);
            d=input('digite el di�metro en pulgadas: ');disp(i); Vd(i)=d*0.0254;
            
        end
        centinelas=true;
        H=input('digite la cabeza total: ');
        denominador=0;
        for j=1:ntubos
            denominador=denominador+(Vl(j)/(Vd(j))^5);
        end
        
        H1=H*(Vl(1)/(Vd(1)^5))/denominador;
        while centinelas
            VH(1)=H1;
            
           % semilla=0.1;
           % caudal=@(VQ) ecuacion(VQ,Vd(1),Vl(1),H1,g,v,Vks(1),Vkm(1));
           % x=fzero(caudal,semilla);
           
           hf1=H1;
            vr=-2*(sqrt(2*g*Vd(1)*hf1)/sqrt(Vl(1)))*log10((Vks(1)/(3.7*Vd(1)))+(2.51*sqrt(Vl(1))*v/(Vd(1)*sqrt(2*g*Vd(1)*hf1))));
            VQ(1)=vr*pi*(Vd(1)^2)/4;
            for i=2:ntubos
                VQ(i)=VQ(i-1)-Vq(i-1);
            end
            
            for i=2:ntubos;
                
                V(i)=(4*VQ(i)/(pi*(Vd(i))^2));
                Re(i)=((V(i)*Vd(i))/v);
                ff(i)=((0.25)/log10((Vks(i)/3.7*Vd(i))+(5.74/(Re(i))^0.9))^2);
                hf(i)=(ff(i)*Vl(i)*V(i)*V(i))/(2*Vd(i)*g);
                hm(i)=(Vkm(i)*V(i)^2)/(2*g);
                VH(i)=hf(i)+hm(i);
                
            end
            
            sumaH=sum(VH);
            errorh=H-sumaH;
            tol=0.00001;
            if abs(errorh)<= tol
                disp('los caudales de las tuberias son: ');disp(VQ)
                centinelas=false;
                
            else
                DH=(errorh)*(Vl(1)/(Vd(1)^5))/denominador;
                H1=H1+DH;
            end
            
        end
        
    end
    
end

%diseño
if opcion==2
    
    H=input('digitar la altura de la cabeza total: ');
    tt=input('digital el valor de teta en grados: ');
    teta=(tt*pi)/180;
    
    Qt=sum(Vq);
    lt=sum(Vl);
    
    for i=1:ntubos
        Q(i)=sum(Vq(i:ntubos));
        
        hf(i)=(H*Vl(1)*cos(teta))/(lt*cos(teta));
        
        ds=input('digite un diametro comercial peque�o en pulgadas: ');
        tol=0.001;
        hf1=hf(i);
        centinela=true;
        while centinela
            Cr1=hf1;
            d=ds*0.0254;
            V=-2*(sqrt(2*g*d*hf1)/sqrt(l))*log10((ks/(3.7*d))+(2.51*sqrt(l)*v/(d*sqrt(2*g*d*hf1))));
            Qr=V*pi*d*d/4;
            if  Qr>=Q
                cumple=1;
            else
                cumple=0;
                ds=ds+2;
            end
            hm=km*V^2/(2*g);
            if cumple==1
                hf2=H-hm;
                Cr2=hf2;
                hf1=hf2;
                errorf=abs(Cr1-Cr2);
                if errorf<tol
                        centinela=false;
                 
                end
            end
        end
        
        
        Vd(i)=ds;
    end
    disp('Los diametros son: ');disp(Vd);
    
end
%potencia
if opcion==3
    z=input('digite la altura: ');
    Q=input('digite el caudal de llegada: ');
    n=input('digite la eficiencia: ');
    
    for i=1:ntubos
        disp('tuberia');disp(i);
        d=input('digite el diametro: ');Vd(i)=d*0.0254;
        
    end
    for i=1:ntubos;
        VQ(i)=Q+sum(Vq(ntubos:-1:i));
        V(i)=(4*VQ(i)/(pi*(Vd(i))^2));
        Re(i)=(V(i)*Vd(i))/(v);
        ff(i)=((0.25)/log10((Vks(i)/3.7*Vd(i))+(5.74/(Re(i))^0.9))^2);
        hf(i)=(ff(i)*Vl(i)*V(i)*V(i))/(2*Vd(i)*g);
        hm(i)=(Vkm(i)*V(i)^2)/(2*g);
        VH(i)=hf(i)+hm(i);
    end
    sumap=sum(hf)+sum(hm)+z;
    disp('la cabeza total: ');disp(sumap);
    pot=(sumap+z)*p*(sum(Vq)+Q)*g/n;
    
    disp('la potencia es: ');disp(pot);
end
end
function f=ecuacion(Q,d,l,H,g,v,ks,km)
f=((-2*sqrt((2*g*d)*(H-((km/g)*((8*Q^2)/(pi^2*d^4))))))/sqrt(l))*(log10((ks/(3.7*d))+((2.51*sqrt(l)*v)/(sqrt(2*g*d*(H-(km/g)*(8*Q^2/pi^2*d^4))))*d)))-(4*Q/(pi*d*d));
end