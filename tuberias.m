clc;clear all;close all;
Tipo=input('Digite el tipo de tuberia 1=Tuberia Simple , 2=Tuberias en serie :  ');
opcion=input('Digite tipo de problema 1=Comprobacion(Q), 2=Diseño(D), 3=Potencia(H): ');
%Tuberia Simple
if Tipo==1
    ks=input('digite la rugosidad en metros: ');
    km=input('digite el coeficiente de perdidas menores en metros: ');
    vs=input('digite la viscosidad cinematica en m^2/seg: ');
    l=input('digite la longitud de la tuberia en metros: ');
    g=9.81;
    %Comprobacion tuberia simple
    if opcion==1
        H=input('digite la cabeza de altura en metros: ');
        d=input('digite  el diametro en metros: ');
        A=((d^2)/4)*pi;
        hf1=H;
        hf=0;
        cont=1;
        e=1;
        while e>0.0001
            if(cont==1)
                hf=H;
            end
        hf1=hf;    
        V=(-2*sqrt(2*g*d*hf1))/(sqrt(l))*log10(((ks)/(d*3.7))+((2.51*vs*sqrt(l))/(d*sqrt(2*g*d*hf1))));
        hf=H-km*(V^2)/(2*g);
        cont=cont+1;
        e=abs(hf-hf1);
        end
        Q=V*A;
        hm=H-hf;
        disp('las perdidas por friccion en metros son');disp(hf);
        disp('las perdidad menores por accesorios en metros son');disp(hm);
        disp('el caudal es en m^3/seg');disp(Q);
        disp('numero de iteraciones');disp(cont);
    end
    %Diseño tuberia simple
    if opcion==2
        Qd=input('digite el caudal de diseño en m^3/seg: ');
        H=input('digite la altura de cabeza en metros: ');
        hf1=H;
        Q=0;
        d=0.0508;
        g=9.81;
        while Q<Qd
            V=(-2*sqrt(2*g*d*hf1))/(sqrt(l))*log10(((ks)/(d*3.7))+((2.51*vs*sqrt(l))/(d*sqrt(2*g*d*hf1))));
            Q=V*((d^2)/4)*pi;
            if Q<Qd
                d=d+0.0508;
            end
        end
        e=1;
        hf=0;
        cont=1;
        while e>0.0001
            if(cont==1)
                hf=H;
            end
        hf1=hf;    
        V=(-2*sqrt(2*g*d*hf1))/(sqrt(l))*log10(((ks)/(d*3.7))+((2.51*vs*sqrt(l))/(d*sqrt(2*g*d*hf1))));
        hf=H-km*(V^2)/(2*g);
        cont=0;
        e=abs(hf-hf1);
        end
        Q=V*((d^2)/4)*pi;
        hm=H-hf;
        disp('las perdidas por friccion en metros son');disp(hf);
        disp('las perdidad menores por accesorios en metros son');disp(hm);
        disp('el caudal mayor al de diseño obtenido en m^3/seg es');disp(Q);
        disp('el diametro utilizado es');disp(d);
    end
    %Calculo de potencia tuberia simple
    if opcion==3
        d=input('digite el diametro de la tuberia en metros: ');
        Q=input('digite  el caudal requerido en m^3/seg: ');
        n=input('digite la eficiencia de la bomba en porcentaje: ');
        z=input('digite la diferencia de altura a vencer en metros:  ');
        densidad=999.1;
        V=Q*4/((d^2)*pi);
        Re=V*d/vs;
        e=1;
        f1=0.02;
        while e>0.00001
            f=1/(-2*log10((ks/(3.7*d))+2.51/(Re*sqrt(f1))))^2;
            e=abs(f-f1);
            f1=f;
        end
        hm=km*V^2/(2*g);
        hf=f*l*V^2/(d*2*g);
        H=hf+hm+z;
        pot=(n*densidad*Q*H*g/100)/1000;
        disp('la potencia requerida en kW es');disp(pot);
    end
end
if Tipo==2
    %tuberias en serie
    ntubos=input('digite el numero de tuberías en serie: ');
    vs=input('digite la viscosidad cinematica en m^2/seg: ');
    g=9.81;
    p=999.1;
    for i=1:ntubos
        disp('TUBERIA');disp(i);
        l(i)=input(' digite la longitud en metros: ');
        ks(i)=input('digite la rugosidad en metros: ');
        km(i)=input('digite coeficiente de perdidas menores en metros: ');
        q(i)=input('digite el caudal lateral en m^3/seg: ');
    end
    %Comprobacion tuberias en serie
    if opcion==1
        H=input('Digite la cabeza total en metros: ');
        for i=1:ntubos
            fprintf('Digite diametro de tuberia %i en metros',i)
            d(i)=input(' : ');
        end
        den=0;
        num=(l(1)/d(1)^5);
        for i=1:ntubos
            den=den+(l(i)/d(i)^5);
        end
        hf(1)=(num/den)*H;
        V(1)=(-2*sqrt(2*g*d(1)*hf(1)))/(sqrt(l(1)))*log10(((ks(1))/(d(1)*3.7))+((2.51*vs*sqrt(l(1)))/(d(1)*sqrt(2*g*d(1)*hf(1)))));
        for i=1:ntubos
        A(i)=((d(i)^2)/4)*pi;
        end
        Q(1)=V(1)*A(1);
        hm(1)=km(1)*(V(1)^2)/(2*g);
        for i=2:ntubos
            Q(i)=Q(i-1)-q(i-1);
            V(i)=Q(i)/A(i);
            Re(i)=V(i)*d(i)/vs;
            f(i)=0.25/(log10(ks(i)/(d(i)*3.7)+(5.74/(Re(i)^0.9))))^2;
            hf(i)=f(i)*l(i)*V(i)^2/(2*d(i)*g);
            hm(i)=km(i)*(V(i)^2)/(2*g);
        end
        Ht=0;
        for i=1:ntubos
            Ht=Ht+hf(i)+hm(i);
        end
        deltah=num/den*(H-Ht);
        while deltah>0.00001
        hf(1)=deltah+hf(1);
        V(1)=(-2*sqrt(2*g*d(1)*hf(1)))/(sqrt(l(1)))*log10(((ks(1))/(d(1)*3.7))+((2.51*vs*sqrt(l(1)))/(d(1)*sqrt(2*g*d(1)*hf(1)))));
        Q(1)=V(1)*A(1);
        hm(1)=km(1)*(V(1)^2)/(2*g);
        for i=2:ntubos
            Q(i)=Q(i-1)-q(i-1);
            V(i)=Q(i)/A(i);
            Re(i)=V(i)*d(i)/vs;
            f(i)=0.25/(log10(ks(i)/(d(i)*3.7)+(5.74/(Re(i)^0.9))))^2;
            hf(i)=f(i)*l(i)*V(i)^2/(2*d(i)*g);
            hm(i)=km(i)*(V(i)^2)/(2*g);
        end
        Ht=0;
        for i=1:ntubos
            Ht=Ht+hf(i)+hm(i);
        end
        deltah=(num/den)*(H-Ht);
        end
        for i=1:ntubos
            fprintf('el caudal en tuberia %i en m^3/seg es',i)
            disp(Q(i))
            fprintf('las perdidas menores en tuberia %i en metros son',i)
            disp(hm(i))
            fprintf('las perdidas por friccion en tuberia %i en metros son',i)
            disp(hf(i))
        end
    end
    if opcion==2
         Qr=input('digite el caudal requerido: ');
         H=input('digitar la altura de la cabeza total en metros: ');
         for i=1:ntubos
         fprintf('digital el valor de teta %i en grados',i);
         tt(i)=input(' : ');
         teta(i)=(tt(i)*pi)/180;
         end
         lt=sum(l);
         for i=1:ntubos
             hf(i)=(H*l(i)*cos(teta(i)))/(lt*cos(teta(i)));
         end
         Qd(1)=Qr;
         for i=2:ntubos
             Qd(i)=Qr-q(i-1);
         end
         for i=1:ntubos
             Q(i)=0;
             d(i)=0.0508;
         end
         for i=1:ntubos
         while Q(i)<Qd(i)
              V(i)=(-2*sqrt(2*g*d(i)*hf(i)))/(sqrt(l(i)))*log10(((ks(i))/(d(i)*3.7))+((2.51*vs*sqrt(l(i)))/(d(i)*sqrt(2*g*d(i)*hf(i)))));
              Q(i)=V(i)*(d(i)^2)*pi/4;
              if Q(i)<Qd(i)
                  d(i)=d(i)+0.0508;
              end
         end
         end
         for i=1:ntubos
             hm(i)=km(i)*V(i)^2/(2*g);
             hf(i)=hf(i)-hm(i);
         end
         hft=sum(hf);
         disp('Las perdidas por friccion totales son: ');disp(hft)
         hmt=sum(hm);
         disp('las perdidas menores totales son: ');disp(hmt)
         disp('los diametros respectivos de cada tuberia en metros son: ');disp(d)
    end
    if opcion==3
        z=input('digite la diferencia de altura a vencer en metros: ');
        Qf=input('digite el caudal de llegada en m^3/seg: ');
        n=input('digite la eficiencia de la bomba en porcentaje: ');
        for i=1:ntubos
            fprintf('Digite diametro de tuberia %i en metros',i)
            d(i)=input(' : ');
        end
        Q(1)=Qf+q(1);
        for i=2:ntubos
            Q(i)=Q(i-1)-q(i-1);
        end
        for i=1:ntubos
            V(i)=Q(i)/((d(i)^2)*(pi/4));
            Re(i)=V(i)*d(i)/vs;
            f(i)=0.25/(log10(ks(i)/(d(i)*3.7)+(5.74/(Re(i)^0.9))))^2;
            hf(i)=f(i)*l(i)*V(i)^2/(2*d(i)*g);
            hm(i)=km(i)*(V(i)^2)/(2*g);
        end
        hft=sum(hf);
        hmt=sum(hm);
        Ht=z+hft+hmt;
        pot=(((999.1*g*Ht*Q(1)*100)/n)/1000);
        disp('las perdidas menores totales en metros son:');disp(hmt)
        disp('las perdidas totales por friccion en metros son:');disp(hft)
        disp('la potencia en kW es:');disp(pot)
    end
end
        
        


         
             
             
             
        
        
            
        
        
    
    
        

        
        

        
        
            
            
        
        
        
        
        
        
