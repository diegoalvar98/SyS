/*____________________________________________________________________________________________

Simulador de síntesis de apertura y algoritmo CLEAN.
Programa creado en la versión de scilab 6.0.1 con el paquete IPCV para el procesado de imágenes
Autores: Diego Álvarez Ortega e Iván Reyes Rodríguez                                    2018
_______________________________________________________________________________________________*/

//Limpiamos todas las variables
clear;
scf()

/*
______________________________
----- PARAMETROS -------------
______________________________
*/
 
maxxy=200         //Tamaño en km de la matriz de las antenas
minxy=-maxxy
lambda=106*10^-3  //Longitud de onda, en metros
step=1          //Paso, en km/division
angulars=20       //tamaño angular de la imagen observada, en arcosegundos
delta =%pi/2      //Declinación de la fuente (0-horizonte, pi/2 polo norte). Array fijo en latitud 90
umbral=20        //Umbral entre 0 y 255 para detener el algoritmo CLEAN
modelo="cygnus.jpg"//ruta de la imagen
cleansize=150      //Tamaño de la gaussiana con la que se genera la imagen CLEAN
N=100               //Iteraciones de clean
gain=0.2           //Ganancia para el cleaan, entre 0.1 y 0.3

/*___________________________
----  EJEMPLOS DE ARRAYS  ----
_____________________________*/
//Array complejo
//xant=[4 1.5 2.2 3 3.5 4 5.7 6 7 8 8 9 4 4.5].*10
//yant=[6 0.5   -0.5 0 6   0.8 0 -2 0 0.1 3 0 3.5 -3.6].*10

//4 radios
/*
xant=[0 4 6 3 7  10 6] //Coordenadas en km de las antenas, deben ser multiplos de step
yant=[0 0 0 2 -1 10 2]  
*/


xant=[-2 -1 0 1 2 0 0  0  0]
yant=[0   0 0 0 0 1 2 -1 -2]


n=size(xant)(2) //numero de antenas
/*___________________________
----     ALGORITMO      ----
_____________________________*/

nxant=xant./step-minxy/step+1 
nyant=yant./step-minxy/step+1

tx=minxy:step:maxxy
ty=minxy:step:maxxy

nx=tx./step-minxy/step+1
ny=ty./step-minxy/step+1


map=zeros(size(tx)(2),size(ty)(2))


for k=1:1:n
        map(nxant(k),nyant(k))=1    //Generamos una matriz cuadrada con 1 en los puntos donde hay antena
end

subplot(2,4,1)
plot(xant,yant,'xb')                //Se imprime la disposicion espacial de las antenas
title("Antenas")
xlabel(gca(),"km")
ylabel(gca(),"km")
gca().data_bounds=[min(xant)-0.5, max(xant)+0.5,min(yant)-0.5, max(yant)+0.5]

ft=fft(map)
c=fftshift(fftshift(ifft(conj(ft).*ft),1),2); //Autocorrelacion de la apertura, es la función de muestreo


//Calculamos ahora las coordenadas uv de cada punto muestreado
NFFT=size(c)(2)+1
[nxv,nyv]=find(abs(c)>0.01)     
xv=(nxv-NFFT/2)*step*1000/lambda/1000000//n-km-m-lambdas-Klambda
yv=(nyv-NFFT/2)*step*1000/lambda/1000000

//Este bucle aplica la rotación de la tierra a distintas horas y acumula los puntos muestreados
nmuestreos=uint8(size(xv)(2))

xvtotal=0
yvtotal=0
xyv=[xv ;yv]
clear('xyvnuevas')
for t=-6:0.2:6
    anguloh=t/12*%pi+%pi/2
    for i=1:1:nmuestreos
    xyvnuevas=[sin(anguloh) ,cos(anguloh) ; -sin(delta)*cos(anguloh),sin(delta)*sin(anguloh)]*xyv(:,i)
    xvtotal=[xvtotal xyvnuevas(1)]
    yvtotal=[yvtotal xyvnuevas(2)]
    end
end


//Funcion de muestreo final
xv=xvtotal  
yv=yvtotal

//Solo es necesario eliminar las baseline en 0
nmuestreos=size(xv)(2)
cuenta=0
for k=1:1:nmuestreos
        if ((xv(k-cuenta)==0.)&(yv(k-cuenta)==0.)) then
                xv(k-cuenta)=[ ]
                yv(k-cuenta)=[ ] //(si no se salta un indice.)
                cuenta=cuenta+1;
        end
end

//Creamos una matriz con la función de muestreo
nmuestreos=uint8(size(xv)(2))
c=zeros(size(tx)(2),size(ty)(2))
nxv=uint8(xv*10000000*lambda/1000/step+NFFT/2)
nyv=uint8(yv*10000000*lambda/1000/step+NFFT/2)

for k=1:1:size(nxv)(2)
       if (nxv(k)>0)&(nxv(k)<=size(tx)(2))&(nyv(k)>0)&(nyv(k)<=size(tx)(2)) then
        c(nxv(k),nyv(k))=1
        end
end
//Sgrayplot(nx,ny,c)  //prueba para ver que se esta mapeando bien la funcion de muestreo


subplot(2,4,2)
plot(xv,yv,'.r',"MarkerSize",3)
gca().axes_reverse=["on" "on" "off"]
title("S(u,v)")
xlabel(gca(),"U(kλ)")
ylabel(gca(),"V(kλ)")

psf=fftshift(fftshift(fft(c,-1),1),2); //LA función de transferencia es la fft de S(u,v)
subplot(2,4,3)

//Prueba y error, sobre todo error, pero funciona
//xdb=(nx-NFFT/2).*lambda*0.00029*100000
xdb=(nx-NFFT/2).*lambda*0.00029*100000/NFFT*15/step
ydb=(ny-NFFT/2).*lambda*0.00029*100000/NFFT*15/step
Sgrayplot(xdb,ydb,abs(psf)) 

gca().axes_reverse=["on" "off" "off"]
gcf().color_map = jetcolormap(64);
title("PSF")
xlabel(gca(),"as")
ylabel(gca(),"as")
/*__________________________________________________________________________
-------     Obtención de la respuesta a una imagen ejemplo          --------
____________________________________________________________________________*/
//Imagen redimensionada:
subplot(2,4,5)
object=imread(modelo);
//calculamos las dimensiones de la matriz basados en la resolucion de la psf
//objsize=int(angulars/lambda/0.00029/100000)
objsize=angulars*step/15*NFFT/(100000 *lambda*0.00029)
objectr=rgb2gray(imresize(object,[objsize,objsize]))
imshow(objectr)
title("Funcion Brillo I")

//Mostramos la TF de la imagen, correspondiente a la funcion visibilidad sin muestrear
subplot(2,4,6)
objectTF=fft(double(objectr))
imshow(uint8(abs(fftshift(fftshift(objectTF,2),1))/max(abs(objectTF))*255))
title("V(Sin muestreo)")

//Convolucion, mostramos psf y dirtyimage
subplot(2,4,7)
psfim=uint8(abs(psf)/max(abs(psf))*255)
dirtybeam=imcrop(psfim',[(NFFT-objsize)/2,(NFFT-objsize)/2,objsize ,objsize ])
dirtyconv=conv2(double(objectr),double(psfim),'same')//dirtybeam o psfim
dirtyimage=uint8(abs(dirtyconv)/max(abs(dirtyconv))*255)
imshow(dirtybeam)
title("PSF reescalado")

subplot(2,4,8)
title("Dirty Image")
imshow(dirtyimage)

subplot(2,4,4)
title("S·V")
sampleo=ifft(double(dirtyimage))
imshow(uint8(abs(fftshift(fftshift(sampleo,2),1))/max(abs(sampleo))*255))


/*
figura=scf()
ax=gca()
Sgrayplot(xdb,ydb,abs(psf)) 
ax.axes_reverse=["on" "off" "off"]
figura.color_map = graycolormap(64)
ax.data_bounds=[-angulars/2 angulars/2,-angulars/2 angulars/2]


xs2png(figura,'PSFraw.png');
*/

//CLEAN


scf()

h=double(psfim)./255; //psfim tarda mas pero es mejor, dirtybeam es mas rapido
//hclean=conv2(h,fspecial('gaussian', [objsize, objsize],iqr(dirtybeam)),"same");
hclean=fspecial('gaussian', [objsize, objsize],cleansize*iqr(h)); //Es una gaussiana ajustada al dirty beam
//h=[1/2 1/2 1/2; 1/2 1 1/2;1/2 1/2 1/2];
//I=rgb2gray(imread('galaxy.png'));
Im=double(dirtyimage);
corre=25.*ones(Im);
//Imclean=zeros(Im);
Mpuntos=zeros(Im);
Mfinal=zeros(Im);
z=0;
subplot(1,2,1)
title("Clean")
subplot(1,2,2)
title("Residue")
for k=1:+1:N
    Paso=zeros(Im);
    [m,p]=max(Im);
    iprima=p(1);jprima=p(2);  
    Paso(iprima,jprima)=Im(iprima,jprima)*(gain);
    //imshow(uint8(conv2(Paso,h,"same")))
    //Imclean=Imclean+conv2(Paso,h,"same");
    Im=Im-conv2(Paso,h,"same");
   //      for r=1:+1:256 Este de aqui era un algoritmo para que no pudiese sustraer por debajo de cero, pero daba unos efectos un poco raros.
   //           for s=1:+1:256
   //               if Im(r,s)<0 then Mfiltradacorre(r,s)=Mfiltradacorre(r,s)+Im(r,s);Im(r,s)=0; end end end
     Mfinal=Mfinal+Paso;  
     //imshow(Mpuntos); 
     
    
     mconv=conv2(Mfinal,hclean,"same")
     Mclean=uint8(mconv*255/max(mconv))
     subplot(1,2,1)
     imshow(Mclean);
     subplot(1,2,2)
     title("Residue (N=" + string(k) +")")
     imshow(uint8(max(Im)*(Im-ones(Im)*min(Im))/(max(Im)-min(Im))));
          
     
     if double(max(Im))<umbral then
         disp("Se ha alcanzado el umbral minimo, imagen procesada")
         break
     end
end

scf()
output=double(Mclean)+(max(Im)*(Im-ones(Im)*min(Im))/(max(Im)-min(Im)))
imshow(uint8(255*(output-ones(output)*min(output))/(max(output)-min(output))));

if k==100 then
    disp("Alcanzado numero maximo de iteraciones")
end



