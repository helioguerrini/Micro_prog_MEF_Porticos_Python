# -*- coding: utf-8 -*-
"""
Criado em Jan 14 quinta-feira 22:26:00 2016
Programa para análise de Pórticos Planos
@autor: Msc Hélio Guerrini Filho
"""
from numpy import array,zeros,pi
import function_Elementos
import function_Givens
############################## Entrada de Dados ###############################
# Matriz de conectividade
# formato: matrizconect = array([[Elemento1,nó_inicial,nó_final],[Elemento2,
# nó_inicial,nó_final],...,[ElementoN,nó_inicial,nó_final]])
matrizconect = array([[1,1,2],[2,2,3],[3,3,4],[4,4,5]])

# Matriz da malha
# malha = array([[no1,x1,y1],...,[noN,xN,yN])
malha = array([[1,0.0,0.0],[2,.125,0.0],[3,.250,0.0],[4,.375,0.0],[5,.5,0.0]])

# Módulo de elasticidade do material das barras
E = 95000000000.0

# Coeficiente de poison
ni = 0.3

#Diametros dos nós
d1,d2,d3,d4 = 0.02, 0.02, 0.02, 0.02

# Geometria e área da seção
matrizdearea = array([[(pi*d1**2)/4] , [(pi*d2**2)/4] , [(pi*d3**2)/4] , [(pi*d4**2)/4]])
momentodeinercia = array([[(pi*d1**4)/64] , [(pi*d2**4)/64] , [(pi*d3**4)/64] , [(pi*d4**4)/64]])

# Condições de contorno
#CondCont = array([graus de liberdade fixos])
CondCont = array([1,2,3])

# Carregamento
# Forca = array([[(+-)intensidade,gdl global],)
#Forca = array([[10000.,4.],[5000.,9.]])
#################### Montagem da matriz global de rigidez #####################
#Declaração da matriz de zeros
Nglobal = 3*len(malha)
K = zeros([Nglobal,Nglobal],float)
for ps in range(len(matrizconect)):
    #Obtém elementos e nós correntes    
    #el = matrizconect[ps][0] #Elemento corrente
    no1, no2 = matrizconect[ps][1], matrizconect[ps][2] #Nós do elemento
    #Obtém coordenadas dos nós
    x1, y1 = malha[no1-1][1], malha[no1-1][2]
    x2, y2 = malha[no2-1][1], malha[no2-1][2]
    #Obtém as propriedades de seção
    Area = matrizdearea[ps][0]
    Inercia = momentodeinercia[ps][0]
    #Obtém a matriz de rigidez do elemento corrente
    k = function_Elementos.ElemPort2D(E,Area,Inercia,x1,y1,x2,y2)
    #print k
    ## Montagem da matriz global
    for i in range(6):
        if i<=2:
            ig = (3*no1 - 3) + i
        else:
            ig = (3*no2 - 6) + i
        for j in range(6):
            if j<=2:
                jg = (3*no1 - 3) + j
            else:
                jg = (3*no2 - 6) + j
            #print i, ig, j, jg
            K[ig][jg] = K[ig][jg] + k[i][j]
#print K
############# Determinação da matriz global de rigidez parcial ################
##################### segundo as condições de contorno ########################
Nparcial = Nglobal - len(CondCont)
Kparcial = zeros([Nparcial,Nparcial],float) #matriz parcial
Gdllivres = zeros([Nparcial,1],float) #matriz parcial
iii = 0 # Contador auxiliar
for i in range(Nglobal):
    VarLogic = 0
    for j in range(len(CondCont)):
        CC = CondCont[j]
        if i+1 == CC:
            VarLogic = 1
    if VarLogic == 0:
        Gdllivres[iii] = Gdllivres[iii] + (i+1)
        iii = iii + 1

#print Gdllivres
            
for i in range(Nparcial):
    Gdli = int(Gdllivres[i])
    for j in range(Nparcial):
        Gdlj = int(Gdllivres[j])
        Kparcial[i][j] = Kparcial[i][j] + K[(Gdli-1)][(Gdlj-1)]
#print Kparcial
        
#################### Determinação dos deslocamentos nodais ####################
VetForca = array([[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[100.],[0.]])
deslocnod = function_Givens.sl(Kparcial,VetForca)
print deslocnod

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            

