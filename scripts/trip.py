#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CompararVcf.py
Programa que dado una carpeta crea una matriz con todos los valores de TP,
FP, TN, Sensibilidad, Precisión y FDR. Para cada uno de los VC empleados.

Para ello la estructura de los archivos es:
NombreMuestra/snv(Numero mutaciones)/Truth|TumorVcf

@author: pamonlan
"""

from pprint import pprint
from pdbio.vcfdataframe import VcfDataFrame
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tabulate import tabulate
import numpy as np
from multiprocessing import Process, Value, Array, Queue
import time as t
import argparse

def read_vcf(vcf_path):
    """
    Función que lee un fichero VCF y devuelve un DataFrame de la clase panda
    con la información de las variantes
    """

    vcf = VcfDataFrame(path=vcf_path)
    vcf.sort()
    vcf_df = vcf.df
    vcf_df = vcf_df[['#CHROM', 'POS', 'ID', 'REF', 'ALT']]

    return vcf_df


def diferencias(ref_df, call_df):
    """
    Función que dado dos DF de panda, los cuales contienen la información del VCF
    permite obtener el número de variantes comunes, y diferentes de cada archivo.

    """
    compartive_df = pd.merge(ref_df, call_df,
                              indicator=True,
                              how='outer')
    #print(compartive_df)
    TP = len(compartive_df[compartive_df['_merge'] == "both"])
        #diff_df.to_csv('treatment_out/variant_calling/diff.tp.csv')
    FN = len(compartive_df[compartive_df['_merge'] == "left_only"])
        #diff_df.to_csv('treatment_out/variant_calling/diff.fn.csv')
    FP = len(compartive_df[compartive_df['_merge'] == "right_only"])
        #diff_df.to_csv('treatment_out/variant_calling/diff.fp.csv') 
            
    return TP, FN, FP

def metricas(ref_df,call_df, TumorPath):
    TP, FN, FP = diferencias(ref_df, call_df,)

    try:
        Precision = round(TP / ( TP + FP ),2)
        Sensibilidad = round(TP / ( TP + FN ), 2)
        FDR = round(FP / ( TP + FP ),2)
        F1 = round(( 2 * Sensibilidad * Precision) / ( Sensibilidad + Precision ), 2)
    
    except Exception as Excpt:
        #print("# Error %s" %Excpt)
        #print(TumorPath)
        Precision, Sensibilidad, FDR, F1 = None, None, None, None
    
    return TP, FN, FP, Precision, Sensibilidad, FDR, F1


def graficas (lResultados, sFile):
    
    lCol = ["Muestra", sFile, "Replica", "Variant Caller", "TP","FN", "FP", \
            "Precision", "Sensibilidad","FDR", "F1"]
        
    df_total = pd.DataFrame(lResultados, \
                            columns = lCol)
        
    #df_total = df_total.apply(pd.to_numeric)
    
    # Representación de métricas
    
    fig, (px, sx, fx) = plt.subplots(1,3,figsize=(15,8)) #incluir hilos y procesos en 1 imagen

    ## Precisión vs AF
    sns.pointplot(x=sFile, y="Precision", hue="Variant Caller", data = df_total, ci="sd", \
    capsize=.1, errwidth=1, ax=px)

    ### Sensibilidad vs AF
    sns.pointplot(x=sFile, y="Sensibilidad", hue="Variant Caller", data = df_total, ci="sd", \
    capsize=.1, errwidth=1, ax=sx)
        
    ## F1Score vs AF
    sns.pointplot(x=sFile, y="F1", hue="Variant Caller", data = df_total, ci="sd", \
    capsize=.1, errwidth=1, ax=fx)

    plt.show()
    
    return df_total
    

def print_table(dfResultados, sFile):
    
    ### DATOS BRUTOS
    print(tabulate(dfResultados.sort_values(by=[sFile,"Variant Caller"]), \
                   headers='keys', tablefmt='latex',showindex=False))

    ### DATOS MEDIA + STD
    df_p_stats = dfResultados.groupby([sFile,"Variant Caller"])["Precision",\
                                    "Sensibilidad","FDR", "F1"].agg(["mean"])
    
    print(tabulate(df_p_stats, headers='keys',tablefmt='latex'))
    
    result = "./"+str(sFile)+"sTablaResultadosRaw.csv"
    stats = "./"+str(sFile)+"Tabla_resultados_media.csv"
    dfResultados.to_csv(result)
    df_p_stats.to_csv(stats)


def obtener_metricas(lPath,lMetadata,qResultados):
    lLista = len(lPath)
    
    for i in range (lLista):
        #Obtenemos la tupla con los Paths de los VCFS
        tPath = lPath[i]
        #Obtenemos la lista con la información del análisis
        lData = lMetadata[i]
        TruthPath = tPath[0]
        TumorPath = tPath[1]
        
        try:
        #Leemos los VCFS
            TruthVCF = read_vcf(TruthPath)
            TumorVCF = read_vcf(TumorPath)
            
        except Exception as excpt:
            print("Error leyendo el archivo: %s" % TumorPath)
            print(excpt)
        
        else:
            #Obtenemos las métricas
            TP, FN, FP, Precision, Sensibilidad, FDR, F1 = \
                        metricas(TruthVCF, TumorVCF, TumorPath)
            lTemp = [TP, FN, FP, Precision, Sensibilidad, FDR, F1]
            
            #Concatenamos las listas las listas
            lData.extend(lTemp)
            
            #Lo almacenamos en la cola de resultados
            qResultados.put(lData)
        
        #Analizado
        
        
     
def leer_cola(cola):
    """
    Función que permite obtener todos los objetos introducidos en una cola.
    
    Args:
        cola (Queue) Cola de proceosos sobre la que almacena los resultados
    
    Returns:
        lResultados (list) Lista con las posiciones de las cadenas que han pasado el filtro
    """
    
    lResultados = []
    while not cola.empty():
        lResultados.append(cola.get())
    
    return lResultados


def dividir_procesos(lPath, lMetadata, np):
    npro =  np #Se almacena el número de procesos
    
    # Se calcula el tamaño de la subcadena por procesos
    sub_len = len(lPath) // npro
    
    # Calcular las posiciones de mas que tiene que procesar el ultimo proceso
    resto =  len(lPath) % npro
    
    lista_proc = []
    qResultados = Queue()  # Cola que guarda las posiciones coincidentes
    
    # Bucle que inicia cada proceso. 
    for p in range(npro):
        
        #Se comprueba que no sea el último proceso
        if p != (npro - 1):
            # Se calcula la posición inicial y finalpara el proceso
            inicio = p * sub_len
            final = inicio + sub_len
            
            #Se añade el proceso en la lista
            lista_proc.append(Process(target=obtener_metricas,
                                      args=(lPath[inicio:final], \
                                            lMetadata[inicio:final], qResultados)))
            lista_proc[p].start()
            
        # Si el proceso es el ultimo debe llegar hasta el final de la cadena
        else:
            inicio = p * sub_len
                #Se calcula la posicion inicial y final teniendo en cuenta el resto
            final = inicio + sub_len + resto
            lista_proc.append(Process(target=obtener_metricas,
                                      args=(lPath[inicio:final], \
                                            lMetadata[inicio:final], qResultados)))
            lista_proc[p].start()

    # Se espera a que acaben todos los procesos.
    for i in range(npro):
        lista_proc[i].join()

    #Se obtiene una lista con los resultados
    lResultados = leer_cola(qResultados)
    
    return lResultados


def argumentos():
    """
    Función que recibe los parámetros introducidos como argumento a la hora
    de llamar al programa desde la terminal. Facilita su automatización
    
    Returns:
        args (argparse) Objeto que almacena los argumentos
    """
    
    parser=argparse.ArgumentParser (
            description ='''            
'''
            )
    
    parser.add_argument(
            "-s","--sample",
            dest="sample",
            action="store",
            required=True,
            help="Ruta al archivo a analizar"
            )
    
    
    parser.add_argument(
            "-t","--type",
            dest="type",
            action="store",
            required=True,
            help="Tipo de análisis a efectuar: AF, Depth, SNV"
            )
    
    parser.add_argument(
            "-p","--process",
            dest="process",
            action="store",
            required=True,
            help="Numero de procesos"
            )

  
    args = parser.parse_args()
    
    return args

    
def main():
    """
    Función que obtiene los TP, FP, FP y calcula la precisión, sensibilidad y FDR
    """    
   
    tCalling = ("mutect2","pisces","octopus","somvarius","sinvict","outlyzer","shaed")
        
    tIntersect = ("M2Out", "M2O", "M2P", "M2Sin", "M2Som", "OcOut", "SinvPi",\
                  "SomOut", "OcSinv", "sinsom", "OcSom", "OP", \
                      "POut","PSom", "SinOut")
    lPath = []
    lMetadata = []
    
    tAf_exome = (0.02, 0.05, 0.075, 0.1, 0.15, 0.3, 0.4, 0.6 )
    tAf = (0.005, 0.01, 0.02, 0.05, 0.075, 0.1, 0.15, 0.3, 0.4, 0.6 )
    tSnv = (10, 20, 40, 60, 80)
    tDepth = (10, 50, 100, 150, 200, 250, 300)
    
    argu = argumentos()
    sVal = argu.type
    sVal = sVal.upper()
    sample = argu.sample
    np = int(argu.process)
    nMax = 10
    
    if sVal == "AF":
        tList = tAf
        sFile = "af."
        
    
    elif sVal == "DEPTH":
        tList = tDepth
        sFile = "depth."
        nMax = 5
    
    else:
        tList = tSnv
        sFile = "snv."
    
    for muestra in [sample,]:
        for j in tList:
            sPath = sFile+str(j)
            
            for rep in range(1,nMax+1):
                
                if sVal == "DEPTH":
                    NombSnvRep = "%s.snv.80.%s.%i" %(muestra,str(j),rep)
                    PathTruthVCF = "%s/%s/TruthVcf/%s.truth.vcf" \
                    %(muestra,sPath,NombSnvRep)
                
                else:
                    NombSnvRep = "%s.%s.%i" %(muestra,sPath,rep)
                    PathTruthVCF = "%s/%s/TruthVcf/%s.truth.vcf" \
                    %(muestra,sPath,NombSnvRep)                
                
                for vc in tCalling:
                    if vc == "octopus" or vc == "pisces":
                        PathTumorVCF = "%s/%s/TumorVcf/%s/%s.%s.pass.vcf" \
                        %(muestra,sPath,NombSnvRep,NombSnvRep,vc)
                        
                    elif vc == "shaed":
                        PathTumorVCF = "%s/%s/TumorVcf/%s/%s.%s.somatic.vcf" \
                        %(muestra,sPath,NombSnvRep,NombSnvRep,vc)
                    else:
                        PathTumorVCF = "%s/%s/TumorVcf/%s/%s.%s.normal.snps.vcf" \
                        %(muestra,sPath,NombSnvRep,NombSnvRep,vc)
                        
                    tPath = (PathTruthVCF, PathTumorVCF)
                    lData = [muestra, j, rep, vc]
                        
                    lPath.append(tPath)
                    lMetadata.append(lData)
                        
    
                    
   
    lResultados = dividir_procesos(lPath, lMetadata, np)
    
    dfResultados = graficas (lResultados, sFile)
    
    print_table(dfResultados, sFile)

if __name__ == "__main__":
    main()
