#--------Función 1
import csv
import os
from Bio import Entrez
from Bio import SeqIO

def source(especie, tamaño):
    Entrez.email = "elenajpsantoso@gmail.com"
    handle = Entrez.esearch(db="nucleotide", term=especie, retmax=tamaño)
    record = Entrez.read(handle)
    ids = record["IdList"]

    especies_contadas = {}

    for id in ids:
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        organism = record.annotations["organism"]
        source = record.annotations.get("source", "N/A")

        if (organism, source) in especies_contadas:
            especies_contadas[(organism, source)] += 1
        else:
            especies_contadas[(organism, source)] = 1

    # Verificar si el directorio "results" existe
    if not os.path.exists("results"):
        os.makedirs("results")

    # Guardar los resultados en el archivo CSV existente o crearlo si no existe
    archivo_csv = "results/source.csv"
    existe_archivo = os.path.exists(archivo_csv)
    modo_apertura = "a" if existe_archivo else "w"

    with open(archivo_csv, modo_apertura, newline="") as csvfile:
        writer = csv.writer(csvfile)
        if not existe_archivo:
            writer.writerow(["Organismo Fuente", "Especie", "Frecuencia"])

        for (organismo, fuente), frecuencia in especies_contadas.items():
            writer.writerow([fuente, organismo, frecuencia])

    return especies_contadas

#--------Función 2

import requests
from Bio import Entrez, SeqIO
from Bio.SeqUtils import ProtParam
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def sequences(accession_number):
    Entrez.email = "elenajpsantos@gmail.com"

    # Obtener la secuencia de ADN a partir del número de acceso
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    dna_sequence = str(record.seq)
    handle.close()

    # Traducción de la secuencia de ADN
    protein_sequence = record.seq.translate(to_stop=True)

    # Separación de péptidos iniciando en metionina
    peptides = []
    start = 0
    while True:
        try:
            peptide = protein_sequence[start:].split("M", 1)[1]
            peptides.append(peptide)
            start += len(protein_sequence[start:]) - len(peptide) + 1
        except IndexError:
            break

    # Calcular el peso molecular e índice de inestabilidad de los péptidos
    data = []
    for peptide in peptides:
        params = ProtParam.ProteinAnalysis(str(peptide))
        molecular_weight = params.molecular_weight()
        instability_index = params.instability_index()
        data.append((peptide, molecular_weight, instability_index))

    # Crear un DataFrame con los resultados
    df = pd.DataFrame(data, columns=['Peptide', 'Molecular Weight', 'Instability Index'])

    # Guardar los resultados en un archivo CSV
    df.to_csv("results/glupeptides.csv", index=False)

    # Generar el jointplot
    sns.set(style="ticks")
    g = sns.jointplot(data=df, x='Molecular Weight', y='Instability Index', kind='scatter', color='purple', s=50)
    g.fig.suptitle('Péptidos: Peso Molecular vs. Índice de Inestabilidad', fontsize=14)
    g.set_axis_labels('Peso Molecular', 'Índice de Inestabilidad')
    plt.tight_layout()
    plt.savefig("results/glupeptides.png")
    plt.close()

    return df