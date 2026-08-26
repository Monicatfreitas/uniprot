import time
import requests


def buscar_ensembl(gene):
    """Busca o gene no Ensembl para obter ID e Biotipo (ex: lncRNA, miRNA)."""
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?content-type=application/json"
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            return data.get("id", "NA"), data.get("biotype", "NA")
    except Exception:
        pass
    return "NA", "NA"


def buscar_rnacentral(gene):
    """Busca o gene no RNAcentral (agregador de miRBase, LNCipedia, etc)."""
    url = f"https://rnacentral.org/api/v1/rna/?q={gene}%20AND%20TAXON:%229606%22"
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            results = data.get("results", [])
            if results:
                urs = results[0].get("rnacentral_id", "NA")
                rna_type = results[0].get("rna_type", "NA")
                return urs, rna_type
    except Exception:
        pass
    return "NA", "NA"


# 1. Carrega os genes não encontrados no UniProt
try:
    with open("genes_nao_identificados.txt", "r", encoding="utf-8") as f:
        genes_faltantes = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    print(
        "Erro: O arquivo 'genes_nao_identificados.txt' não foi encontrado."
    )
    exit()

print(f"Iniciando busca para {len(genes_faltantes)} genes no Ensembl e RNAcentral...\n")

# 2. Executa a busca e salva os resultados em TSV
with open("resultados_rnas.tsv", "w", encoding="utf-8") as out:
    out.write(
        "Gene_Original\tID_Ensembl\tBiotipo_Ensembl\tID_RNAcentral\tTipo_RNAcentral\tStatus\n"
    )

    for gene in genes_faltantes:
        ensembl_id, ensembl_biotype = buscar_ensembl(gene)
        rnacentral_id, rnacentral_type = buscar_rnacentral(gene)

        # Verifica se foi encontrado em pelo menos uma das plataformas
        if ensembl_id != "NA" or rnacentral_id != "NA":
            status = "Encontrado"
            print(
                f"[ENCONTRADO] {gene} | Ensembl: {ensembl_biotype} | RNAcentral: {rnacentral_type}"
            )
        else:
            status = "Não Encontrado"
            print(f"[AUSENTE] {gene} não possui registro nessas bases.")

        out.write(
            f"{gene}\t{ensembl_id}\t{ensembl_biotype}\t{rnacentral_id}\t{rnacentral_type}\t{status}\n"
        )
        time.sleep(0.2)

print("\n" + "=" * 40)
print("PROCESSO CONCLUÍDO!")
print("Resultados salvos em 'resultados_rnas.tsv'")