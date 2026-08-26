import csv
import time
import requests


def carregar_genes_unicos():
    """Lê os dois arquivos de resultado do UniProt, filtra os 'Encontrados' e remove duplicadas."""
    genes_encontrados = set()

    arquivos = ["resultados_uniprot.tsv", "recuperados_uniprot.tsv"]

    for arq in arquivos:
        try:
            with open(arq, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for line in reader:
                    # Verifica se o gene foi encontrado/recuperado e não é NA
                    status = line.get("Status", "")
                    gene = line.get("Gene_Original", "").strip()

                    if status in [
                        "Encontrado",
                        "Recuperado",
                    ] and line.get("ID_UniProt") not in ["NA", ""]:
                        if gene:
                            genes_encontrados.add(gene)
        except FileNotFoundError:
            print(f"[AVISO] Arquivo '{arq}' não encontrado. Pulando...")

    return list(genes_encontrados)


def buscar_snps_missense_ensembl(gene):
    """Consulta o Ensembl para resgatar os rsIDs do dbSNP que causam mutação missense."""
    url = f"https://rest.ensembl.org/overlap/id/{gene}?feature=variation;content-type=application/json"

    # Tenta via símbolo de gene humano primeiro
    url_symbol = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?content-type=application/json"

    snps_missense = []

    try:
        # 1. Pega o ID do Ensembl (ENSG) a partir do simbolo
        r_lookup = requests.get(url_symbol, timeout=10)
        if r_lookup.status_code == 200:
            ensg_id = r_lookup.json().get("id")

            # 2. Busca todas as variações ligadas ao gene
            url_vars = f"https://rest.ensembl.org/overlap/id/{ensg_id}?feature=variation;content-type=application/json"
            r_vars = requests.get(url_vars, timeout=10)

            if r_vars.status_code == 200:
                variacoes = r_vars.json()

                for v in variacoes:
                    # Filtra apenas variantes do dbSNP (rsID) e que sejam missense
                    rs_id = v.get("id", "")
                    consequence = v.get("consequence_type", "")

                    if rs_id.startswith("rs") and consequence == "missense_variant":
                        snps_missense.append(
                            {
                                "rs_id": rs_id,
                                "consequence": consequence,
                                "clinical_significance": ", ".join(
                                    v.get("clinical_significance", [])
                                ),
                            }
                        )

    except Exception as e:
        print(f"[ERRO] Falha ao processar {gene}: {e}")

    return snps_missense


# --- EXECUÇÃO PRINCIPAL ---
genes_unicos = carregar_genes_unicos()
print(
    f"Total de genes únicos unificados para busca: {len(genes_unicos)}"
)

with open(
    "dbsnp_missense_resultados.tsv", "w", encoding="utf-8"
) as out:
    out.write("Gene\trsID_dbSNP\tConsequencia\tSignificancia_Clinica\n")

    for index, gene in enumerate(genes_unicos, start=1):
        print(
            f"[{index}/{len(genes_unicos)}] Buscando SNPs missense para {gene}..."
        )
        snps = buscar_snps_missense_ensembl(gene)

        if snps:
            print(f" -> {len(snps)} variantes missense encontradas!")
            for snp in snps:
                out.write(
                    f"{gene}\t{snp['rs_id']}\t{snp['consequence']}\t{snp['clinical_significance']}\n"
                )
        else:
            print(f" -> Nenhuma variante missense (rsID) encontrada.")
            out.write(f"{gene}\tNA\tNA\tNA\n")

        time.sleep(0.3)

print("\n" + "=" * 40)
print("BUSCA DBSNP CONCLUÍDA!")
print("Tabela gerada com sucesso: 'dbsnp_missense_resultados.tsv'")