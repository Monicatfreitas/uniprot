import csv
import time
import requests


def carregar_genes_validos():
    """Lê os resultados do UniProt e recupera apenas os genes codificantes válidos (428 genes)."""
    genes_validos = set()
    arquivos = ["resultados_uniprot.tsv", "recuperados_uniprot.tsv"]

    for arq in arquivos:
        try:
            with open(arq, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for line in reader:
                    status = line.get("Status", "")
                    gene = line.get("Gene_Original", "").strip()
                    uniprot_id = line.get("ID_UniProt", "")

                    if status in ["Encontrado", "Recuperado"] and uniprot_id not in [
                        "NA",
                        "",
                    ]:
                        if gene:
                            genes_validos.add(gene)
        except FileNotFoundError:
            pass

    return sorted(list(genes_validos))


def buscar_alphamissense_gene(gene):
    """Consulta a API do Ensembl VEP para obter predições do AlphaMissense no gene."""
    # Lookup do ID do gene no Ensembl
    url_symbol = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?content-type=application/json"
    resultados = []

    try:
        r_lookup = requests.get(url_symbol, timeout=8)
        if r_lookup.status_code == 200:
            ensg_id = r_lookup.json().get("id")

            # Consulta as variantes associadas com o plugin/anotação VEP
            url_vars = f"https://rest.ensembl.org/overlap/id/{ensg_id}?feature=variation;content-type=application/json"
            r_vars = requests.get(url_vars, timeout=8)

            if r_vars.status_code == 200:
                variacoes = r_vars.json()

                for v in variacoes:
                    rs_id = v.get("id", "")
                    consequence = v.get("consequence_type", "")

                    if rs_id.startswith("rs") and consequence == "missense_variant":
                        # Resgata atributos estendidos (incluindo pontuações do AlphaMissense se disponíveis)
                        am_score = v.get("alphamissense_score", "N/A")
                        am_class = v.get("alphamissense_class", "N/A")
                        peptide_shift = v.get("peptide_allele_string", "N/A")

                        resultados.append(
                            {
                                "rs_id": rs_id,
                                "proteina_mudanca": peptide_shift,
                                "am_score": am_score,
                                "am_class": am_class,
                            }
                        )
    except Exception as e:
        pass

    return resultados


# --- EXECUÇÃO PRINCIPAL ---
genes = carregar_genes_validos()
print(f"Iniciando busca AlphaMissense para {len(genes)} genes codificantes...\n")

with open("alphamissense_resultados.tsv", "w", encoding="utf-8") as out:
    out.write(
        "Gene\trsID_dbSNP\tProteina_Mudanca\tAlphaMissense_Score\tAlphaMissense_Classificacao\tStatus\n"
    )

    for index, gene in enumerate(genes, start=1):
        dados = buscar_alphamissense_gene(gene)

        if dados:
            print(
                f"[{index}/{len(genes)}] {gene}: {len(dados)} variantes com anotação recuperadas."
            )
            for d in dados:
                out.write(
                    f"{gene}\t{d['rs_id']}\t{d['proteina_mudanca']}\t{d['am_score']}\t{d['am_class']}\tEncontrado\n"
                )
        else:
            print(f"[{index}/{len(genes)}] {gene}: Sem anotações AlphaMissense.")
            out.write(f"{gene}\tNA\tNA\tNA\tNA\tSem_Dados\n")

        time.sleep(0.2)

print("\n" + "=" * 40)
print("PROCESSO CONCLUÍDO!")
print("Tabela gerada: 'alphamissense_resultados.tsv'")