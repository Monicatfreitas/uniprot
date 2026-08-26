import time
import requests

# 1. Carrega apenas os genes que falharam na primeira tentativa
try:
    with open("genes_nao_identificados.txt", "r", encoding="utf-8") as f:
        genes_faltantes = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    print(
        "Erro: O arquivo 'genes_nao_identificados.txt' não foi encontrado."
    )
    exit()

genes_ainda_ausentes = []
url = "https://rest.uniprot.org/uniprotkb/search"

print(f"Repesquisando {len(genes_faltantes)} genes com busca ampliada...\n")

with open("recuperados_uniprot.tsv", "w", encoding="utf-8") as out:
    out.write(
        "Gene_Original\tID_UniProt\tNome_Proteina\tOrganismo\tStatus\n"
    )

    for gene in genes_faltantes:
        try:
            # ESTRATÉGIA AMPLIADA: Busca o termo livremente dentro do organismo humano
            query = f'("{gene}") AND (organism_id:9606)'
            params = {
                "query": query,
                "fields": "accession,protein_name,organism_name",
                "format": "json",
            }

            response = requests.get(url, params=params)

            if response.status_code == 200:
                data = response.json()
                results = data.get("results", [])

                if results:
                    uniprot_id = results[0]["primaryAccession"]

                    desc = results[0].get("proteinDescription", {})
                    rec_name = (
                        desc.get("recommendedName")
                        or desc.get("submissionNames", [{}])[0]
                        or desc.get("submittedName", [{}])[0]
                    )
                    protein_name = rec_name.get("fullName", {}).get(
                        "value", "Nome Indisponível"
                    )
                    organism = results[0].get("organism", {}).get(
                        "scientificName", "Homo sapiens"
                    )

                    out.write(
                        f"{gene}\t{uniprot_id}\t{protein_name}\t{organism}\tRecuperado\n"
                    )
                    print(f"[RECUPERADO] {gene} -> {uniprot_id}")
                else:
                    out.write(f"{gene}\tNA\tNA\tNA\tConfirmado Ausente\n")
                    print(f"[AINDA AUSENTE] {gene}")
                    genes_ainda_ausentes.append(gene)
            else:
                print(f"[ERRO API] {gene}")
                genes_ainda_ausentes.append(gene)

        except Exception as e:
            print(f"[ERRO] {e} em {gene}")
            genes_ainda_ausentes.append(gene)

        time.sleep(0.2)

# Atualiza a lista dos que realmente não existem no UniProt
with open(
    "genes_definitivamente_ausentes.txt", "w", encoding="utf-8"
) as f_out:
    for g in genes_ainda_ausentes:
        f_out.write(f"{g}\n")

print("\n" + "=" * 40)
print(f"Total repesquisado: {len(genes_faltantes)}")
print(
    f"Recuperados com sucesso: {len(genes_faltantes) - len(genes_ainda_ausentes)}"
)
print(f"Permanecem não encontrados: {len(genes_ainda_ausentes)}")
print("Resultados salvos em 'recuperados_uniprot.tsv'")