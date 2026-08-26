import time
import requests

# 1. Carrega os genes do arquivo genes.txt
with open("genes.txt", "r", encoding="utf-8") as f:
    genes = [line.strip() for line in f if line.strip()]

genes_faltantes = []

url = "https://rest.uniprot.org/uniprotkb/search"

# 2. Executa a busca e grava o arquivo TSV
with open("resultados_uniprot.tsv", "w", encoding="utf-8") as out:
    out.write(
        "Gene_Original\tID_UniProt\tNome_Proteina\tOrganismo\tStatus\n"
    )

    for gene in genes:
        try:
            # Query flexível com taxID humano (9606)
            query = f"(gene:{gene}) AND (organism_id:9606)"
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
                    protein_name = results[0]["proteinDescription"][
                        "recommendedName"
                    ]["fullName"]["value"]
                    organism = results[0]["organism"]["scientificName"]

                    out.write(
                        f"{gene}\t{uniprot_id}\t{protein_name}\t{organism}\tEncontrado\n"
                    )
                    print(f"[SUCESSO] {gene} -> {uniprot_id}")
                else:
                    out.write(f"{gene}\tNA\tNA\tNA\tNão Encontrado\n")
                    print(f"[AVISO] {gene} não encontrado.")
                    genes_faltantes.append(gene)
            else:
                print(f"[ERRO] Falha na API para {gene}")
                genes_faltantes.append(gene)

        except Exception as e:
            print(f"[ERRO] {e} em {gene}")
            genes_faltantes.append(gene)

        time.sleep(0.2)

# 3. Salva os genes faltantes em um novo arquivo TXT
with open("genes_nao_identificados.txt", "w", encoding="utf-8") as f_out:
    for g in genes_faltantes:
        f_out.write(f"{g}\n")

# 4. Imprime o resumo final no terminal
print("\n" + "=" * 40)
print("RELATÓRIO DE EXECUÇÃO")
print("=" * 40)
print(f"Total de genes processados: {len(genes)}")
print(f"Total de não encontrados: {len(genes_faltantes)}")
print(
    "A lista dos genes ausentes foi salva em: 'genes_nao_identificados.txt'\n"
)