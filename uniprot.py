import requests
import time

print("Lendo o arquivo de genes...")
with open("genes.txt", "r") as f:
    genes = [line.strip() for line in f if line.strip()]

print(f"{len(genes)} genes encontrados. Consultando o UniProt...")

with open("resultados_uniprot.tsv", "w") as out:
    out.write("Gene_Pesquisado\tUniProt_ID\tProtein_Name\tOrganism\tLength\n")

    for gene in genes:
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+AND+organism_id:9606&size=1"
        
        try:
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                results = data.get("results", [])
                
                if results:
                    entry = results[0]
                    uniprot_id = entry.get("primaryAccession", "N/A")
                    protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "N/A")
                    organism = entry.get("organism", {}).get("scientificName", "N/A")
                    length = str(entry.get("sequence", {}).get("length", "N/A"))
                    
                    out.write(f"{gene}\t{uniprot_id}\t{protein_name}\t{organism}\t{length}\n")
                    print(f"[SUCESSO] {gene} -> {uniprot_id}")
                else:
                    out.write(f"{gene}\tNA\tNA\tNA\tNA\n")
                    print(f"[AVISO] {gene} não encontrado.")
            else:
                print(f"[ERRO] Falha na API para {gene}")
        except Exception as e:
            print(f"[ERRO] {e} em {gene}")
            
        time.sleep(0.3)

print("\nPronto! Arquivo 'resultados_uniprot.tsv' gerado com sucesso.")