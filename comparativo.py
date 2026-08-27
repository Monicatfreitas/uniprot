import os
import re
import urllib.parse
import urllib.request
import json
import pandas as pd

AA_PROPS = {
    "Arg": "Básico/Carga Positiva", "Lys": "Básico/Carga Positiva", "His": "Básico/Carga Positiva",
    "Asp": "Ácido/Carga Negativa", "Glu": "Ácido/Carga Negativa",
    "Ser": "Polar Neutro", "Thr": "Polar Neutro", "Asn": "Polar Neutro", "Gln": "Polar Neutro",
    "Cys": "Especial/Reativo", "Gly": "Especial/Flexível", "Pro": "Especial/Rígido",
    "Ala": "Apolar/Hidrofóbico", "Val": "Apolar/Hidrofóbico", "Ile": "Apolar/Hidrofóbico",
    "Leu": "Apolar/Hidrofóbico", "Met": "Apolar/Hidrofóbico",
    "Phe": "Apolar/Aromático", "Tyr": "Apolar/Aromático", "Trp": "Apolar/Aromático",
}

# Mapeamento para conversão de código de 1 letra para 3 letras (ex: R175H -> Arg175His)
AA_1_TO_3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "E": "Glu", "Q": "Gln",
    "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"
}

def buscar_nomes_proteinas_uniprot(lista_genes):
    mapa_proteinas = {}
    genes_unicos = [g for g in set(lista_genes) if pd.notna(g) and str(g).strip() not in ["", "nan", "NA"]]
    print(f"Buscando nomes oficiais de proteínas no UniProt para {len(genes_unicos)} genes...")
    
    for gene in genes_unicos:
        try:
            query = f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true"
            url = f"https://rest.uniprot.org/uniprotkb/search?query={urllib.parse.quote(query)}&fields=protein_name&size=1"
            req = urllib.request.Request(url, headers={'User-Agent': 'Python/BioinformaticsScript'})
            with urllib.request.urlopen(req, timeout=5) as response:
                data = json.loads(response.read().decode())
                results = data.get("results", [])
                if results:
                    nome_rec = results[0].get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")
                    if nome_rec:
                        mapa_proteinas[gene] = nome_rec
                        continue
            mapa_proteinas[gene] = "-"
        except Exception:
            mapa_proteinas[gene] = "-"
    return mapa_proteinas

def extrair_e_classificar_aminoacido(row, colunas_candidatas):
    """Procura a variante em várias colunas possíveis e extrai a mudança e seu impacto."""
    var_str = "-"
    for col in colunas_candidatas:
        if col in row and pd.notna(row[col]):
            v = str(row[col]).strip()
            if v not in ["-", "nan", "NA", "", "None"]:
                var_str = v
                break

    if var_str == "-":
        return "-", "-", "Não Informado"

    # Caso 1: Código de 3 letras (ex: p.Arg175His ou Arg175His)
    match_3 = re.search(r"(?:p\.)?([A-Z][a-z]{2})\d+([A-Z][a-z]{2})", var_str)
    if match_3:
        aa_ref, aa_alt = match_3.group(1), match_3.group(2)
        prop_ref = AA_PROPS.get(aa_ref, "Desconhecido")
        prop_alt = AA_PROPS.get(aa_alt, "Desconhecido")
        tipo = f"Conservativa ({prop_ref})" if prop_ref == prop_alt else f"Não-Conservativa ({prop_ref} -> {prop_alt})"
        return var_str, f"{aa_ref} -> {aa_alt}", tipo

    # Caso 2: Código de 1 letra (ex: R175H ou p.R175H)
    match_1 = re.search(r"(?:p\.)?([A-Z])\d+([A-Z])", var_str)
    if match_1:
        a1, a2 = match_1.group(1), match_1.group(2)
        aa_ref = AA_1_TO_3.get(a1, a1)
        aa_alt = AA_1_TO_3.get(a2, a2)
        prop_ref = AA_PROPS.get(aa_ref, "Desconhecido")
        prop_alt = AA_PROPS.get(aa_alt, "Desconhecido")
        tipo = f"Conservativa ({prop_ref})" if prop_ref == prop_alt else f"Não-Conservativa ({prop_ref} -> {prop_alt})"
        return var_str, f"{aa_ref} -> {aa_alt} ({a1}>{a2})", tipo

    return var_str, "-", "Substituição Missense"

def gerar_tabela_completa_proteina_tsv():
    print("Iniciando reconstrução dinâmica de aminoácidos e atributos...")
    pasta_resultados = "resultados"
    os.makedirs(pasta_resultados, exist_ok=True)

    try:
        df_dbsnp = pd.read_csv("dbsnp_missense_resultados.tsv", sep="\t")
        df_am = pd.read_csv("alphamissense_resultados.tsv", sep="\t")
    except FileNotFoundError as e:
        print(f"[ERRO] Arquivo não encontrado: {e}")
        return

    for df in [df_dbsnp, df_am]:
        if "rsID_dbSNP" in df.columns:
            df["rsID_dbSNP"] = df["rsID_dbSNP"].astype(str).str.strip()
        if "Gene" in df.columns:
            df["Gene"] = df["Gene"].astype(str).str.strip()

    dbsnp_validos = df_dbsnp[df_dbsnp["rsID_dbSNP"].notna() & (~df_dbsnp["rsID_dbSNP"].isin(["NA", "nan", ""]))].copy()
    am_validos = df_am[df_am["rsID_dbSNP"].notna() & (~df_am["rsID_dbSNP"].isin(["NA", "nan", ""]))].copy()

    df_merged = pd.merge(dbsnp_validos, am_validos, on=["Gene", "rsID_dbSNP"], how="outer", suffixes=("_dbSNP", "_AlphaMissense"))

    # 1. Nome da Proteína via UniProt
    mapa_nomes = buscar_nomes_proteinas_uniprot(df_merged["Gene"].unique())
    df_merged["Nome_Oficial_Proteina"] = df_merged["Gene"].map(mapa_nomes).fillna("-")

    # 2. Varredura de Colunas Candidatas a conter a variante de proteína
    colunas_candidatas = [
        c for c in df_merged.columns 
        if any(k in c.lower() for k in ["protein_variant", "hgvs.p", "aa_change", "prot_change", "protein", "variant", "mutation"])
    ]

    res_aa = df_merged.apply(lambda r: extrair_e_classificar_aminoacido(r, colunas_candidatas), axis=1)
    df_merged["Proteina_Mudanca"] = [r[0] for r in res_aa]
    df_merged["Troca_Aminoacido"] = [r[1] for r in res_aa]
    df_merged["Tipo_Modificacao_FisicoQuimica"] = [r[2] for r in res_aa]

    # Presença de Fonte e Links
    em_dbsnp = df_merged["rsID_dbSNP"].isin(dbsnp_validos["rsID_dbSNP"])
    em_am = df_merged["rsID_dbSNP"].isin(am_validos["rsID_dbSNP"])
    df_merged["Presenca_Fonte"] = "Desconhecido"
    df_merged.loc[em_dbsnp & em_am, "Presenca_Fonte"] = "Em Ambos"
    df_merged.loc[em_dbsnp & ~em_am, "Presenca_Fonte"] = "Apenas dbSNP"
    df_merged.loc[~em_dbsnp & em_am, "Presenca_Fonte"] = "Apenas AlphaMissense"

    df_merged["Link_dbSNP"] = df_merged["rsID_dbSNP"].apply(
        lambda x: f"https://www.ncbi.nlm.nih.gov/snp/{str(x).strip()}" if str(x).startswith("rs") else "-"
    )

    # Reordenação
    cols_prioridade = [
        "Gene", "Nome_Oficial_Proteina", "rsID_dbSNP", "Proteina_Mudanca", 
        "Troca_Aminoacido", "Tipo_Modificacao_FisicoQuimica", "Link_dbSNP", "Presenca_Fonte"
    ]
    cols_existentes = [c for c in cols_prioridade if c in df_merged.columns]
    outras_cols = [c for c in df_merged.columns if c not in cols_existentes]
    df_merged = df_merged[cols_existentes + outras_cols]

    output_path = os.path.join(pasta_resultados, "tabela_comparativa_detalhada_proteinas.tsv")
    df_merged.to_csv(output_path, sep="\t", index=False)
    print(f"[SUCESSO] Tabela atualizada e reconstruída em: '{output_path}'!")

if __name__ == "__main__":
    gerar_tabela_completa_proteina_tsv()