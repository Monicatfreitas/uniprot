import os
import re
import urllib.parse
import urllib.request
import json
import pandas as pd

AA_PROPS = {
    "Arg": "Básico/Carga Positiva",
    "Lys": "Básico/Carga Positiva",
    "His": "Básico/Carga Positiva",
    "Asp": "Ácido/Carga Negativa",
    "Glu": "Ácido/Carga Negativa",
    "Ser": "Polar Neutro",
    "Thr": "Polar Neutro",
    "Asn": "Polar Neutro",
    "Gln": "Polar Neutro",
    "Cys": "Especial/Reativo",
    "Gly": "Especial/Flexível",
    "Pro": "Especial/Rígido",
    "Ala": "Apolar/Hidrofóbico",
    "Val": "Apolar/Hidrofóbico",
    "Ile": "Apolar/Hidrofóbico",
    "Leu": "Apolar/Hidrofóbico",
    "Met": "Apolar/Hidrofóbico",
    "Phe": "Apolar/Aromático",
    "Tyr": "Apolar/Aromático",
    "Trp": "Apolar/Aromático",
}


def buscar_nomes_proteinas_uniprot(lista_genes):
    """Consulta a API do UniProtKB para obter o nome oficial da proteína por gene (Humano)."""
    mapa_proteinas = {}
    genes_unicos = [
        g
        for g in set(lista_genes)
        if pd.notna(g) and str(g).strip() not in ["", "nan", "NA"]
    ]

    print(f"Buscando nomes oficiais de proteínas no UniProt para {len(genes_unicos)} genes...")

    # Consulta em blocos para otimizar a requisição
    for gene in genes_unicos:
        try:
            query = f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true"
            url = f"https://rest.uniprot.org/uniprotkb/search?query={urllib.parse.quote(query)}&fields=protein_name&size=1"

            req = urllib.request.Request(
                url, headers={"User-Agent": "Python/BioinformaticsScript"}
            )
            with urllib.request.urlopen(req, timeout=5) as response:
                data = json.loads(response.read().decode())
                results = data.get("results", [])

                if results:
                    nome_rec = (
                        results[0]
                        .get("proteinDescription", {})
                        .get("recommendedName", {})
                        .get("fullName", {})
                        .get("value")
                    )
                    if nome_rec:
                        mapa_proteinas[gene] = nome_rec
                        continue

            mapa_proteinas[gene] = "-"
        except Exception:
            mapa_proteinas[gene] = "-"

    return mapa_proteinas


def analisar_tipo_modificacao(prot_var):
    if str(prot_var) in ["-", "nan", "NA", "", "None"]:
        return "-", "Desconhecido"

    match_3 = re.search(r"p\.([A-Z][a-z]{2})\d+([A-Z][a-z]{2})", str(prot_var))
    if match_3:
        aa_ref, aa_alt = match_3.group(1), match_3.group(2)
        prop_ref = AA_PROPS.get(aa_ref, "Desconhecido")
        prop_alt = AA_PROPS.get(aa_alt, "Desconhecido")

        if prop_ref == prop_alt:
            tipo = f"Conservativa ({prop_ref})"
        else:
            tipo = f"Não-Conservativa ({prop_ref} -> {prop_alt})"

        return f"{aa_ref} -> {aa_alt}", tipo

    match_1 = re.search(r"([A-Z])\d+([A-Z])", str(prot_var))
    if match_1:
        return f"{match_1.group(1)} -> {match_1.group(2)}", "Substituição Missense"

    return "-", "Substituição Missense"


def gerar_tabela_completa_proteina_tsv():
    print("Iniciando processamento exaustivo com API UniProt...")

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

    dbsnp_validos = df_dbsnp[
        df_dbsnp["rsID_dbSNP"].notna()
        & (~df_dbsnp["rsID_dbSNP"].isin(["NA", "nan", ""]))
    ].copy()

    am_validos = df_am[
        df_am["rsID_dbSNP"].notna()
        & (~df_am["rsID_dbSNP"].isin(["NA", "nan", ""]))
    ].copy()

    df_merged = pd.merge(
        dbsnp_validos,
        am_validos,
        on=["Gene", "rsID_dbSNP"],
        how="outer",
        suffixes=("_dbSNP", "_AlphaMissense"),
    )

    # 1. Busca dinâmica do Nome Oficial da Proteína no UniProt
    mapa_nomes = buscar_nomes_proteinas_uniprot(df_merged["Gene"].unique())
    df_merged["Nome_Oficial_Proteina"] = df_merged["Gene"].map(mapa_nomes).fillna("-")

    # 2. Resgate de Variantes de Proteína
    col_var_prot = [
        c
        for c in df_merged.columns
        if any(
            k in c.lower()
            for k in [
                "hgvs.p",
                "protein_variant",
                "aa_change",
                "prot_change",
                "protein",
            ]
        )
    ]
    var_prot_col = (
        col_var_prot[0] if col_var_prot else "Proteina_Alteracao_Original"
    )

    if var_prot_col not in df_merged.columns:
        df_merged[var_prot_col] = "-"

    # 3. Análise da Modificação
    trocas_e_tipos = df_merged[var_prot_col].apply(analisar_tipo_modificacao)
    df_merged["Troca_Aminoacido"] = [t[0] for t in trocas_e_tipos]
    df_merged["Tipo_Modificacao_FisicoQuimica"] = [t[1] for t in trocas_e_tipos]

    # Status e Links
    em_dbsnp = df_merged["rsID_dbSNP"].isin(dbsnp_validos["rsID_dbSNP"])
    em_am = df_merged["rsID_dbSNP"].isin(am_validos["rsID_dbSNP"])

    df_merged["Presenca_Fonte"] = "Desconhecido"
    df_merged.loc[em_dbsnp & em_am, "Presenca_Fonte"] = "Em Ambos"
    df_merged.loc[em_dbsnp & ~em_am, "Presenca_Fonte"] = "Apenas dbSNP"
    df_merged.loc[~em_dbsnp & em_am, "Presenca_Fonte"] = "Apenas AlphaMissense"

    df_merged["Link_dbSNP"] = df_merged["rsID_dbSNP"].apply(
        lambda x: (
            f"https://www.ncbi.nlm.nih.gov/snp/{str(x).strip()}"
            if str(x).startswith("rs")
            else "-"
        )
    )

    # Reordenação
    cols_prioridade = [
        "Gene",
        "Nome_Oficial_Proteina",
        "rsID_dbSNP",
        var_prot_col,
        "Troca_Aminoacido",
        "Tipo_Modificacao_FisicoQuimica",
        "Link_dbSNP",
        "Presenca_Fonte",
    ]

    cols_existentes = [c for c in cols_prioridade if c in df_merged.columns]
    outras_cols = [c for c in df_merged.columns if c not in cols_existentes]
    df_merged = df_merged[cols_existentes + outras_cols]

    output_path = os.path.join(
        pasta_resultados, "tabela_comparativa_detalhada_proteinas.tsv"
    )
    df_merged.to_csv(output_path, sep="\t", index=False)
    print(
        f"[SUCESSO] Tabela salva com nomes oficiais de proteínas em: '{output_path}'!"
    )


if __name__ == "__main__":
    gerar_tabela_completa_proteina_tsv()