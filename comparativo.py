import os
import re
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
    print("Iniciando mapeamento detalhado de proteína e modificações...")

    pasta_resultados = "resultados"
    os.makedirs(pasta_resultados, exist_ok=True)

    try:
        df_dbsnp = pd.read_csv("dbsnp_missense_resultados.tsv", sep="\t")
        df_am = pd.read_csv("alphamissense_resultados.tsv", sep="\t")
    except FileNotFoundError as e:
        print(f"[ERRO] Arquivo não encontrado: {e}")
        return

    # Normalização de chaves
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

    # 1. Nome da Proteína / ID UniProt
    col_uniprot = [
        c
        for c in df_merged.columns
        if any(
            k in c.lower()
            for k in ["uniprot", "protein_name", "protein_id", "transcript"]
        )
    ]
    if col_uniprot:
        df_merged["Nome_Proteina_UniProt"] = df_merged[col_uniprot[0]].fillna(
            "-"
        )
    else:
        df_merged["Nome_Proteina_UniProt"] = df_merged["Gene"].apply(
            lambda g: f"Proteína de {g}" if pd.notna(g) else "-"
        )

    # 2. Busca da coluna de variante da proteína com fallback seguro
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

    # 3. Aplicação da Análise Físico-Química
    trocas_e_tipos = df_merged[var_prot_col].apply(analisar_tipo_modificacao)
    df_merged["Troca_Aminoacido"] = [t[0] for t in trocas_e_tipos]
    df_merged["Tipo_Modificacao_FisicoQuimica"] = [t[1] for t in trocas_e_tipos]

    # Status de Fonte
    em_dbsnp = df_merged["rsID_dbSNP"].isin(dbsnp_validos["rsID_dbSNP"])
    em_am = df_merged["rsID_dbSNP"].isin(am_validos["rsID_dbSNP"])

    df_merged["Presenca_Fonte"] = "Desconhecido"
    df_merged.loc[em_dbsnp & em_am, "Presenca_Fonte"] = "Em Ambos"
    df_merged.loc[em_dbsnp & ~em_am, "Presenca_Fonte"] = "Apenas dbSNP"
    df_merged.loc[~em_dbsnp & em_am, "Presenca_Fonte"] = "Apenas AlphaMissense"

    # Link Direto
    df_merged["Link_dbSNP"] = df_merged["rsID_dbSNP"].apply(
        lambda x: (
            f"https://www.ncbi.nlm.nih.gov/snp/{str(x).strip()}"
            if str(x).startswith("rs")
            else "-"
        )
    )

    # Reordenação Prioritária
    cols_prioridade = [
        "Gene",
        "Nome_Proteina_UniProt",
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

    # Salva no arquivo final em resultados/
    output_path = os.path.join(
        pasta_resultados, "tabela_comparativa_detalhada_proteinas.tsv"
    )
    df_merged.to_csv(output_path, sep="\t", index=False)
    print(
        f"[SUCESSO] Arquivo salvo sem erros em '{output_path}' ({len(df_merged)} linhas)!"
    )


if __name__ == "__main__":
    gerar_tabela_completa_proteina_tsv()