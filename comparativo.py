import pandas as pd


def consolidar_e_exibir_tabela():
    print("Iniciando consolidação e gerando tabela comparativa...\n")

    try:
        df_dbsnp = pd.read_csv("dbsnp_missense_resultados.tsv", sep="\t")
        df_am = pd.read_csv("alphamissense_resultados.tsv", sep="\t")
    except FileNotFoundError as e:
        print(f"[ERRO] Arquivo não encontrado: {e}")
        return

    # Normalização dos identificadores
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

    # Merge completo
    df_merged = pd.merge(
        dbsnp_validos,
        am_validos,
        on=["Gene", "rsID_dbSNP"],
        how="outer",
        suffixes=("_dbSNP", "_AlphaMissense"),
    )

    # Status de presença
    em_dbsnp = df_merged["rsID_dbSNP"].isin(dbsnp_validos["rsID_dbSNP"])
    em_am = df_merged["rsID_dbSNP"].isin(am_validos["rsID_dbSNP"])

    df_merged["Presenca_Fonte"] = "Desconhecido"
    df_merged.loc[em_dbsnp & em_am, "Presenca_Fonte"] = "Em Ambos"
    df_merged.loc[em_dbsnp & ~em_am, "Presenca_Fonte"] = "Apenas dbSNP"
    df_merged.loc[~em_dbsnp & em_am, "Presenca_Fonte"] = "Apenas AlphaMissense"

    # Seleciona 5 colunas chave para montar a amostra
    cols_chaves = ["Gene", "rsID_dbSNP", "Presenca_Fonte"]
    outras_cols = [c for c in df_merged.columns if c not in cols_chaves]

    # Para a amostra do terminal, pega linhas que tenham dados preenchidos
    df_preenchidos = df_merged[df_merged["Presenca_Fonte"] == "Em Ambos"].dropna(
        subset=outras_cols[:2]
    )

    # Se não houver totalmente preenchidos, usa o merged normal
    if df_preenchidos.empty:
        df_preenchidos = df_merged

    df_view = df_preenchidos[cols_chaves + outras_cols[:2]].head(8).copy()
    df_view = df_view.fillna("-")

    # Renomeia colunas para nomes curtos para caber no terminal
    df_view.columns = [str(c)[:12] for c in df_view.columns]

    # Força impressão sem quebras
    pd.set_option("display.max_columns", None)
    pd.set_option("display.expand_frame_repr", False)
    pd.set_option("display.width", 1000)

    print("=" * 75)
    print("AMOSTRA DA TABELA (COM CABEÇALHOS E COLUNAS ALINHADISAS)")
    print("=" * 75)
    print(df_view.to_string(index=False))
    print("=" * 75)

    # Exportação completa do dataset de 197.219 linhas
    output_name = "tabela_comparativa_completa.tsv"
    df_merged.to_csv(output_name, sep="\t", index=False)
    print(
        f"\n[SUCESSO] Dataset completo ({len(df_merged)} linhas) salvo em '{output_name}'!"
    )


if __name__ == "__main__":
    consolidar_e_exibir_tabela()