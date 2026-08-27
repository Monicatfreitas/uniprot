import os
import pandas as pd


def gerar_tabela_completa_e_exaustiva():
    print("Iniciando consolidação exaustiva de todos os atributos...")

    # 1. Definição e criação da pasta de saída
    pasta_resultados = "resultados"
    os.makedirs(pasta_resultados, exist_ok=True)

    try:
        df_dbsnp = pd.read_csv("dbsnp_missense_resultados.tsv", sep="\t")
        df_am = pd.read_csv("alphamissense_resultados.tsv", sep="\t")
    except FileNotFoundError as e:
        print(f"[ERRO] Arquivo não encontrado: {e}")
        return

    # Normalização de IDs
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

    # Merge Completo (Outer Merge)
    df_merged = pd.merge(
        dbsnp_validos,
        am_validos,
        on=["Gene", "rsID_dbSNP"],
        how="outer",
        suffixes=("_dbSNP", "_AlphaMissense"),
    )

    # 2. Presença de Fonte
    em_dbsnp = df_merged["rsID_dbSNP"].isin(dbsnp_validos["rsID_dbSNP"])
    em_am = df_merged["rsID_dbSNP"].isin(am_validos["rsID_dbSNP"])

    df_merged["Presenca_Fonte"] = "Desconhecido"
    df_merged.loc[em_dbsnp & em_am, "Presenca_Fonte"] = "Em Ambos"
    df_merged.loc[em_dbsnp & ~em_am, "Presenca_Fonte"] = "Apenas dbSNP"
    df_merged.loc[~em_dbsnp & em_am, "Presenca_Fonte"] = "Apenas AlphaMissense"

    # Mapeamento dinâmico de colunas
    col_score = [
        c
        for c in df_merged.columns
        if "am_pathogenicity" in c.lower() or "score" in c.lower()
    ]
    col_am_class = [c for c in df_merged.columns if "am_class" in c.lower()]
    col_clin = [
        c
        for c in df_merged.columns
        if "clinical" in c.lower() or "significance" in c.lower()
    ]

    score_col = col_score[0] if col_score else None
    am_class_col = col_am_class[0] if col_am_class else None
    clin_col = col_clin[0] if col_clin else None

    # 3. Categorização do Score
    def categorizar_score(val):
        try:
            s = float(val)
            if s >= 0.564:
                return (
                    "Patogênico Forte (>=0.74)"
                    if s >= 0.74
                    else "Provável Patogênico (0.564-0.74)"
                )
            elif s <= 0.34:
                return "Provável Benigno (<=0.34)"
            else:
                return "Zona Ambígua (0.34-0.564)"
        except (ValueError, TypeError):
            return "Sem Score"

    if score_col:
        df_merged["AM_Score_Faixa"] = df_merged[score_col].apply(
            categorizar_score
        )

    # 4. Padrão Biológico e Alerta
    def classificar_status_e_alerta(row):
        fonte = row["Presenca_Fonte"]
        if fonte != "Em Ambos":
            return fonte, False

        val_clin = str(row[clin_col]).lower() if clin_col else ""
        val_am = str(row[am_class_col]).lower() if am_class_col else ""
        s_am = row[score_col] if score_col else None

        c_path = any(
            x in val_clin for x in ["pathogenic", "likely_pathogenic"]
        )
        c_benign = any(
            x in val_clin for x in ["benign", "likely_benign"]
        )
        c_vus = "uncertain" in val_clin or "vus" in val_clin or val_clin == ""

        am_path = "likely_pathogenic" in val_am or "pathogenic" in val_am
        am_benign = "likely_benign" in val_am or "benign" in val_am

        if (c_path and am_benign) or (c_benign and am_path):
            return "Discordância Direta (ClinVar vs AM)", True

        if c_vus and s_am is not None:
            try:
                score_num = float(s_am)
                if score_num >= 0.564:
                    return (
                        "Reclassificação: VUS ClinVar -> AM Patogênico",
                        False,
                    )
                elif score_num <= 0.34:
                    return "Reclassificação: VUS ClinVar -> AM Benigno", False
            except (ValueError, TypeError):
                pass

        if c_path and am_path:
            return "Concordante (Patogênico)", False
        if c_benign and am_benign:
            return "Concordante (Benigno)", False

        return "Outros / Inconclusivos", False

    res = df_merged.apply(classificar_status_e_alerta, axis=1)
    df_merged["Padrao_Comparacao"] = [r[0] for r in res]
    df_merged["Alerta_Conflito_Direto"] = [r[1] for r in res]

    # Reordenação
    cols_prioridade = [
        "Gene",
        "rsID_dbSNP",
        "Presenca_Fonte",
        "Padrao_Comparacao",
        "Alerta_Conflito_Direto",
    ]
    if score_col:
        cols_prioridade.extend([score_col, "AM_Score_Faixa"])

    outras_colunas = [c for c in df_merged.columns if c not in cols_prioridade]
    df_merged = df_merged[cols_prioridade + outras_colunas]

    # 5. Exportação para dentro do novo diretório
    output_path = os.path.join(
        pasta_resultados, "tabela_comparativa_completa.tsv"
    )
    df_merged.to_csv(output_path, sep="\t", index=False)
    print(
        f"[SUCESSO] Tabela salva em: '{output_path}' ({len(df_merged)} linhas e {len(df_merged.columns)} colunas)!"
    )


if __name__ == "__main__":
    gerar_tabela_completa_e_exaustiva()