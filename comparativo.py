import os
from typing import Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", palette="muted")
plt.rcParams.update({
    "font.sans-serif": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "figure.titlesize": 14,
    "figure.titleweight": "bold"
})

PASTA_RESULTADOS = "resultados"
PASTA_FIGURAS_GENES = os.path.join(PASTA_RESULTADOS, "figuras_por_gene")
PASTA_RELATORIOS_GENES = os.path.join(PASTA_RESULTADOS, "relatorios_por_gene")

ARQUIVOS_DBSNP = [
    "comparacao_dbsnp_alphamissense.tsv",
    "dbsnp_missense_resultados.tsv",
    "dbsnp_alphamissense_unificado.tsv",
    "tabela_comparativa_completa.tsv"
]

ARQUIVOS_EXTRAS = {
    "alphamissense": ["alphamissense_resultados.tsv"],
    "uniprot": ["resultados_uniprot.tsv"]
}

# Base do link do AlphaFold (AlphaFold DB usa o UniProt accession na URL)
URL_BASE_ALPHAFOLD = "https://alphafold.ebi.ac.uk/entry/"


def limpar_e_converter_numerico(serie: pd.Series) -> pd.Series:
    """Converte e limpa qualquer formato decimal (ponto ou vírgula)."""
    s = serie.astype(str).str.replace(',', '.').str.strip()
    s = s.replace(['nan', 'none', 'n/a', '-', 'null', '', '<na>', 'None'], np.nan)
    return pd.to_numeric(s, errors="coerce")


def buscar_e_carregar_tsv(lista_nomes: list, obrigatorio: bool = True) -> Optional[pd.DataFrame]:
    """Procura o arquivo tanto na raiz quanto dentro de resultados/."""
    for nome in lista_nomes:
        caminhos_para_testar = [nome, os.path.join(PASTA_RESULTADOS, nome)]
        for caminho in caminhos_para_testar:
            if os.path.exists(caminho):
                try:
                    df = pd.read_csv(caminho, sep="\t", engine="python", on_bad_lines="skip")
                    df.columns = df.columns.str.strip()
                    print(f"[OK] Arquivo carregado com sucesso: {caminho}")
                    return df
                except Exception as e:
                    print(f"[AVISO] Falha ao ler {caminho}: {e}")
                    continue

    if obrigatorio:
        print(f"[ERRO] Nenhum arquivo válido encontrado para a lista: {lista_nomes}")
    return None


def tratar_colunas_base(df: pd.DataFrame) -> pd.DataFrame:
    """Normaliza o nome do gene (mantendo em maiúsculas) e limpa os rsIDs."""
    if "Gene_Original" in df.columns and "Gene" not in df.columns:
        df = df.rename(columns={"Gene_Original": "Gene"})

    for col in df.columns:
        if col.lower() in ["gene", "gene_name", "symbol"]:
            df[col] = df[col].astype(str).str.strip().str.upper()
        elif "rs" in col.lower() or "dbsnp" in col.lower():
            df['rs_key'] = df[col].astype(str).str.replace("rs", "", case=False).str.strip()

    return df


def _detectar_coluna(df: pd.DataFrame, chaves: list, excluir_prefixo: Optional[str] = None) -> Optional[str]:
    """Procura a primeira coluna cujo nome (em minúsculas) contenha alguma das chaves."""
    for c in df.columns:
        cl = c.lower()
        if excluir_prefixo and cl.startswith(excluir_prefixo):
            continue
        if any(k in cl for k in chaves):
            return c
    return None


def carregar_dados_dbsnp():
    # 1. Carrega a tabela principal
    df_dbsnp = buscar_e_carregar_tsv(ARQUIVOS_DBSNP, obrigatorio=True)
    if df_dbsnp is None:
        return None

    df_dbsnp = tratar_colunas_base(df_dbsnp)
    col_gene = next((c for c in df_dbsnp.columns if "gene" in c.lower()), "Gene")

    # 2. Carrega o arquivo do AlphaMissense
    df_am = buscar_e_carregar_tsv(ARQUIVOS_EXTRAS["alphamissense"], obrigatorio=False)

    if df_am is not None:
        df_am = tratar_colunas_base(df_am)

        candidatos_score = [c for c in df_am.columns if any(k in c.lower() for k in ["pathogenicity", "score", "am_score", "alphamissense"])]
        col_score_am = candidatos_score[0] if candidatos_score else None

        # DIAGNÓSTICO 1: mostra TODAS as colunas candidatas a "score" no arquivo do
        # AlphaMissense. Se houver mais de uma, a escolhida (a primeira) pode não ser
        # a certa — nesse caso, ajuste manualmente qual usar.
        print(f"[DIAGNÓSTICO] Colunas candidatas a score no AlphaMissense: {candidatos_score}")
        print(f"[DIAGNÓSTICO] Coluna escolhida para score: {col_score_am!r}")

        if col_score_am:
            amostra_bruta = df_am[col_score_am].dropna().astype(str).head(5).tolist()
            print(f"[DIAGNÓSTICO] Amostra de valores BRUTOS (antes de converter) em {col_score_am!r}: {amostra_bruta}")
            df_am[col_score_am] = limpar_e_converter_numerico(df_am[col_score_am])
            print(f"[DIAGNÓSTICO] Após conversão numérica: {df_am[col_score_am].notna().sum()} válidos de {len(df_am)} no arquivo AlphaMissense")

        usa_rs_key = 'rs_key' in df_dbsnp.columns and 'rs_key' in df_am.columns
        chave_merge = [col_gene, 'rs_key'] if usa_rs_key else [col_gene]

        # DIAGNÓSTICO 2: compara os valores de chave dos dois lados ANTES do merge.
        genes_dbsnp = set(df_dbsnp[col_gene].dropna().unique())
        genes_am = set(df_am[col_gene].dropna().unique())
        print(f"[DIAGNÓSTICO] Genes únicos no dbSNP: {len(genes_dbsnp)} | no AlphaMissense: {len(genes_am)} | em comum: {len(genes_dbsnp & genes_am)}")
        if usa_rs_key:
            rs_dbsnp = set(df_dbsnp['rs_key'].dropna().unique())
            rs_am = set(df_am['rs_key'].dropna().unique())
            print(f"[DIAGNÓSTICO] rs_key únicos no dbSNP: {len(rs_dbsnp)} | no AlphaMissense: {len(rs_am)} | em comum: {len(rs_dbsnp & rs_am)}")
            print(f"[DIAGNÓSTICO] Exemplo rs_key dbSNP: {list(rs_dbsnp)[:5]} | Exemplo rs_key AlphaMissense: {list(rs_am)[:5]}")

        colunas_para_substituir = [c for c in df_am.columns if c in df_dbsnp.columns and c not in chave_merge]
        if colunas_para_substituir:
            print(f"[INFO] Substituindo colunas placeholder do dbSNP pelos dados reais do AlphaMissense: {colunas_para_substituir}")
            df_dbsnp = df_dbsnp.drop(columns=colunas_para_substituir)

        if usa_rs_key:
            df_dbsnp = pd.merge(df_dbsnp, df_am, on=[col_gene, 'rs_key'], how='left')
        else:
            df_dbsnp = pd.merge(df_dbsnp, df_am, on=col_gene, how='left')

        if col_score_am and col_score_am in df_dbsnp.columns:
            print(f"[DIAGNÓSTICO] Após o merge com AlphaMissense: {df_dbsnp[col_score_am].notna().sum()} de {len(df_dbsnp)} linhas com score preenchido")
    else:
        print("[AVISO] Arquivo do AlphaMissense não encontrado — seguindo só com dbSNP.")

    # 3. Carrega o arquivo do UniProt (código UniProt, posição no gene/proteína, etc.)
    #    Esse merge não existia antes, por isso UniProt_ID/posições/AlphaFold ficavam
    #    sempre "None" na tabela final, mesmo a coluna já existindo no esquema.
    df_uniprot = buscar_e_carregar_tsv(ARQUIVOS_EXTRAS["uniprot"], obrigatorio=False)

    if df_uniprot is not None:
        df_uniprot = tratar_colunas_base(df_uniprot)

        col_uniprot_id = _detectar_coluna(df_uniprot, ["uniprot", "accession", "swissprot"])
        col_pos_proteina = _detectar_coluna(df_uniprot, ["protein_pos", "posicao_proteina", "aa_position", "protein_start", "residue"])
        col_pos_gene = _detectar_coluna(df_uniprot, ["genomic_pos", "posicao_genomica", "posicao_gene", "gene_pos", "chrom_pos", "cds_pos"])
        col_mudanca_proteina = _detectar_coluna(df_uniprot, ["proteina_mudanca", "protein_change", "aa_change", "hgvs_p", "hgvsp"])

        # DIAGNÓSTICO 3: mostra o que foi detectado no arquivo de UniProt, já que os
        # nomes de coluna exatos desse arquivo não estavam disponíveis no momento do
        # ajuste — confira estas linhas na saída e me avise se algo bateu errado.
        print(f"[DIAGNÓSTICO] Colunas do arquivo UniProt: {list(df_uniprot.columns)}")
        print(f"[DIAGNÓSTICO] UniProt_ID detectado: {col_uniprot_id!r} | Posição na proteína: {col_pos_proteina!r} | "
              f"Posição no gene: {col_pos_gene!r} | Mudança na proteína: {col_mudanca_proteina!r}")

        usa_rs_key_up = 'rs_key' in df_dbsnp.columns and 'rs_key' in df_uniprot.columns
        chave_merge_up = [col_gene, 'rs_key'] if usa_rs_key_up else [col_gene]

        colunas_para_substituir_up = [c for c in df_uniprot.columns if c in df_dbsnp.columns and c not in chave_merge_up]
        if colunas_para_substituir_up:
            print(f"[INFO] Substituindo colunas placeholder do dbSNP pelos dados reais do UniProt: {colunas_para_substituir_up}")
            df_dbsnp = df_dbsnp.drop(columns=colunas_para_substituir_up)

        if usa_rs_key_up:
            df_dbsnp = pd.merge(df_dbsnp, df_uniprot, on=[col_gene, 'rs_key'], how='left')
        else:
            df_dbsnp = pd.merge(df_dbsnp, df_uniprot, on=col_gene, how='left')

        if col_uniprot_id and col_uniprot_id in df_dbsnp.columns:
            print(f"[DIAGNÓSTICO] Após o merge com UniProt: {df_dbsnp[col_uniprot_id].notna().sum()} de {len(df_dbsnp)} linhas com UniProt_ID preenchido")

            # 4. Gera o link para o AlphaFold a partir do UniProt_ID
            df_dbsnp["Link_AlphaFold"] = df_dbsnp[col_uniprot_id].apply(
                lambda uid: f"{URL_BASE_ALPHAFOLD}{str(uid).strip()}" if pd.notna(uid) and str(uid).strip().lower() not in ("", "nan", "none") else np.nan
            )
        else:
            print("[AVISO] Não foi possível localizar a coluna de UniProt_ID no arquivo de UniProt — link do AlphaFold não gerado. "
                  "Confira o nome real da coluna na lista impressa acima e ajuste `_detectar_coluna` se necessário.")
    else:
        print("[AVISO] Arquivo de UniProt (resultados_uniprot.tsv) não encontrado — UniProt_ID, "
              "posições e link do AlphaFold ficarão vazios.")

    return df_dbsnp


def gerar_figuras_e_relatorios(df_merged: pd.DataFrame):
    os.makedirs(PASTA_FIGURAS_GENES, exist_ok=True)
    os.makedirs(PASTA_RELATORIOS_GENES, exist_ok=True)

    col_gene = next((c for c in df_merged.columns if "gene" in c.lower()), None)
    if not col_gene:
        print("[ERRO] Coluna referente a 'Gene' não foi encontrada.")
        return

    col_score_am = next((c for c in df_merged.columns if any(k in c.lower() for k in ["pathogenicity", "score", "am_score", "alphamissense", "am_pathogenicity"])), None)
    col_class_dbsnp = next(
        (c for c in df_merged.columns
         if not c.lower().startswith("alphamissense")
         and any(k in c.lower() for k in ["classificacao", "clinical", "significance", "clinvar", "significancia"])),
        None
    )
    col_pos = next((c for c in df_merged.columns if any(k in c.lower() for k in ["pos", "protein_start", "position", "posicao"])), None)
    col_rs = next((c for c in df_merged.columns if "rs" in c.lower() or "dbsnp" in c.lower()), None)

    if col_score_am:
        df_merged[col_score_am] = limpar_e_converter_numerico(df_merged[col_score_am])
        print(f"[INFO] Usando coluna de score: {col_score_am!r} — {df_merged[col_score_am].notna().sum()} valores válidos de {len(df_merged)}")
    if col_class_dbsnp:
        print(f"[DIAGNÓSTICO] Coluna de classificação (dbSNP/ClinVar) escolhida: {col_class_dbsnp!r} — {df_merged[col_class_dbsnp].notna().sum()} valores válidos de {len(df_merged)}")
    else:
        print("[DIAGNÓSTICO] Nenhuma coluna de classificação dbSNP/ClinVar encontrada.")
    if col_pos:
        df_merged[col_pos] = limpar_e_converter_numerico(df_merged[col_pos])

    todos_os_genes = df_merged[col_gene].dropna().unique()
    total_genes = len(todos_os_genes)

    print(f"\n -> Iniciando processamento de {total_genes} genes encontrados a partir dos dados do dbSNP...")

    resumo_geral = []

    for idx, gene in enumerate(todos_os_genes, 1):
        if str(gene).strip().lower() in ['nan', '', 'none', 'null']:
            continue

        df_gene = df_merged[df_merged[col_gene] == gene].copy()
        if df_gene.empty:
            continue

        gene_nome = str(gene).upper()
        nome_limpo = "".join(c for c in gene_nome if c.isalnum() or c in ('_', '-'))

        caminho_tsv = os.path.join(PASTA_RELATORIOS_GENES, f"tabela_{nome_limpo}.tsv")
        df_gene.to_csv(caminho_tsv, sep="\t", index=False)

        scores = df_gene[col_score_am].dropna() if col_score_am else pd.Series(dtype=float)
        score_medio = scores.mean() if not scores.empty else np.nan
        score_max = scores.max() if not scores.empty else np.nan
        score_min = scores.min() if not scores.empty else np.nan

        cobertura_pct = round(100 * len(scores) / len(df_gene), 1) if len(df_gene) > 0 else 0.0

        resumo_geral.append({
            "Gene": gene_nome,
            "Total_Variantes": len(df_gene),
            "Variantes_Com_Score_AM": len(scores),
            "Cobertura_AM_%": cobertura_pct,
            "Score_Medio_AM": round(score_medio, 4) if not np.isnan(score_medio) else "N/A",
            "Score_Max_AM": round(score_max, 4) if not np.isnan(score_max) else "N/A",
            "Score_Min_AM": round(score_min, 4) if not np.isnan(score_min) else "N/A",
            "Variantes_Patogenicas": (scores >= 0.56).sum() if not scores.empty else 0,
            "Variantes_Benignas": (scores <= 0.34).sum() if not scores.empty else 0,
        })

        fig = plt.figure(figsize=(16, 10), dpi=300)
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])
        fig.suptitle(f"Análise dbSNP / AlphaMissense — Gene {gene_nome} ({idx}/{total_genes})", fontsize=16, fontweight="bold", y=0.98)

        df_valid = df_gene.dropna(subset=[col_score_am]).copy() if col_score_am else pd.DataFrame()

        # Gráfico 1: Posição vs Score
        ax1 = fig.add_subplot(gs[0, :])
        if not df_valid.empty:
            eixo_x = df_valid[col_pos] if (col_pos and df_valid[col_pos].dropna().count() > 0) else np.arange(len(df_valid))

            tem_hue_valido = (
                col_class_dbsnp is not None
                and col_class_dbsnp in df_valid.columns
                and df_valid[col_class_dbsnp].notna().any()
            )
            hue_var = df_valid[col_class_dbsnp] if tem_hue_valido else None

            plot_kwargs = {
                "data": df_valid,
                "x": eixo_x,
                "y": col_score_am,
                "s": 100,
                "alpha": 0.85,
                "ax": ax1
            }
            if hue_var is not None:
                plot_kwargs["hue"] = hue_var
                plot_kwargs["style"] = hue_var
                plot_kwargs["palette"] = "tab10"

            sns.scatterplot(**plot_kwargs)

            if col_rs and col_rs in df_valid.columns:
                top_variants = df_valid.nlargest(3, col_score_am)
                for _, row in top_variants.iterrows():
                    pos_x = row[col_pos] if (col_pos and not np.isnan(row[col_pos])) else row.name
                    rs_val = str(row[col_rs]).replace('rs', '').strip()
                    ax1.annotate(
                        f"rs{rs_val}\n({row[col_score_am]:.2f})",
                        (pos_x, row[col_score_am]),
                        textcoords="offset points", xytext=(0, 8),
                        ha='center', fontsize=8, weight='bold',
                        bbox=dict(boxstyle="round,pad=0.2", fc="yellow", alpha=0.6)
                    )

            ax1.axhline(0.56, color="#d9534f", linestyle="--", linewidth=1.5, label="Corte Patogênico (>= 0.56)")
            ax1.axhline(0.34, color="#5cb85c", linestyle="--", linewidth=1.5, label="Corte Benigno (<= 0.34)")
            ax1.set_title("Pontuação de Patogenicidade por Posição na Sequência", fontweight="bold")
            ax1.set_ylabel("AlphaMissense Score")
            ax1.set_xlabel("Posição do Aminoácido")
            ax1.set_ylim(-0.05, 1.08)
            ax1.legend(loc="upper right", frameon=True, facecolor="white")
        else:
            ax1.text(0.5, 0.5, "Sem dados numéricos de pontuação AlphaMissense", ha="center", va="center", fontsize=12)

        # Gráfico 2: Histograma
        ax2 = fig.add_subplot(gs[1, 0])
        if not df_valid.empty:
            sns.histplot(df_valid[col_score_am], kde=True, bins=15, color="#337ab7", ax=ax2)
            if not np.isnan(score_medio):
                ax2.axvline(score_medio, color="black", linestyle=":", label=f"Média: {score_medio:.2f}")
                ax2.legend()
            ax2.set_title("Distribuição Numérica dos Scores", fontweight="bold")
            ax2.set_xlabel("AlphaMissense Score")
            ax2.set_ylabel("Frequência de Variantes")
        else:
            ax2.text(0.5, 0.5, "Sem dados suficientes para histograma", ha="center", va="center", fontsize=12)

        # Gráfico 3: Boxplot por classificação
        ax3 = fig.add_subplot(gs[1, 1])
        if not df_valid.empty and col_class_dbsnp and col_class_dbsnp in df_valid.columns and df_valid[col_class_dbsnp].dropna().count() > 0:
            df_box = df_valid.dropna(subset=[col_class_dbsnp])
            sns.boxplot(
                data=df_box,
                x=col_class_dbsnp,
                y=col_score_am,
                hue=col_class_dbsnp,
                palette="Set2",
                legend=False,
                ax=ax3
            )
            sns.stripplot(
                data=df_box,
                x=col_class_dbsnp,
                y=col_score_am,
                color="black",
                alpha=0.5,
                jitter=0.2,
                ax=ax3
            )
            ax3.set_title("Validação: Scores vs Classificação dbSNP / ClinVar", fontweight="bold")
            ax3.set_ylabel("AlphaMissense Score")
            ax3.set_xlabel("Classificação dbSNP")
            ax3.tick_params(axis='x', rotation=25)
        else:
            ax3.text(0.5, 0.5, "Sem cruzamento com classificações dbSNP", ha="center", va="center", fontsize=12)

        plt.tight_layout()
        caminho_figura = os.path.join(PASTA_FIGURAS_GENES, f"figura_gene_{nome_limpo}.png")
        plt.savefig(caminho_figura, dpi=300, bbox_inches="tight")
        plt.close()

    df_resumo = pd.DataFrame(resumo_geral)
    df_resumo.to_csv(os.path.join(PASTA_RESULTADOS, "resumo_estatistico_genes.tsv"), sep="\t", index=False)

    caminho_unificado = os.path.join(PASTA_RESULTADOS, "dbsnp_alphamissense_unificado.tsv")
    df_merged.to_csv(caminho_unificado, sep="\t", index=False)

    print(f"\n[SUCESSO] Processamento concluído!")
    print(f" -> Tabelas salvas em: '{PASTA_RELATORIOS_GENES}/'")
    print(f" -> Gráficos salvos em: '{PASTA_FIGURAS_GENES}/'")
    print(f" -> Resumo consolidado: '{PASTA_RESULTADOS}/resumo_estatistico_genes.tsv'")
    print(f" -> Tabela unificada (usada pela interface): '{caminho_unificado}'")


if __name__ == "__main__":
    df_dbsnp = carregar_dados_dbsnp()
    if df_dbsnp is not None:
        gerar_figuras_e_relatorios(df_dbsnp)