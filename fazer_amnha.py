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

def carregar_dados_dbsnp():
    # 1. Carrega a tabela principal
    df_dbsnp = buscar_e_carregar_tsv(ARQUIVOS_DBSNP, obrigatorio=True)
    if df_dbsnp is None:
        return None

    df_dbsnp = tratar_colunas_base(df_dbsnp)

    # 2. Carrega o arquivo do AlphaMissense
    df_am = buscar_e_carregar_tsv(ARQUIVOS_EXTRAS["alphamissense"], obrigatorio=False)

    if df_am is not None:
        df_am = tratar_colunas_base(df_am)

        col_score_am = next((c for c in df_am.columns if any(k in c.lower() for k in ["pathogenicity", "score", "am_score", "alphamissense"])), None)

        if col_score_am:
            df_am[col_score_am] = limpar_e_converter_numerico(df_am[col_score_am])

        col_gene = next((c for c in df_dbsnp.columns if "gene" in c.lower()), "Gene")

        usa_rs_key = 'rs_key' in df_dbsnp.columns and 'rs_key' in df_am.columns
        chave_merge = [col_gene, 'rs_key'] if usa_rs_key else [col_gene]

        # CORREÇÃO: a tabela dbSNP já vem com colunas placeholder (ex: AlphaMissense_Score
        # vazia) que têm o MESMO NOME das colunas reais da tabela AlphaMissense. Sem isso,
        # o merge cria uma coluna "_AM" com os dados de verdade escondida atrás da vazia,
        # e o restante do script sempre lê a vazia.
        colunas_para_substituir = [c for c in df_am.columns if c in df_dbsnp.columns and c not in chave_merge]
        if colunas_para_substituir:
            print(f"[INFO] Substituindo colunas placeholder do dbSNP pelos dados reais do AlphaMissense: {colunas_para_substituir}")
            df_dbsnp = df_dbsnp.drop(columns=colunas_para_substituir)

        if usa_rs_key:
            df_merged = pd.merge(df_dbsnp, df_am, on=[col_gene, 'rs_key'], how='left')
        else:
            df_merged = pd.merge(df_dbsnp, df_am, on=col_gene, how='left')

        return df_merged

    return df_dbsnp

def gerar_figuras_e_relatorios(df_merged: pd.DataFrame):
    os.makedirs(PASTA_FIGURAS_GENES, exist_ok=True)
    os.makedirs(PASTA_RELATORIOS_GENES, exist_ok=True)

    col_gene = next((c for c in df_merged.columns if "gene" in c.lower()), None)
    if not col_gene:
        print("[ERRO] Coluna referente a 'Gene' não foi encontrada.")
        return

    col_score_am = next((c for c in df_merged.columns if any(k in c.lower() for k in ["pathogenicity", "score", "am_score", "alphamissense", "am_pathogenicity"])), None)
    col_class_dbsnp = next((c for c in df_merged.columns if any(k in c.lower() for k in ["classificacao", "clinical", "significance", "clinvar"])), None)
    col_pos = next((c for c in df_merged.columns if any(k in c.lower() for k in ["pos", "protein_start", "position", "posicao"])), None)
    col_rs = next((c for c in df_merged.columns if "rs" in c.lower() or "dbsnp" in c.lower()), None)

    if col_score_am:
        df_merged[col_score_am] = limpar_e_converter_numerico(df_merged[col_score_am])
        print(f"[INFO] Usando coluna de score: {col_score_am!r} — {df_merged[col_score_am].notna().sum()} valores válidos de {len(df_merged)}")
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

        resumo_geral.append({
            "Gene": gene_nome,
            "Total_Variantes": len(df_gene),
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
            hue_var = df_valid[col_class_dbsnp] if (col_class_dbsnp and col_class_dbsnp in df_valid.columns) else None

            sns.scatterplot(
                data=df_valid, x=eixo_x, y=col_score_am, hue=hue_var,
                style=hue_var, s=100, alpha=0.85, palette="tab10", ax=ax1
            )

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
            sns.boxplot(data=df_valid.dropna(subset=[col_class_dbsnp]), x=col_class_dbsnp, y=col_score_am, palette="Set2", ax=ax3)
            sns.stripplot(data=df_valid.dropna(subset=[col_class_dbsnp]), x=col_class_dbsnp, y=col_score_am, color="black", alpha=0.5, jitter=0.2, ax=ax3)
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

    print(f"\n[SUCESSO] Processamento concluído!")
    print(f" -> Tabelas salvas em: '{PASTA_RELATORIOS_GENES}/'")
    print(f" -> Gráficos salvos em: '{PASTA_FIGURAS_GENES}/'")
    print(f" -> Resumo consolidado: '{PASTA_RESULTADOS}/resumo_estatistico_genes.tsv'")

if __name__ == "__main__":
    df_dbsnp = carregar_dados_dbsnp()
    if df_dbsnp is not None:
        gerar_figuras_e_relatorios(df_dbsnp)
