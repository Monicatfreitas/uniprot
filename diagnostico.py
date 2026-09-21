"""
Versão de DIAGNÓSTICO do seu pipeline. Roda exatamente a mesma lógica de carregamento
e merge, mas imprime no console o que está acontecendo em cada etapa: quantas linhas
sobraram, quais colunas foram detectadas, e quantos valores de score não são nulos.

Rode este arquivo e me mande a saída do console — com isso eu identifico exatamente
onde os dados estão sendo perdidos (arquivo não encontrado, merge sem correspondência,
ou coluna de score não detectada).
"""

import os
from typing import Optional
import numpy as np
import pandas as pd

PASTA_RESULTADOS = "resultados"

ARQUIVOS_DBSNP = [
    "comparacao_dbsnp_alphamissense.tsv",
    "dbsnp_missense_resultados.tsv",
    "dbsnp_alphamissense_unificado.tsv",
    "tabela_comparativa_completa.tsv",
]

ARQUIVOS_EXTRAS = {
    "alphamissense": ["alphamissense_resultados.tsv"],
    "uniprot": ["resultados_uniprot.tsv"],
}


def limpar_e_converter_numerico(serie: pd.Series) -> pd.Series:
    s = serie.astype(str).str.replace(',', '.').str.strip()
    s = s.replace(['nan', 'none', 'n/a', '-', 'null', '', '<na>', 'None'], np.nan)
    return pd.to_numeric(s, errors="coerce")


def buscar_e_carregar_tsv(lista_nomes: list, obrigatorio: bool = True) -> Optional[pd.DataFrame]:
    for nome in lista_nomes:
        for caminho in [nome, os.path.join(PASTA_RESULTADOS, nome)]:
            if os.path.exists(caminho):
                try:
                    df = pd.read_csv(caminho, sep="\t", engine="python", on_bad_lines="skip")
                    df.columns = df.columns.str.strip()
                    print(f"[OK] Carregado: {caminho}  |  shape={df.shape}")
                    print(f"     colunas: {list(df.columns)}")
                    return df
                except Exception as e:
                    print(f"[AVISO] Falha ao ler {caminho}: {e}")
    if obrigatorio:
        print(f"[ERRO] Nenhum arquivo encontrado para: {lista_nomes}")
    return None


def tratar_colunas_base(df: pd.DataFrame, rotulo: str) -> pd.DataFrame:
    if "Gene_Original" in df.columns and "Gene" not in df.columns:
        df = df.rename(columns={"Gene_Original": "Gene"})

    for col in df.columns:
        if col.lower() in ["gene", "gene_name", "symbol"]:
            df[col] = df[col].astype(str).str.strip().str.upper()
        elif "rs" in col.lower() or "dbsnp" in col.lower():
            df["rs_key"] = df[col].astype(str).str.replace("rs", "", case=False).str.strip()

    if "rs_key" in df.columns:
        print(f"     [{rotulo}] exemplo de rs_key: {df['rs_key'].dropna().unique()[:5].tolist()}")
    else:
        print(f"     [{rotulo}] NENHUMA coluna de rsID detectada!")

    return df


def diagnosticar():
    print("=" * 70)
    print("ETAPA 1: carregando tabela principal dbSNP")
    print("=" * 70)
    df_dbsnp = buscar_e_carregar_tsv(ARQUIVOS_DBSNP, obrigatorio=True)
    if df_dbsnp is None:
        print(">> PARANDO: nenhum arquivo dbSNP foi encontrado. Confira o nome/local do arquivo.")
        return
    df_dbsnp = tratar_colunas_base(df_dbsnp, "dbSNP")

    print("\n" + "=" * 70)
    print("ETAPA 2: carregando tabela AlphaMissense")
    print("=" * 70)
    df_am = buscar_e_carregar_tsv(ARQUIVOS_EXTRAS["alphamissense"], obrigatorio=False)
    if df_am is None:
        print(">> PARANDO: alphamissense_resultados.tsv não foi encontrado.")
        return
    df_am = tratar_colunas_base(df_am, "AlphaMissense")

    col_score_am = next(
        (c for c in df_am.columns
         if any(k in c.lower() for k in ["pathogenicity", "score", "am_score", "alphamissense"])),
        None,
    )
    print(f"\n     Coluna de score detectada: {col_score_am!r}")
    if col_score_am is None:
        print("     >> PROBLEMA ENCONTRADO: nenhuma coluna bateu com as palavras-chave.")
        print("     >> Copie aqui a lista de colunas acima e eu ajusto a busca.")
        return

    df_am[col_score_am] = limpar_e_converter_numerico(df_am[col_score_am])
    print(f"     Valores não-nulos em {col_score_am!r} ANTES do merge: "
          f"{df_am[col_score_am].notna().sum()} de {len(df_am)}")

    print("\n" + "=" * 70)
    print("ETAPA 3: merge dbSNP + AlphaMissense")
    print("=" * 70)
    col_gene = next((c for c in df_dbsnp.columns if "gene" in c.lower()), "Gene")
    print(f"     Coluna de gene usada no merge: {col_gene!r}")

    usa_rs_key = "rs_key" in df_dbsnp.columns and "rs_key" in df_am.columns
    print(f"     Merge vai usar rs_key além do gene? {usa_rs_key}")

    if usa_rs_key:
        df_merged = pd.merge(df_dbsnp, df_am, on=[col_gene, "rs_key"], how="left", suffixes=("", "_AM"))
    else:
        df_merged = pd.merge(df_dbsnp, df_am, on=col_gene, how="left", suffixes=("", "_AM"))

    print(f"     Shape depois do merge: {df_merged.shape}")

    col_score_final = next(
        (c for c in df_merged.columns
         if any(k in c.lower() for k in ["pathogenicity", "score", "am_score", "alphamissense", "am_pathogenicity"])),
        None,
    )
    print(f"     Coluna de score no df final: {col_score_final!r}")

    if col_score_final:
        n_validos = pd.to_numeric(
            df_merged[col_score_final].astype(str).str.replace(",", ".").str.strip(),
            errors="coerce",
        ).notna().sum()
        print(f"     Valores não-nulos em {col_score_final!r} DEPOIS do merge: {n_validos} de {len(df_merged)}")
        if n_validos == 0:
            print("\n     >> PROBLEMA ENCONTRADO: o merge não casou nenhuma linha.")
            print("     >> Provável causa: rs_key ou nome do gene com formatos diferentes entre as tabelas.")
            print(f"     >> Exemplos de {col_gene} no dbSNP: {df_dbsnp[col_gene].dropna().unique()[:5].tolist()}")
            gene_am_col = next((c for c in df_am.columns if 'gene' in c.lower()), None)
            if gene_am_col:
                print(f"     >> Exemplos de gene no AlphaMissense: {df_am[gene_am_col].dropna().unique()[:5].tolist()}")
        else:
            print("\n     >> Merge OK, os dados deveriam aparecer nos gráficos.")
    else:
        print("     >> PROBLEMA ENCONTRADO: coluna de score não existe no df final (merge pode ter falhado).")


if __name__ == "__main__":
    diagnosticar()
