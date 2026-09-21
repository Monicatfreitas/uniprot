import os
import re
import pandas as pd

# ---------------------------------------------------------------- Caminhos ---

DIR_BASE = os.path.dirname(os.path.abspath(__file__))
DIR_RESULTADOS = os.path.join(DIR_BASE, "resultados")

CANDIDATOS_ENTRADA = [
    os.path.join(DIR_RESULTADOS, "dbsnp_alphamissense_unificado.tsv"),
    os.path.join(DIR_BASE, "dbsnp_alphamissense_unificado.tsv"),
    os.path.join(DIR_BASE, "dbsnp_missense_resultados.tsv"),
]

CAMINHO_ALPHAMISSENSE = os.path.join(DIR_BASE, "AlphaMissense_hg38.tsv.gz")
CAMINHO_SAIDA = os.path.join(DIR_RESULTADOS, "dbsnp_alphamissense_anotado.tsv")

VAZIOS = {"", "none", "nan", "n/a", "na", "null", "-", "<na>"}

# ----------------------------------------------------------------- Helpers ---


def normalizar_rsid(valor):
    if valor is None or (isinstance(valor, float) and pd.isna(valor)):
        return None
    texto = str(valor).strip().lower()
    if texto in VAZIOS:
        return None
    match = re.fullmatch(r"(?:rs)?(\d+)", texto)
    return f"rs{match.group(1)}" if match else None


def localizar_coluna(df, opcoes):
    """Localiza a coluna ignorando maiúsculas, minúsculas e caracteres especiais."""
    for col in df.columns:
        col_limpa = col.strip().lower()
        for opcao in opcoes:
            if col_limpa == opcao.lower():
                return col
    return None


# ---------------------------------------------------------------- Pipeline ---


def preencher_posicoes():
    entrada = next((c for c in CANDIDATOS_ENTRADA if os.path.exists(c)), None)
    if not entrada:
        print("Erro: Nenhum arquivo de entrada encontrado.")
        return

    if not os.path.exists(CAMINHO_ALPHAMISSENSE):
        print(f"Erro: O arquivo '{CAMINHO_ALPHAMISSENSE}' não foi encontrado.")
        print(
            "Baixe-o executando: wget"
            " https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
        )
        return

    print("Carregando tabela de entrada...")
    df_entrada = pd.read_csv(entrada, sep="\t", low_memory=False)

    col_rs = localizar_coluna(
        df_entrada,
        ["rs_key", "rsid", "rs_id", "dbsnp", "dbsnp_id", "variation_id", "rs"],
    )
    col_cdna = localizar_coluna(
        df_entrada,
        ["Posição no Gene / cDNA", "Posicao_Gene", "posicao_gene", "cdna"],
    )
    col_genomica = localizar_coluna(
        df_entrada,
        [
            "Posição Genômica / Chr",
            "Posicao_Genomica",
            "posicao_genomica",
            "genomica",
        ],
    )

    if not col_rs:
        print("Erro: Coluna de rsID não encontrada no arquivo de entrada.")
        return

    if not col_cdna:
        col_cdna = "Posição no Gene / cDNA"
        df_entrada[col_cdna] = None
    if not col_genomica:
        col_genomica = "Posição Genômica / Chr"
        df_entrada[col_genomica] = None

    df_entrada["_rsid"] = df_entrada[col_rs].map(normalizar_rsid)

    print("Lendo posições do AlphaMissense local...")
    df_am = pd.read_csv(
        CAMINHO_ALPHAMISSENSE,
        sep="\t",
        comment="#",
        compression="gzip",
        low_memory=False,
    )

    # Remove cerquilha (#) do cabeçalho se houver
    df_am.columns = df_am.columns.str.replace("^#", "", regex=True).str.strip()

    col_am_rs = localizar_coluna(df_am, ["dbsnp_id", "rsid", "rs_id"])
    col_am_chrom = localizar_coluna(df_am, ["CHROM", "chrom", "chr"])
    col_am_pos = localizar_coluna(df_am, ["POS", "pos", "position"])
    col_am_tx = localizar_coluna(
        df_am, ["transcript_id", "transcript", "tx_id"]
    )

    if not col_am_rs or not col_am_chrom or not col_am_pos:
        print(
            "Erro: Colunas primárias não encontradas no AlphaMissense."
            f" Colunas disponíveis: {list(df_am.columns)}"
        )
        return

    df_am = df_am.dropna(subset=[col_am_rs]).copy()
    df_am["_rsid"] = df_am[col_am_rs].map(normalizar_rsid)
    df_am = df_am.dropna(subset=["_rsid"])

    # Gera os dados genômicos (chrX:POS) e obtém transcrição
    df_am["val_genomica"] = (
        "chr"
        + df_am[col_am_chrom].astype(str)
        + ":"
        + df_am[col_am_pos].astype(str)
    )
    df_am["val_cdna"] = df_am[col_am_tx].astype(str) if col_am_tx else None

    # Mapeamentos em memória
    mapa_genomico = (
        df_am.drop_duplicates(subset=["_rsid"])
        .set_index("_rsid")["val_genomica"]
        .to_dict()
    )
    mapa_cdna = (
        df_am.drop_duplicates(subset=["_rsid"])
        .set_index("_rsid")["val_cdna"]
        .to_dict()
        if col_am_tx
        else {}
    )

    print("Preenchendo apenas os valores nulos/None das duas colunas...")
    novos_genomica = df_entrada["_rsid"].map(mapa_genomico)
    novos_cdna = df_entrada["_rsid"].map(mapa_cdna)

    mask_genomica_vazio = df_entrada[col_genomica].isna() | df_entrada[
        col_genomica
    ].astype(str).str.strip().str.lower().isin(VAZIOS)
    mask_cdna_vazio = df_entrada[col_cdna].isna() | df_entrada[
        col_cdna
    ].astype(str).str.strip().str.lower().isin(VAZIOS)

    df_entrada.loc[mask_genomica_vazio, col_genomica] = novos_genomica[
        mask_genomica_vazio
    ]
    df_entrada.loc[mask_cdna_vazio, col_cdna] = novos_cdna[mask_cdna_vazio]

    df_entrada.drop(columns=["_rsid"], inplace=True)

    os.makedirs(DIR_RESULTADOS, exist_ok=True)
    df_entrada.to_csv(CAMINHO_SAIDA, sep="\t", index=False)

    print(
        f"\nConcluído com sucesso! {len(df_entrada):,} linhas atualizadas e"
        f" salvas em:\n{CAMINHO_SAIDA}"
    )


if __name__ == "__main__":
    preencher_posicoes()