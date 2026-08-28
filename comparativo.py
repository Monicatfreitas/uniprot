from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
import re
import urllib.parse
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Otimização de Sessão HTTP com pool de conexões e retentativas automáticas
session = requests.Session()
adapter = HTTPAdapter(
    pool_connections=50,
    pool_maxsize=50,
    max_retries=Retry(
        total=3, backoff_factor=0.2, status_forcelist=[500, 502, 503, 504]
    ),
)
session.mount("http://", adapter)
session.mount("https://", adapter)
session.headers.update({"User-Agent": "Python-Genomics-Pipeline/2.0"})

# Desativa avisos de certificados SSL não verificados
requests.packages.urllib3.disable_warnings()

# Mapeamento de propriedades dos aminoácidos
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

AA_1_TO_3 = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "E": "Glu",
    "Q": "Gln",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
}


def consultar_gene_uniprot_unidade(gene):
    """Consulta individual de um único gene no UniProt."""
    try:
        query = f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true"
        url = f"https://rest.uniprot.org/uniprotkb/search?query={urllib.parse.quote(query)}&fields=accession,protein_name&size=1"
        resp = session.get(url, timeout=10, verify=False)
        if resp.status_code == 200:
            data = resp.json()
            results = data.get("results", [])
            if results:
                uniprot_id = results[0].get("primaryAccession", "-")
                nome_rec = (
                    results[0]
                    .get("proteinDescription", {})
                    .get("recommendedName", {})
                    .get("fullName", {})
                    .get("value", "-")
                )
                return gene, {"uniprot_id": uniprot_id, "nome_proteina": nome_rec}
    except Exception:
        pass
    return gene, {"uniprot_id": "-", "nome_proteina": "-"}


def buscar_dados_uniprot(lista_genes, max_workers=20):
    """Consulta a API REST do UniProt em paralelo com Múltiplas Threads."""
    mapa_uniprot = {}
    genes_unicos = [
        str(g).strip()
        for g in set(lista_genes)
        if pd.notna(g) and str(g).strip() not in ["", "nan", "NA", "-"]
    ]

    print(
        f"Buscando IDs e nomes no UniProt para {len(genes_unicos)} gene(s) em paralelo..."
    )

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(consultar_gene_uniprot_unidade, g): g
            for g in genes_unicos
        }
        for future in as_completed(futures):
            gene, res = future.result()
            mapa_uniprot[gene] = res

    return mapa_uniprot


def consultar_ncbi_rsid_unidade(rsid):
    """Consulta individual de um rsID na API do dbSNP NCBI mantendo buscas regex originais."""
    clean_id = rsid.lower().replace("rs", "").strip()
    if not clean_id.isdigit():
        return rsid, "-"

    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{clean_id}"
    hgvs_encontrado = "-"

    try:
        resp = session.get(url, timeout=5, verify=False)
        if resp.status_code == 200:
            data = resp.json()
            primary = data.get("primary_snapshot_data", {})

            for allele in primary.get("allele_annotations", []):
                for assembly in allele.get("assembly_annotation", []):
                    for hgvs in assembly.get("hgvs_genomic", []):
                        if "p." in hgvs:
                            hgvs_encontrado = hgvs
                            break
                    if hgvs_encontrado != "-":
                        break
                if hgvs_encontrado != "-":
                    break

            if hgvs_encontrado == "-":
                json_str = json.dumps(primary)
                matches_3 = re.findall(
                    r"p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}", json_str
                )
                if matches_3:
                    hgvs_encontrado = matches_3[0]
                else:
                    matches_1 = re.findall(r"p\.[A-Z]\d+[A-Z]", json_str)
                    if matches_1:
                        hgvs_encontrado = matches_1[0]
    except Exception:
        pass

    return rsid, hgvs_encontrado


def buscar_aminoácidos_ncbi_batch(lista_rsids, max_workers=30):
    """Consulta paralelizada de alta velocidade para centenas de milhares de rsIDs."""
    rsids_unicos = list(
        set([
            str(r).strip()
            for r in lista_rsids
            if pd.notna(r) and str(r).startswith("rs")
        ])
    )
    cache_rsid = {}
    total = len(rsids_unicos)

    print(
        f"Iniciando busca paralela otimizada ({max_workers} threads) para {total} rsIDs únicos..."
    )

    processed_count = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(consultar_ncbi_rsid_unidade, rsid): rsid
            for rsid in rsids_unicos
        }

        for future in as_completed(futures):
            processed_count += 1
            rsid, hgvs = future.result()
            cache_rsid[rsid] = hgvs

            # Feedback periódico idêntico ao console original
            if processed_count % 2000 == 0 or processed_count == total:
                print(
                    f"  [{processed_count}/{total}] Processando {rsid}...",
                    flush=True,
                )

    print("\nConsultas ao NCBI concluídas!")
    return cache_rsid


def extrair_detalhes_troca_aminoacido(var_str):
    """Desmembra e classifica a troca de aminoácido."""
    if var_str == "-" or not var_str or pd.isna(var_str):
        return {
            "HGVS_Proteina": "-",
            "AA_Ref": "-",
            "AA_Alt": "-",
            "Posicao": "-",
            "Troca_Formatada": "-",
            "Tipo_Modificacao": "Não Informado",
        }

    match_3 = re.search(
        r"(?:p\.)?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})",
        str(var_str),
        re.IGNORECASE,
    )
    if match_3:
        aa_ref = match_3.group(1).capitalize()
        posicao = match_3.group(2)
        aa_alt = match_3.group(3).capitalize()

        prop_ref = AA_PROPS.get(aa_ref, "Desconhecido")
        prop_alt = AA_PROPS.get(aa_alt, "Desconhecido")

        tipo = (
            f"Conservativa ({prop_ref})"
            if prop_ref == prop_alt
            else f"Não-Conservativa ({prop_ref} -> {prop_alt})"
        )
        return {
            "HGVS_Proteina": var_str,
            "AA_Ref": aa_ref,
            "AA_Alt": aa_alt,
            "Posicao": posicao,
            "Troca_Formatada": f"{aa_ref}{posicao}{aa_alt}",
            "Tipo_Modificacao": tipo,
        }

    match_1 = re.search(r"(?:p\.)?([A-Z])(\d+)([A-Z])", str(var_str))
    if match_1:
        a1, posicao, a2 = match_1.group(1), match_1.group(2), match_1.group(3)
        aa_ref = AA_1_TO_3.get(a1, a1)
        aa_alt = AA_1_TO_3.get(a2, a2)

        prop_ref = AA_PROPS.get(aa_ref, "Desconhecido")
        prop_alt = AA_PROPS.get(aa_alt, "Desconhecido")

        tipo = (
            f"Conservativa ({prop_ref})"
            if prop_ref == prop_alt
            else f"Não-Conservativa ({prop_ref} -> {prop_alt})"
        )
        return {
            "HGVS_Proteina": var_str,
            "AA_Ref": aa_ref,
            "AA_Alt": aa_alt,
            "Posicao": posicao,
            "Troca_Formatada": f"{aa_ref}{posicao}{aa_alt} ({a1}>{a2})",
            "Tipo_Modificacao": tipo,
        }

    return {
        "HGVS_Proteina": var_str,
        "AA_Ref": "-",
        "AA_Alt": "-",
        "Posicao": "-",
        "Troca_Formatada": "-",
        "Tipo_Modificacao": "Substituição Missense",
    }


def gerar_tabela_completa_proteina_tsv():
    print("Iniciando processo de consolidação e busca online...")
    pasta_resultados = "resultados"
    os.makedirs(pasta_resultados, exist_ok=True)

    try:
        df_dbsnp = pd.read_csv("dbsnp_missense_resultados.tsv", sep="\t")
        df_am = pd.read_csv("alphamissense_resultados.tsv", sep="\t")
    except FileNotFoundError as e:
        print(f"[ERRO] Arquivo de entrada não encontrado: {e}")
        return

    # Trata espaços e tipos
    for df in [df_dbsnp, df_am]:
        if "rsID_dbSNP" in df.columns:
            df["rsID_dbSNP"] = df["rsID_dbSNP"].astype(str).str.strip()
        if "Gene" in df.columns:
            df["Gene"] = df["Gene"].astype(str).str.strip()

    # Outer Merge Preservando TODAS as colunas das ferramentas (scores, classificações, etc.)
    df_merged = pd.merge(
        df_dbsnp,
        df_am,
        on=["Gene", "rsID_dbSNP"],
        how="outer",
        suffixes=("_dbSNP", "_AlphaMissense"),
    )

    # 1. Mapear ID UniProt e Nome Oficial
    mapa_uniprot = buscar_dados_uniprot(
        df_merged["Gene"].dropna().unique(), max_workers=20
    )
    df_merged["UniProt_ID"] = df_merged["Gene"].map(
        lambda g: mapa_uniprot.get(g, {}).get("uniprot_id", "-")
    )
    df_merged["Nome_Oficial_Proteina"] = df_merged["Gene"].map(
        lambda g: mapa_uniprot.get(g, {}).get("nome_proteina", "-")
    )
    df_merged["Link_UniProt"] = df_merged["UniProt_ID"].apply(
        lambda x: f"https://www.uniprot.org/uniprotkb/{x}" if x != "-" else "-"
    )

    # 2. Reconciliação do HGVS de Proteína (Tenta a partir dos inputs ou consulta NCBI)
    col_prot_input = None
    for c in ["protein_variant", "hgvs_p", "proteina_mudanca", "HGVS_p"]:
        if c in df_merged.columns:
            col_prot_input = c
            break

    mapa_rsid_hgvs = buscar_aminoácidos_ncbi_batch(
        df_merged["rsID_dbSNP"], max_workers=30
    )

    def obter_hgvs_final(row):
        hgvs_ncbi = mapa_rsid_hgvs.get(str(row["rsID_dbSNP"]).strip(), "-")
        if hgvs_ncbi != "-":
            return hgvs_ncbi
        if (
            col_prot_input
            and pd.notna(row[col_prot_input])
            and str(row[col_prot_input]).strip() != ""
        ):
            return row[col_prot_input]
        return "-"

    df_merged["HGVS_Proteina_Consolidado"] = df_merged.apply(
        obter_hgvs_final, axis=1
    )

    # 3. Detalhamento e Classificação Fisico-Química das Trocas
    detalhes = [
        extrair_detalhes_troca_aminoacido(m)
        for m in df_merged["HGVS_Proteina_Consolidado"]
    ]

    df_merged["Proteina_HGVS"] = [d["HGVS_Proteina"] for d in detalhes]
    df_merged["AA_Referencia"] = [d["AA_Ref"] for d in detalhes]
    df_merged["AA_Alterado"] = [d["AA_Alt"] for d in detalhes]
    df_merged["Posicao_Proteina"] = [d["Posicao"] for d in detalhes]
    df_merged["Troca_Aminoacido"] = [d["Troca_Formatada"] for d in detalhes]
    df_merged["Tipo_Modificacao_FisicoQuimica"] = [
        d["Tipo_Modificacao"] for d in detalhes
    ]

    # 4. Análise de Comparação de Presença (dbSNP vs AlphaMissense)
    set_dbsnp = set(df_dbsnp["rsID_dbSNP"].dropna())
    set_am = set(df_am["rsID_dbSNP"].dropna())

    def determinar_fonte(rsid):
        no_dbsnp = rsid in set_dbsnp
        no_am = rsid in set_am
        if no_dbsnp and no_am:
            return "Em Ambos (dbSNP + AlphaMissense)"
        elif no_dbsnp:
            return "Apenas dbSNP"
        elif no_am:
            return "Apenas AlphaMissense"
        return "Desconhecido"

    df_merged["Comparativo_Fontes"] = df_merged["rsID_dbSNP"].apply(
        determinar_fonte
    )
    df_merged["Link_dbSNP"] = df_merged["rsID_dbSNP"].apply(
        lambda x: (
            f"https://www.ncbi.nlm.nih.gov/snp/{str(x).strip()}"
            if str(x).startswith("rs")
            else "-"
        )
    )

    # 5. Organização Inteligente de Colunas
    cols_bloco_identificacao = [
        "Gene",
        "UniProt_ID",
        "Nome_Oficial_Proteina",
        "rsID_dbSNP",
        "Comparativo_Fontes",
    ]

    cols_bloco_mutacao = [
        "Proteina_HGVS",
        "AA_Referencia",
        "AA_Alterado",
        "Posicao_Proteina",
        "Troca_Aminoacido",
        "Tipo_Modificacao_FisicoQuimica",
    ]

    cols_bloco_links = [
        "Link_UniProt",
        "Link_dbSNP",
    ]

    cols_scores_original = [
        c
        for c in df_merged.columns
        if c
        not in cols_bloco_identificacao + cols_bloco_mutacao + cols_bloco_links
        and c not in ["HGVS_Proteina_Consolidado"]
    ]

    ordem_final = (
        cols_bloco_identificacao
        + cols_bloco_mutacao
        + cols_scores_original
        + cols_bloco_links
    )
    df_merged = df_merged[ordem_final]

    output_path = os.path.join(
        pasta_resultados, "tabela_comparativa_detalhada_proteinas.tsv"
    )
    df_merged.to_csv(output_path, sep="\t", index=False)
    print(f"\n[SUCESSO] Tabela unificada gerada em: '{output_path}'!")


if __name__ == "__main__":
    gerar_tabela_completa_proteina_tsv()