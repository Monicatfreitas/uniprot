import os
import re
import json
import urllib.request
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import networkx as nx

# Configuração da página
st.set_page_config(page_title="Consolidado por Gene", layout="wide")

URL_BASE_ALPHAFOLD = "https://alphafold.ebi.ac.uk/entry/"

MAPA_TRADUCAO_CLINICA = {
    'uncertain significance': 'Sig. Incerta',
    'likely pathogenic': 'Prov. Patogênica',
    'pathogenic': 'Patogênica',
    'likely benign': 'Prov. Benigna',
    'benign': 'Benigna',
    'conflicting interpretations of pathogenicity': 'Interpr. Conflitantes',
    'not provided': 'Não Informado'
}

# Paleta única usada tanto no CSS quanto nos gráficos Plotly
COR_FUNDO = "#020617"
COR_SUPERFICIE = "#0F172A"
COR_BORDA = "#1E293B"
COR_ACENTO = "#38BDF8"
COR_TEXTO = "#E2E8F0"
COR_SECUNDARIA = "#64748B"
COR_ALERTA = "#F87171"
COR_OK = "#34D399"

st.markdown(f"""
<style>
    /* --- IMPORTS & TIPOGRAFIA DE ALTA PRECISÃO --- */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@500;600&display=swap');

    :root {{
        --bg: {COR_FUNDO};
        --surface: {COR_SUPERFICIE};
        --border: {COR_BORDA};
        --accent: {COR_ACENTO};
        --text: {COR_TEXTO};
        --muted: {COR_SECUNDARIA};
    }}

    html, body, [class*="css"] {{
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif !important;
        font-size: 15px !important;
        color: var(--text) !important;
    }}

    .stApp {{
        background-color: var(--bg) !important;
    }}

    /* --- ESTRUTURA E ESPAÇAMENTO DO CANVAS --- */
    .block-container {{
        padding-top: 3rem !important;
        padding-bottom: 3rem !important;
        max-width: 1400px !important;
    }}

    /* --- CABEÇALHO DO APP --- */
    .app-header {{
        display: flex;
        justify-content: space-between;
        align-items: flex-end;
        flex-wrap: wrap;
        gap: 12px;
        margin-bottom: 1.5rem !important;
        padding-bottom: 1.1rem !important;
        border-bottom: 1px solid var(--border) !important;
    }}
    .app-name {{
        font-size: 1.6rem !important;
        font-weight: 700 !important;
        color: #F8FAFC !important;
        letter-spacing: -0.01em !important;
        margin: 0 0 4px 0 !important;
        line-height: 1.2 !important;
    }}
    .app-desc {{
        font-size: 0.92rem !important;
        color: var(--muted) !important;
        margin: 0 !important;
        max-width: 62ch;
        line-height: 1.5 !important;
    }}
    .app-sources {{
        font-size: 0.78rem !important;
        color: #475569 !important;
        text-align: right;
        white-space: nowrap;
        padding-bottom: 3px;
    }}

    /* --- BUSCA DE GENE --- */
    .gene-search-label {{
        font-size: 0.85rem !important;
        font-weight: 600 !important;
        color: #CBD5E1 !important;
        margin: 0 0 0.5rem 0 !important;
    }}
    div[data-testid="stSelectbox"] {{
        max-width: 480px;
    }}

    /* --- CARTÃO DE IDENTIFICAÇÃO DO GENE --- */
    .header-card {{
        background: linear-gradient(135deg, #0F172A 0%, #0B1220 100%);
        border: 1px solid var(--border);
        border-left: 3px solid var(--accent);
        border-radius: 12px;
        padding: 24px 28px;
        margin: 0.5rem 0 1.75rem 0;
    }}
    .gene-label {{
        font-size: 0.78rem;
        font-weight: 600;
        color: var(--muted);
        letter-spacing: 0.04em;
        margin-bottom: 6px;
    }}
    .gene-title-container {{
        display: flex;
        align-items: baseline;
        flex-wrap: wrap;
        gap: 14px;
    }}
    .gene-title {{
        font-family: 'JetBrains Mono', monospace !important;
        font-size: 2.2rem !important;
        font-weight: 600 !important;
        color: #F8FAFC !important;
        letter-spacing: -0.02em !important;
        margin: 0 !important;
        padding: 0 !important;
        line-height: 1.1 !important;
    }}
    .uniprot-badge {{
        font-family: 'JetBrains Mono', monospace;
        font-size: 0.8rem;
        font-weight: 500;
        color: var(--accent);
        background: rgba(56, 189, 248, 0.1);
        border: 1px solid rgba(56, 189, 248, 0.3);
        border-radius: 999px;
        padding: 4px 12px;
        white-space: nowrap;
    }}
    .protein-subtitle {{
        color: #94A3B8;
        font-size: 0.95rem;
        margin-top: 12px;
        line-height: 1.5;
        max-width: 70ch;
    }}
    .protein-subtitle b {{
        color: #CBD5E1;
        font-weight: 600;
    }}

    /* --- TÍTULOS E CABEÇALHOS DE SEÇÃO --- */
    .section-header {{
        font-size: 1.25rem !important;
        font-weight: 600 !important;
        color: #F8FAFC !important;
        letter-spacing: -0.02em !important;
        margin-top: 2rem !important;
        margin-bottom: 1.2rem !important;
        padding-bottom: 0.6rem !important;
        border-bottom: 1px solid var(--border) !important;
        display: flex !important;
        align-items: center !important;
        gap: 8px !important;
    }}
    .sub-header {{
        font-size: 0.95rem !important;
        font-weight: 600 !important;
        color: #CBD5E1 !important;
        margin: 0.4rem 0 1rem 0 !important;
    }}

    /* --- ABAS DE NAVEGAÇÃO MODERNAS (TABS) --- */
    .stTabs [data-baseweb="tab-list"] {{
        gap: 32px !important;
        border-bottom: 1px solid var(--border) !important;
        background-color: transparent !important;
        padding-bottom: 0px !important;
        margin-bottom: 1.5rem !important;
        flex-wrap: wrap !important;
    }}
    .stTabs [data-baseweb="tab"] {{
        height: 42px !important;
        padding: 0 4px 12px 4px !important;
        background-color: transparent !important;
        border: none !important;
        color: var(--muted) !important;
        font-size: 1rem !important;
        font-weight: 500 !important;
        transition: color 0.2s ease !important;
    }}
    .stTabs [data-baseweb="tab"]:hover {{
        color: #94A3B8 !important;
    }}
    .stTabs [aria-selected="true"] {{
        background-color: transparent !important;
        color: var(--accent) !important;
        font-weight: 600 !important;
        border-bottom: 2px solid var(--accent) !important;
        border-radius: 0 !important;
    }}
    .stTabs [data-baseweb="tab-highlight"] {{
        background-color: transparent !important;
    }}

    /* --- CARDS DE ANOTAÇÕES FUNCIONAIS (GRID REFINADO) --- */
    .annotation-card {{
        background: var(--surface);
        border: 1px solid var(--border);
        border-radius: 10px;
        padding: 20px 22px;
        margin-bottom: 20px;
        transition: border-color 0.2s ease;
    }}
    .annotation-card:hover {{
        border-color: #334155;
    }}
    .annotation-title {{
        font-size: 0.9rem !important;
        letter-spacing: 0.01em;
        color: var(--accent);
        font-weight: 600;
        margin-bottom: 14px;
        padding-bottom: 8px;
        border-bottom: 1px solid var(--border);
    }}
    .annotation-list {{
        list-style-type: none !important;
        padding-left: 0 !important;
        margin: 0 !important;
        max-height: 260px;
        overflow-y: auto;
    }}
    .annotation-item {{
        color: #CBD5E1 !important;
        font-size: 0.92rem !important;
        padding: 6px 0 !important;
        line-height: 1.5 !important;
        display: flex;
        align-items: flex-start;
        gap: 10px;
    }}
    .annotation-item::before {{
        content: "•";
        color: var(--accent);
        font-weight: bold;
        font-size: 1.1rem;
        line-height: 1.3;
    }}
    .annotation-empty {{
        color: var(--muted) !important;
        font-size: 0.9rem !important;
        padding: 6px 0;
    }}

    /* --- SEQUÊNCIA FASTA --- */
    .fasta-box {{
        background: var(--surface);
        border: 1px solid var(--border);
        border-radius: 10px;
        padding: 18px 20px;
        font-family: 'JetBrains Mono', monospace;
        font-size: 0.85rem;
        line-height: 1.9;
        color: #94A3B8;
        overflow-x: auto;
        white-space: pre;
    }}
    .fasta-pos {{
        color: #475569;
        user-select: none;
    }}
    .fasta-hit {{
        background: rgba(248, 113, 113, 0.2);
        color: #FCA5A5;
        font-weight: 600;
        border-radius: 3px;
        padding: 1px 2px;
    }}

    /* --- COMPONENTES NATIVOS DO STREAMLIT --- */
    div[data-testid="stMetric"] {{
        background-color: var(--surface) !important;
        border: 1px solid var(--border) !important;
        border-radius: 10px !important;
        padding: 16px 20px !important;
    }}
    div[data-testid="stMetricLabel"] {{
        color: var(--muted) !important;
        font-size: 0.82rem !important;
        font-weight: 600 !important;
        letter-spacing: 0.01em;
    }}
    div[data-testid="stMetricValue"] {{
        color: #F8FAFC !important;
        font-size: 1.6rem !important;
        font-weight: 700 !important;
        font-family: 'JetBrains Mono', monospace !important;
    }}

    hr, div[data-testid="stDivider"] hr {{
        border-color: var(--border) !important;
    }}

    div[data-testid="stCaptionContainer"] p {{
        color: var(--muted) !important;
        font-size: 0.85rem !important;
    }}

    div[data-baseweb="select"] > div {{
        background-color: var(--surface) !important;
        border-color: var(--border) !important;
        border-radius: 8px !important;
        min-height: 46px !important;
        font-size: 0.98rem !important;
    }}
    div[data-baseweb="select"]:focus-within > div {{
        border-color: var(--accent) !important;
        box-shadow: 0 0 0 1px var(--accent) !important;
    }}

    div[data-testid="stDataFrame"] {{
        border: 1px solid var(--border) !important;
        border-radius: 10px !important;
        overflow: hidden !important;
    }}

    div[data-testid="stAlert"] {{
        background-color: var(--surface) !important;
        border: 1px solid var(--border) !important;
        border-left: 3px solid var(--accent) !important;
        border-radius: 8px !important;
        color: #CBD5E1 !important;
    }}
</style>
""", unsafe_allow_html=True)


def aplicar_tema(fig, altura=None):
    """Aplica o mesmo tema escuro do CSS a qualquer figura Plotly."""
    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="Inter, sans-serif", size=12, color=COR_TEXTO),
        title_font=dict(family="Inter, sans-serif", color="#F8FAFC"),
        hoverlabel=dict(
            bgcolor=COR_SUPERFICIE,
            bordercolor=COR_BORDA,
            font=dict(family="Inter, sans-serif", color=COR_TEXTO, size=12)
        ),
        legend=dict(bgcolor="rgba(0,0,0,0)")
    )
    fig.update_xaxes(gridcolor="rgba(148,163,184,0.10)", zerolinecolor=COR_BORDA,
                     linecolor=COR_BORDA, title_font=dict(size=12, color=COR_SECUNDARIA))
    fig.update_yaxes(gridcolor="rgba(148,163,184,0.10)", zerolinecolor=COR_BORDA,
                     linecolor=COR_BORDA, title_font=dict(size=12, color=COR_SECUNDARIA))
    if altura:
        fig.update_layout(height=altura)
    return fig


def traduzir_significancia(val):
    if pd.isna(val):
        return val
    val_str = str(val).strip().lower()
    for eng, pt in MAPA_TRADUCAO_CLINICA.items():
        val_str = val_str.replace(eng, pt)
    val_str = val_str.title()
    if len(val_str) > 25:
        return "Conflitante / Múltipla"
    return val_str


@st.cache_data(show_spinner="Carregando base de dados...")
def carregar_dados_otimizado():
    dir_base = os.path.dirname(os.path.abspath(__file__))
    caminhos = [
        os.path.join(dir_base, "resultados", "dbsnp_alphamissense_unificado.tsv"),
        os.path.join(dir_base, "comparacao_dbsnp_alphamissense.tsv"),
        os.path.join(dir_base, "tabela_comparativa_completa.tsv"),
        os.path.join(dir_base, "dbsnp_missense_resultados.tsv"),
        os.path.join(dir_base, "alphamissense_resultados.tsv"),
        os.path.join(dir_base, "recuperados_uniprot.tsv")
    ]
    caminho_escolhido = next((c for c in caminhos if os.path.exists(c)), None)
    if not caminho_escolhido:
        return pd.DataFrame()

    df = pd.read_csv(caminho_escolhido, sep="\t", on_bad_lines="skip", low_memory=False)
    df.columns = df.columns.str.strip()

    if "Gene_Original" in df.columns and "Gene" not in df.columns:
        df.rename(columns={"Gene_Original": "Gene"}, inplace=True)

    col_rs = next((c for c in df.columns if 'rs' in c.lower() or 'dbsnp' in c.lower()), None)
    col_change = next((c for c in df.columns if any(k in c.lower() for k in ["proteina_mudanca", "protein_change", "aa_change"])), None)

    if col_rs and col_change:
        df.drop_duplicates(subset=[col_rs, col_change], keep='first', inplace=True)

    col_class = next((c for c in df.columns if any(k in c.lower() for k in ["classificacao", "clinical", "significancia", "clinvar"])), None)
    if col_class:
        df[col_class] = df[col_class].apply(traduzir_significancia)

    return df


@st.cache_data(ttl=86400, show_spinner=False)
def buscar_info_uniprot_api(uniprot_id):
    """ Busca FASTA e anotações funcionais reais via REST UniProt """
    if not uniprot_id or str(uniprot_id).lower() in ["n/a", "none", "nan"]:
        return None, {}

    fasta_seq = None
    anotacoes = {"signal_peptide": [], "signal_range": None, "domains": [], "binding_sites": [], "interactions": []}

    try:
        url_fasta = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
        req = urllib.request.Request(url_fasta, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=5) as resp:
            fasta_seq = resp.read().decode('utf-8')
    except Exception:
        fasta_seq = None

    try:
        url_json = f"https://www.uniprot.org/uniprot/{uniprot_id}.json"
        req = urllib.request.Request(url_json, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=5) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            features = data.get("features", [])
            for f in features:
                ftype = f.get("type")
                desc = f.get("description", "")
                location = f.get("location", {})
                start = location.get("start", {}).get("value")
                end = location.get("end", {}).get("value")
                pos_str = f"({start}-{end})" if start and end else ""

                if ftype == "Signal":
                    anotacoes["signal_peptide"].append(f"Peptídeo Sinal {pos_str} {desc}".strip())
                    if start and end:
                        anotacoes["signal_range"] = (int(start), int(end))
                elif ftype in ["Domain", "Region", "Repeat"]:
                    if desc:
                        anotacoes["domains"].append(f"{desc} {pos_str}".strip())
                elif ftype in ["Binding site", "Active site", "Site", "Metal binding"]:
                    if desc:
                        anotacoes["binding_sites"].append(f"{ftype} {pos_str}: {desc}".strip())

            comments = data.get("comments", [])
            for c in comments:
                if c.get("commentType") == "INTERACTION":
                    for inter in c.get("interactions", []):
                        partner_gene = inter.get("partner1", {}).get("geneName") or inter.get("partner2", {}).get("geneName")
                        partner_acc = inter.get("partner2", {}).get("uniProtKBAccession")
                        if partner_gene:
                            anotacoes["interactions"].append(f"Gene: {partner_gene}")
                        elif partner_acc:
                            anotacoes["interactions"].append(f"UniProt: {partner_acc}")
    except Exception:
        pass

    return fasta_seq, anotacoes


@st.cache_data(ttl=86400, show_spinner="Baixando estrutura 3D...")
def obter_pdb_alphafold(uniprot_id):
    """ Busca o arquivo PDB usando a API oficial da AlphaFold DB """
    if not uniprot_id or str(uniprot_id).lower() in ["n/a", "none", "nan"]:
        return None

    uniprot_clean = str(uniprot_id).strip().upper()

    try:
        url_api = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_clean}"
        req = urllib.request.Request(url_api, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=8) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            if isinstance(data, list) and len(data) > 0:
                pdb_url = data[0].get("pdbUrl")
                if pdb_url:
                    req_pdb = urllib.request.Request(pdb_url, headers={'User-Agent': 'Mozilla/5.0'})
                    with urllib.request.urlopen(req_pdb, timeout=10) as resp_pdb:
                        return resp_pdb.read().decode('utf-8')
    except Exception:
        pass

    urls_fallback = [
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_clean}-F1-model_v4.pdb",
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_clean}-F1-model_v3.pdb",
        f"https://files.rcsb.org/download/{uniprot_clean}.pdb"
    ]

    for url in urls_fallback:
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=8) as resp:
                conteudo = resp.read().decode('utf-8')
                if "ATOM" in conteudo:
                    return conteudo
        except Exception:
            continue

    return None


def extrair_coordenadas_ca(pdb_text):
    """ Extrai resíduos e coordenadas 3D de C-Alpha do arquivo PDB """
    ca_atoms = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_num = int(line[22:26].strip())
            res_name = line[17:20].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            ca_atoms.append({"res_num": res_num, "res_name": res_name, "x": x, "y": y, "z": z})
    return pd.DataFrame(ca_atoms)


def calcular_rede_ligacoes(df_ca, distancia_corte=8.0):
    """ Calcula rede de contatos/ligações entre resíduos baseada na distância espacial """
    if df_ca.empty:
        return nx.Graph(), 0, pd.DataFrame()

    coords = df_ca[['x', 'y', 'z']].values
    n = len(df_ca)

    G = nx.Graph()
    for _, row in df_ca.iterrows():
        G.add_node(row['res_num'], name=f"{row['res_name']}{row['res_num']}")

    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff ** 2, axis=-1))

    for i in range(n):
        for j in range(i + 2, n):
            d = dist_matrix[i, j]
            if d <= distancia_corte:
                r1, r2 = df_ca.iloc[i]['res_num'], df_ca.iloc[j]['res_num']
                G.add_edge(r1, r2, weight=d)

    grados = dict(G.degree())
    df_ligacoes = pd.DataFrame([
        {"Posição (AA)": node, "Resíduo": G.nodes[node]['name'], "Número de Ligações (3D)": deg}
        for node, deg in grados.items()
    ]).sort_values(by="Número de Ligações (3D)", ascending=False)

    return G, G.number_of_edges(), df_ligacoes


def limpar_e_converter(serie):
    s = serie.astype(str).str.replace(',', '.').str.strip()
    s = s.replace(['nan', 'none', 'n/a', '-', 'null', '', '<na>', 'None', 'NONE'], np.nan)
    return pd.to_numeric(s, errors='coerce')


def extrair_posicao_de_mudanca(serie):
    def _extrai_pos(x):
        if not isinstance(x, str) or x.lower() in ['none', 'nan', '', 'null', 'n/a']:
            return np.nan
        m = re.search(r"\d+", x)
        return float(m.group(0)) if m else np.nan
    return serie.apply(_extrai_pos)


def extrair_detalhes_proteina(serie):
    def _extrai(x):
        if not isinstance(x, str) or x.lower() in ['none', 'nan', '', 'null', 'n/a']:
            return pd.Series([np.nan, np.nan, np.nan])

        m = re.search(r"(?:p\.)?([A-Za-z]{1,3})?(\d+)([A-Za-z]{1,3})?", x)
        if m:
            aa_ref = m.group(1) if m.group(1) else np.nan
            pos = float(m.group(2))
            aa_alt = m.group(3) if m.group(3) else np.nan
            return pd.Series([pos, aa_ref, aa_alt])
        return pd.Series([np.nan, np.nan, np.nan])

    res = serie.apply(_extrai)
    res.columns = ['Posicao_Proteina', 'AA_Referencia', 'AA_Alterado']
    return res


def montar_bloco_fasta(sequencia, pos_destacada=None, largura=60):
    """Monta a sequência em blocos numerados, destacando a posição da mutação."""
    linhas = []
    for inicio in range(0, len(sequencia), largura):
        trecho = sequencia[inicio:inicio + largura]
        numero = f'<span class="fasta-pos">{str(inicio + 1).rjust(6)}  </span>'
        if pos_destacada is not None and inicio <= pos_destacada < inicio + largura:
            rel = pos_destacada - inicio
            trecho = (
                trecho[:rel]
                + f'<span class="fasta-hit">{trecho[rel]}</span>'
                + trecho[rel + 1:]
            )
        linhas.append(numero + trecho)
    return '<div class="fasta-box">' + "\n".join(linhas) + "</div>"


# Inicialização
df = carregar_dados_otimizado()

st.markdown("""
<div class="app-header">
    <div>
        <p class="app-name">Painel de Análise Genômica e Estrutural</p>
        <p class="app-desc">Consolida variantes de dbSNP e AlphaMissense por gene e cruza com anotações do UniProt e modelos 3D do AlphaFold.</p>
    </div>
    <div class="app-sources">dbSNP  |  AlphaMissense  |  UniProt  |  AlphaFold DB</div>
</div>
""", unsafe_allow_html=True)

if df.empty:
    st.error("Nenhum arquivo de dados consolidado foi encontrado. Coloque o arquivo .tsv na pasta do projeto ou em /resultados e recarregue a página.")
else:
    col_gene = next((c for c in df.columns if c.lower() in ['gene', 'gene_name', 'symbol']), None)

    if not col_gene:
        st.error("Coluna de identificação do Gene não foi encontrada nos dados. Verifique se existe uma coluna chamada Gene, Gene_Name ou Symbol.")
    else:
        genes_disponiveis = sorted(df[col_gene].dropna().astype(str).str.strip().str.upper().unique())

        st.markdown('<p class="gene-search-label">Buscar gene</p>', unsafe_allow_html=True)
        gene_selecionado = st.selectbox(
            "Buscar gene",
            options=[""] + genes_disponiveis,
            index=0,
            placeholder="Digite para buscar (Ex: ABCB11, OLR1, TP53...)",
            label_visibility="collapsed"
        )

        if gene_selecionado:
            df_gene = df[df[col_gene].astype(str).str.strip().str.upper() == gene_selecionado].copy()

            if not df_gene.empty:
                col_uniprot = next((c for c in df_gene.columns if any(k in c.lower() for k in ['uniprot_id', 'id_uniprot', 'uniprot', 'accession'])), None)
                col_nome_prot = next((c for c in df_gene.columns if any(k in c.lower() for k in ['nome_proteina', 'protein_name', 'proteina'])), None)
                col_rs = next((c for c in df_gene.columns if 'rs' in c.lower() or 'dbsnp' in c.lower()), None)
                col_score = next((c for c in df_gene.columns if any(k in c.lower() for k in ['alphamissense_score', 'am_pathogenicity', 'score', 'pathogenicity'])), None)
                col_class_dbsnp = next((c for c in df_gene.columns if any(k in c.lower() for k in ["classificacao", "clinical", "significancia", "clinvar"])), None)
                col_pos = next((c for c in df_gene.columns if any(k in c.lower() for k in ['posicao_proteina', 'pos_proteina', 'protein_pos', 'aa_pos'])), None)
                col_change = next((c for c in df_gene.columns if any(k in c.lower() for k in ["proteina_mudanca", "protein_change", "aa_change"])), None)

                uniprot_val = "N/A"
                if col_uniprot:
                    valid_uniprots = df_gene[col_uniprot].dropna().astype(str).str.strip()
                    valid_uniprots = valid_uniprots[~valid_uniprots.str.lower().isin(['none', 'nan', '', 'null', 'n/a'])]
                    if not valid_uniprots.empty:
                        uniprot_val = valid_uniprots.iloc[0]

                nome_prot_val = "N/A"
                if col_nome_prot:
                    valid_prots = df_gene[col_nome_prot].dropna().astype(str).str.strip()
                    valid_prots = valid_prots[~valid_prots.str.lower().isin(['none', 'nan', '', 'null', 'n/a'])]
                    if not valid_prots.empty:
                        nome_prot_val = valid_prots.iloc[0]

                if col_change:
                    detalhes_aa = extrair_detalhes_proteina(df_gene[col_change])
                    df_gene['Posicao_Proteina'] = detalhes_aa['Posicao_Proteina']
                    df_gene['AA_Ref'] = detalhes_aa['AA_Referencia']
                    df_gene['AA_Alt'] = detalhes_aa['AA_Alterado']
                    df_gene['Mudanca_Proteina'] = df_gene[col_change]
                else:
                    df_gene['Posicao_Proteina'] = np.nan
                    df_gene['AA_Ref'] = np.nan
                    df_gene['AA_Alt'] = np.nan
                    df_gene['Mudanca_Proteina'] = np.nan

                col_gene_pos = next((c for c in df_gene.columns if any(k in c.lower() for k in ['posicao_gene', 'gene_pos', 'cds_pos', 'pos_gene', 'hgvs_c', 'cdna'])), None)
                col_genomic_pos = next((c for c in df_gene.columns if any(k in c.lower() for k in ['posicao_genomica', 'genomic_pos', 'chr_pos', 'hgvs_g', 'chrom'])), None)
                col_consequence = next((c for c in df_gene.columns if any(k in c.lower() for k in ['consequencia', 'consequence', 'so_term', 'functional_class'])), None)

                df_gene['Posicao_Gene'] = df_gene[col_gene_pos] if col_gene_pos else np.nan
                df_gene['Posicao_Genomica'] = df_gene[col_genomic_pos] if col_genomic_pos else np.nan
                df_gene['Consequencia_Gene'] = df_gene[col_consequence] if col_consequence else np.nan
                df_gene['Link_AlphaFold'] = f"{URL_BASE_ALPHAFOLD}{uniprot_val}" if uniprot_val != "N/A" else None
                df_gene['AM_Score_Num'] = limpar_e_converter(df_gene[col_score]) if col_score else np.nan

                fasta_original, anotacoes_uniprot = buscar_info_uniprot_api(uniprot_val)

                # --- CARTÃO DE IDENTIFICAÇÃO ---
                badge_html = f'<span class="uniprot-badge">UniProt {uniprot_val}</span>' if uniprot_val != "N/A" else ""
                subtitulo_html = f'<div class="protein-subtitle"><b>Proteína</b> — {nome_prot_val}</div>' if nome_prot_val != "N/A" else ""

                st.markdown(f"""
                <div class="header-card">
                    <div class="gene-label">Gene selecionado</div>
                    <div class="gene-title-container">
                        <h1 class="gene-title">{gene_selecionado}</h1>
                        {badge_html}
                    </div>
                    {subtitulo_html}
                </div>
                """, unsafe_allow_html=True)

                # --- MÉTRICAS ---
                c1, c2, c3, c4 = st.columns(4)
                with c1:
                    st.metric("Total de registros", f"{len(df_gene):,}".replace(",", "."))
                with c2:
                    if col_rs:
                        rs_vals = df_gene[col_rs].dropna().astype(str).str.strip()
                        rs_count = rs_vals[~rs_vals.isin(['-', 'nan', '', 'N/A', 'None', '<NA>'])].nunique()
                        st.metric("Variantes dbSNP (rsID)", f"{rs_count:,}".replace(",", "."))
                    else:
                        st.metric("Variantes dbSNP", "0")
                with c3:
                    trocas = df_gene['Mudanca_Proteina'].dropna().astype(str).str.strip()
                    trocas_validas = trocas[~trocas.isin(['-', 'nan', '', 'N/A', 'None', '<NA>', 'null'])].nunique()
                    st.metric("Trocas únicas (AA)", f"{trocas_validas:,}".replace(",", "."))
                with c4:
                    vals_s = df_gene['AM_Score_Num'].dropna()
                    m_score = vals_s.mean() if not vals_s.empty else np.nan
                    st.metric("Score médio AlphaMissense", f"{m_score:.3f}" if not np.isnan(m_score) else "N/A")

                # --- PAINEL PLOTLY ---
                st.markdown('<div class="section-header">Painel visual da proteína</div>', unsafe_allow_html=True)
                df_valid = df_gene.dropna(subset=['AM_Score_Num']).copy()

                if not df_valid.empty:
                    fig = make_subplots(
                        rows=2, cols=2,
                        specs=[[{"colspan": 2}, None], [{}, {}]],
                        subplot_titles=(
                            "Patogenicidade por posição na proteína",
                            "Distribuição dos scores",
                            "Scores por classificação dbSNP / ClinVar"
                        ),
                        vertical_spacing=0.18,
                        horizontal_spacing=0.1
                    )

                    eixo_x = df_valid['Posicao_Proteina'] if df_valid['Posicao_Proteina'].notna().any() else np.arange(len(df_valid))

                    fig.add_trace(
                        go.Scatter(
                            x=eixo_x,
                            y=df_valid['AM_Score_Num'],
                            mode='markers',
                            marker=dict(size=7, color=COR_ACENTO, opacity=0.65,
                                        line=dict(width=0.5, color="rgba(226,232,240,0.35)")),
                            name='Variantes',
                            text=df_valid[col_rs] if col_rs else None,
                            customdata=df_valid['Mudanca_Proteina'] if 'Mudanca_Proteina' in df_valid.columns else None,
                            hovertemplate="<b>rsID:</b> %{text}<br><b>Mudança:</b> %{customdata}<br><b>Posição AA:</b> %{x}<br><b>Score:</b> %{y:.3f}<extra></extra>"
                        ),
                        row=1, col=1
                    )

                    fig.add_hline(y=0.56, line_dash="dash", line_color=COR_ALERTA, line_width=1,
                                  annotation_text="Patogênico ≥ 0,56", annotation_position="top left",
                                  annotation_font=dict(size=11, color=COR_ALERTA), row=1, col=1)
                    fig.add_hline(y=0.34, line_dash="dash", line_color=COR_OK, line_width=1,
                                  annotation_text="Benigno ≤ 0,34", annotation_position="bottom left",
                                  annotation_font=dict(size=11, color=COR_OK), row=1, col=1)

                    fig.add_trace(
                        go.Histogram(x=df_valid['AM_Score_Num'], nbinsx=15,
                                     marker=dict(color=COR_ACENTO, opacity=0.75,
                                                 line=dict(width=1, color=COR_FUNDO)),
                                     name='Frequência'),
                        row=2, col=1
                    )

                    if col_class_dbsnp and col_class_dbsnp in df_valid.columns:
                        df_box = df_valid.dropna(subset=[col_class_dbsnp])
                        for cat in df_box[col_class_dbsnp].unique():
                            sub_df = df_box[df_box[col_class_dbsnp] == cat]
                            fig.add_trace(
                                go.Box(y=sub_df['AM_Score_Num'], name=str(cat), boxpoints='all',
                                       jitter=0.3, pointpos=-1.8,
                                       marker=dict(size=4, color=COR_ACENTO, opacity=0.6),
                                       line=dict(color="#94A3B8", width=1.2),
                                       fillcolor="rgba(56,189,248,0.08)"),
                                row=2, col=2
                            )

                    fig.update_layout(title="", height=640, showlegend=False, margin=dict(l=20, r=20, t=60, b=40))
                    aplicar_tema(fig)
                    for anot in fig.layout.annotations:
                        anot.font.size = 13
                        anot.font.color = "#CBD5E1"
                    fig.update_xaxes(title_text="Posição do aminoácido", row=1, col=1)
                    fig.update_yaxes(title_text="AlphaMissense score", row=1, col=1)
                    fig.update_xaxes(title_text="AlphaMissense score", row=2, col=1)
                    fig.update_yaxes(title_text="Frequência", row=2, col=1)
                    fig.update_xaxes(title_text="Classificação", tickangle=-25, tickfont=dict(size=10), row=2, col=2)
                    fig.update_yaxes(title_text="AlphaMissense score", row=2, col=2)

                    event = st.plotly_chart(fig, use_container_width=True, on_select="rerun")

                    if event and "select" in event and event["select"]["points"]:
                        pt = event["select"]["points"][0]
                        pos_clicada = pt.get("x")

                        if pos_clicada is not None:
                            df_detalhe = df_valid[df_valid['Posicao_Proteina'] == pos_clicada]

                            if not df_detalhe.empty:
                                st.info(f"Sítio selecionado no gráfico: posição {int(pos_clicada)}")
                                cols_det = [col_rs, 'Mudanca_Proteina', 'Posicao_Proteina', col_score, col_class_dbsnp]
                                cols_existentes = [c for c in cols_det if c in df_detalhe.columns]

                                st.dataframe(
                                    df_detalhe[cols_existentes],
                                    use_container_width=True,
                                    hide_index=True,
                                    column_config={
                                        col_score: st.column_config.NumberColumn("AlphaMissense score", format="%.3f"),
                                        "Posicao_Proteina": st.column_config.NumberColumn("Posição (AA)", format="%d")
                                    }
                                )

                # --- MUTAÇÕES POR SÍTIO ---
                if col_pos and df_gene[col_pos].notna().any():
                    posicoes = limpar_e_converter(df_gene[col_pos])
                elif 'Posicao_Proteina' in df_gene.columns and df_gene['Posicao_Proteina'].notna().any():
                    posicoes = df_gene['Posicao_Proteina']
                elif col_change and col_change in df_gene.columns:
                    posicoes = extrair_posicao_de_mudanca(df_gene[col_change])
                else:
                    posicoes = pd.Series(dtype=float)

                posicoes_validas = posicoes.dropna()
                st.markdown(f'<div class="section-header">Mutações por sítio ({len(posicoes_validas)})</div>', unsafe_allow_html=True)

                if not posicoes_validas.empty:
                    contagem_sitio = posicoes_validas.astype(int).value_counts().sort_index()
                    hotspots = contagem_sitio[contagem_sitio > 1].sort_values(ascending=False)

                    fig_sitio = go.Figure()

                    for pos, qtd in contagem_sitio.items():
                        cor_haste = COR_ALERTA if qtd > 1 else '#475569'
                        fig_sitio.add_trace(go.Scatter(
                            x=[pos, pos],
                            y=[0, qtd],
                            mode='lines',
                            line=dict(color=cor_haste, width=1.5 if qtd == 1 else 2.5),
                            hoverinfo='skip',
                            showlegend=False
                        ))

                    fig_sitio.add_trace(go.Scatter(
                        x=contagem_sitio.index,
                        y=contagem_sitio.values,
                        mode='markers',
                        marker=dict(
                            size=[10 if y > 1 else 5 for y in contagem_sitio.values],
                            color=[COR_ALERTA if y > 1 else '#64748B' for y in contagem_sitio.values],
                            line=dict(width=1, color=COR_FUNDO)
                        ),
                        hovertemplate="<b>Sítio (posição):</b> %{x}<br><b>Mutações:</b> %{y}<extra></extra>",
                        showlegend=False
                    ))

                    fig_sitio.update_layout(
                        title="",
                        height=400,
                        margin=dict(l=20, r=20, t=30, b=20),
                        xaxis=dict(title="Posição na proteína (sítio)", showgrid=False, zeroline=True),
                        yaxis=dict(title="Número de mutações", dtick=1, zeroline=True)
                    )
                    aplicar_tema(fig_sitio)

                    st.caption(
                        f"{len(posicoes_validas)} variantes distribuídas em {contagem_sitio.shape[0]} sítios. "
                        f"Em vermelho, os {len(hotspots)} sítios com mais de uma mutação."
                    )
                    st.plotly_chart(fig_sitio, use_container_width=True)

                else:
                    st.info("Nenhuma posição válida encontrada para este gene.")

                # --- ESTRUTURA POR ABAS ---
                st.markdown('<div class="section-header">Detalhamento das alterações</div>', unsafe_allow_html=True)
                tab_proteina, tab_gene, tab_resumo, tab_3d = st.tabs([
                    "Proteína",
                    "Gene, genômica e FASTA",
                    "Atributos e anotações",
                    "Estrutura 3D e rede de contatos"
                ])

                colunas_redundantes = ['AlphaMissense_Classificacao', 'AM_Score_Num']

                with tab_proteina:
                    st.markdown('<div class="sub-header">Mutações com foco na estrutura proteica e nos aminoácidos</div>', unsafe_allow_html=True)
                    cols_preferenciais = [col_rs, 'Mudanca_Proteina', 'Posicao_Proteina', 'AA_Ref', 'AA_Alt', col_score, col_class_dbsnp, 'Link_AlphaFold']
                    cols_finais_prot = [c for c in cols_preferenciais if c and c in df_gene.columns]
                    outras = [c for c in df_gene.columns if c not in cols_finais_prot and c not in colunas_redundantes + ['Posicao_Gene', 'Posicao_Genomica', 'Consequencia_Gene']]

                    st.dataframe(
                        df_gene[cols_finais_prot + outras],
                        use_container_width=True,
                        hide_index=True,
                        column_config={
                            col_score: st.column_config.NumberColumn("AlphaMissense score", format="%.3f"),
                            col_class_dbsnp: st.column_config.TextColumn("Significância clínica"),
                            "Mudanca_Proteina": st.column_config.TextColumn("Mudança na proteína (HGVS p.)"),
                            "Posicao_Proteina": st.column_config.NumberColumn("Posição (AA)", format="%d"),
                            "AA_Ref": st.column_config.TextColumn("AA ref"),
                            "AA_Alt": st.column_config.TextColumn("AA alt"),
                            "Link_AlphaFold": st.column_config.LinkColumn("AlphaFold 3D", display_text="Ver modelo 3D")
                        }
                    )

                with tab_gene:
                    st.markdown('<div class="sub-header">Posição no gene (cDNA/CDS), coordenada genômica e consequência</div>', unsafe_allow_html=True)
                    cols_pref_gene = [col_rs, 'Consequencia_Gene', 'Mudanca_Proteina', col_score, col_class_dbsnp, 'Posicao_Gene', 'Posicao_Genomica']
                    cols_finais_gene = [c for c in cols_pref_gene if c and c in df_gene.columns]
                    outras_gene = [c for c in df_gene.columns if c not in cols_finais_gene and c not in colunas_redundantes + ['Link_AlphaFold', 'Posicao_Proteina', 'AA_Ref', 'AA_Alt']]

                    st.dataframe(
                        df_gene[cols_finais_gene + outras_gene],
                        use_container_width=True,
                        hide_index=True,
                        column_config={
                            col_score: st.column_config.NumberColumn("AlphaMissense score", format="%.3f"),
                            col_class_dbsnp: st.column_config.TextColumn("Significância clínica"),
                            "Posicao_Gene": st.column_config.TextColumn("Posição no gene / cDNA"),
                            "Posicao_Genomica": st.column_config.TextColumn("Posição genômica / Chr"),
                            "Consequencia_Gene": st.column_config.TextColumn("Consequência funcional"),
                            "Mudanca_Proteina": st.column_config.TextColumn("Mudança na proteína")
                        }
                    )

                    st.markdown('<div class="section-header">Sequência FASTA com a mutação destacada</div>', unsafe_allow_html=True)

                    if fasta_original:
                        linhas_fasta = fasta_original.strip().split("\n")
                        seq_limpa = "".join(linhas_fasta[1:]).replace("\n", "").replace(" ", "")
                        df_mutacoes = df_gene.dropna(subset=['Posicao_Proteina', 'AA_Ref', 'AA_Alt']).copy()

                        if not df_mutacoes.empty:
                            opcoes_mut = [
                                f"{row[col_rs]} | {row['Mudanca_Proteina']} (Posição {int(row['Posicao_Proteina'])})"
                                for _, row in df_mutacoes.iterrows()
                            ]

                            mut_selecionada = st.selectbox("Escolha a mutação para destacar no FASTA", opcoes_mut)
                            idx_sel = opcoes_mut.index(mut_selecionada)
                            row_sel = df_mutacoes.iloc[idx_sel]

                            pos = int(row_sel['Posicao_Proteina']) - 1
                            aa_alt = str(row_sel['AA_Alt'])

                            if 0 <= pos < len(seq_limpa):
                                aa_real_fasta = seq_limpa[pos]
                                st.info(f"Posição {pos + 1} — referência {aa_real_fasta} para {aa_alt}. Sequência com {len(seq_limpa)} resíduos.")
                                st.markdown(montar_bloco_fasta(seq_limpa, pos_destacada=pos), unsafe_allow_html=True)
                            else:
                                st.warning(f"A posição {pos + 1} está fora da sequência canônica ({len(seq_limpa)} resíduos). Verifique a isoforma usada no mapeamento.")
                                st.markdown(montar_bloco_fasta(seq_limpa), unsafe_allow_html=True)
                        else:
                            st.caption("Sem mutações mapeadas com posição e aminoácidos para destacar. Sequência canônica abaixo.")
                            st.markdown(montar_bloco_fasta(seq_limpa), unsafe_allow_html=True)
                    else:
                        st.info("Sequência FASTA indisponível no UniProt para este identificador.")

                with tab_resumo:
                    st.markdown('<div class="section-header">Anotações funcionais UniProt</div>', unsafe_allow_html=True)

                    if anotacoes_uniprot:
                        col_a1, col_a2 = st.columns(2, gap="medium")

                        with col_a1:
                            # Card 1: Peptídeo Sinal
                            signal_list = anotacoes_uniprot.get("signal_peptide", [])
                            if signal_list:
                                items_html = "".join([f'<li class="annotation-item">{s}</li>' for s in signal_list])
                                signal_content = f'<ul class="annotation-list">{items_html}</ul>'
                            else:
                                signal_content = '<div class="annotation-empty">Nenhum peptídeo sinal mapeado.</div>'

                            st.markdown(f"""
                            <div class="annotation-card">
                                <div class="annotation-title">Peptídeo sinal</div>
                                {signal_content}
                            </div>
                            """, unsafe_allow_html=True)

                            # Card 2: Domínios e Regiões
                            domains_list = anotacoes_uniprot.get("domains", [])
                            if domains_list:
                                items_html = "".join([f'<li class="annotation-item">{d}</li>' for d in domains_list])
                                domain_content = f'<ul class="annotation-list">{items_html}</ul>'
                            else:
                                domain_content = '<div class="annotation-empty">Nenhum domínio identificado.</div>'

                            st.markdown(f"""
                            <div class="annotation-card">
                                <div class="annotation-title">Domínios e regiões</div>
                                {domain_content}
                            </div>
                            """, unsafe_allow_html=True)

                        with col_a2:
                            # Card 3: Sítios Ativos e de Ligação
                            binding_list = anotacoes_uniprot.get("binding_sites", [])
                            if binding_list:
                                items_html = "".join([f'<li class="annotation-item">{b}</li>' for b in binding_list])
                                binding_content = f'<ul class="annotation-list">{items_html}</ul>'
                            else:
                                binding_content = '<div class="annotation-empty">Nenhum sítio de ligação informado.</div>'

                            st.markdown(f"""
                            <div class="annotation-card">
                                <div class="annotation-title">Sítios ativos e de ligação</div>
                                {binding_content}
                            </div>
                            """, unsafe_allow_html=True)

                            # Card 4: Interações Proteína-Proteína
                            interactions_list = anotacoes_uniprot.get("interactions", [])
                            if interactions_list:
                                items_html = "".join([f'<li class="annotation-item">{inter}</li>' for inter in interactions_list])
                                interaction_content = f'<ul class="annotation-list">{items_html}</ul>'
                            else:
                                interaction_content = '<div class="annotation-empty">Nenhuma interação registrada.</div>'

                            st.markdown(f"""
                            <div class="annotation-card">
                                <div class="annotation-title">Interações proteína-proteína</div>
                                {interaction_content}
                            </div>
                            """, unsafe_allow_html=True)
                    else:
                        st.info("Sem anotações UniProt disponíveis para esta proteína.")

                with tab_3d:
                    if uniprot_val == "N/A":
                        st.info("UniProt ID não informado para carregar os dados de estrutura.")
                    else:
                        pdb_code = obter_pdb_alphafold(uniprot_val)

                        if pdb_code:
                            df_ca = extrair_coordenadas_ca(pdb_code)
                            G, total_edges, df_lig = calcular_rede_ligacoes(df_ca, distancia_corte=8.0)

                            col3d_1, col3d_2 = st.columns([1, 1])

                            with col3d_1:
                                st.markdown('<div class="sub-header">Estrutura 3D — backbone C-alfa (AlphaFold)</div>', unsafe_allow_html=True)

                                sig_range = anotacoes_uniprot.get("signal_range")
                                if sig_range:
                                    st.caption(f"Peptídeo sinal destacado nos resíduos {sig_range[0]}–{sig_range[1]}.")
                                else:
                                    st.caption("Nenhum peptídeo sinal mapeado para esta sequência.")

                                if not df_ca.empty:
                                    fig_3d = go.Figure()

                                    # Conexão contínua
                                    fig_3d.add_trace(go.Scatter3d(
                                        x=df_ca['x'], y=df_ca['y'], z=df_ca['z'],
                                        mode='lines',
                                        line=dict(color='rgba(148,163,184,0.45)', width=2),
                                        hoverinfo='none',
                                        showlegend=False
                                    ))

                                    # Se houver peptídeo sinal
                                    if sig_range:
                                        df_sig = df_ca[(df_ca['res_num'] >= sig_range[0]) & (df_ca['res_num'] <= sig_range[1])]
                                        df_rest = df_ca[(df_ca['res_num'] < sig_range[0]) | (df_ca['res_num'] > sig_range[1])]

                                        if not df_sig.empty:
                                            fig_3d.add_trace(go.Scatter3d(
                                                x=df_sig['x'], y=df_sig['y'], z=df_sig['z'],
                                                mode='markers',
                                                marker=dict(size=4.5, color=COR_ALERTA),
                                                name='Peptídeo sinal',
                                                text=[f"{r['res_name']}{r['res_num']}" for _, r in df_sig.iterrows()],
                                                hovertemplate="<b>Peptídeo sinal</b><br>Resíduo: %{text}<extra></extra>"
                                            ))
                                        if not df_rest.empty:
                                            fig_3d.add_trace(go.Scatter3d(
                                                x=df_rest['x'], y=df_rest['y'], z=df_rest['z'],
                                                mode='markers',
                                                marker=dict(size=3, color=COR_ACENTO),
                                                name='Cadeia proteica',
                                                text=[f"{r['res_name']}{r['res_num']}" for _, r in df_rest.iterrows()],
                                                hovertemplate="<b>Resíduo:</b> %{text}<extra></extra>"
                                            ))
                                    else:
                                        fig_3d.add_trace(go.Scatter3d(
                                            x=df_ca['x'], y=df_ca['y'], z=df_ca['z'],
                                            mode='markers',
                                            marker=dict(size=3, color=COR_ACENTO),
                                            name='Resíduos C-alfa',
                                            text=[f"{r['res_name']}{r['res_num']}" for _, r in df_ca.iterrows()],
                                            hovertemplate="<b>Resíduo:</b> %{text}<extra></extra>"
                                        ))

                                    fig_3d.update_layout(
                                        title="",
                                        height=480,
                                        margin=dict(l=0, r=0, t=0, b=0),
                                        legend=dict(orientation="h", yanchor="bottom", y=0, x=0),
                                        scene=dict(
                                            xaxis=dict(visible=False),
                                            yaxis=dict(visible=False),
                                            zaxis=dict(visible=False),
                                            bgcolor='rgba(0,0,0,0)'
                                        )
                                    )
                                    aplicar_tema(fig_3d)
                                    st.plotly_chart(fig_3d, use_container_width=True)

                            with col3d_2:
                                st.markdown('<div class="sub-header">Rede de contatos espaciais (C-alfa &lt; 8 Å)</div>', unsafe_allow_html=True)

                                mc1, mc2 = st.columns(2)
                                with mc1:
                                    st.metric("Resíduos mapeados", f"{len(df_ca):,}".replace(",", "."))
                                with mc2:
                                    st.metric("Contatos espaciais", f"{total_edges:,}".replace(",", "."))

                                if not df_ca.empty and G.number_of_edges() > 0:
                                    edge_x, edge_y, edge_z = [], [], []
                                    pos_dict = {row['res_num']: (row['x'], row['y'], row['z']) for _, row in df_ca.iterrows()}

                                    for edge in G.edges():
                                        x0, y0, z0 = pos_dict[edge[0]]
                                        x1, y1, z1 = pos_dict[edge[1]]
                                        edge_x.extend([x0, x1, None])
                                        edge_y.extend([y0, y1, None])
                                        edge_z.extend([z0, z1, None])

                                    edge_trace = go.Scatter3d(
                                        x=edge_x, y=edge_y, z=edge_z,
                                        line=dict(width=1, color='rgba(148,163,184,0.18)'),
                                        hoverinfo='none', mode='lines', showlegend=False
                                    )

                                    node_x, node_y, node_z, node_degrees, node_text = [], [], [], [], []
                                    for node in G.nodes():
                                        x, y, z = pos_dict[node]
                                        node_x.append(x)
                                        node_y.append(y)
                                        node_z.append(z)
                                        deg = G.degree(node)
                                        node_degrees.append(deg)
                                        node_text.append(f"{G.nodes[node]['name']} ({deg} ligações 3D)")

                                    node_trace = go.Scatter3d(
                                        x=node_x, y=node_y, z=node_z,
                                        mode='markers',
                                        hoverinfo='text',
                                        text=node_text,
                                        marker=dict(
                                            showscale=True,
                                            colorscale='Plasma',
                                            color=node_degrees,
                                            size=5,
                                            colorbar=dict(thickness=10, title=dict(text='Ligações 3D', side='right'),
                                                          outlinewidth=0, x=1.0)
                                        ),
                                        showlegend=False
                                    )

                                    fig_net = go.Figure(data=[edge_trace, node_trace])
                                    fig_net.update_layout(
                                        title="",
                                        height=480,
                                        margin=dict(l=0, r=0, t=0, b=0),
                                        scene=dict(
                                            xaxis=dict(visible=False),
                                            yaxis=dict(visible=False),
                                            zaxis=dict(visible=False),
                                            bgcolor='rgba(0,0,0,0)'
                                        )
                                    )
                                    aplicar_tema(fig_net)
                                    st.plotly_chart(fig_net, use_container_width=True)

                                st.markdown('<div class="sub-header">Resíduos com maior conectividade (hubs)</div>', unsafe_allow_html=True)
                                if not df_lig.empty:
                                    st.dataframe(
                                        df_lig.head(10),
                                        use_container_width=True,
                                        hide_index=True
                                    )
                        else:
                            st.error("Estrutura não encontrada no AlphaFold DB para este UniProt. Confira o identificador ou tente outro gene.")